#!/usr/bin/env python3


from __future__ import annotations

import os
import re
from html import unescape
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Set, Tuple

try:
    import vtk
except ImportError:  
    from vtkmodules import all as vtk  

from vtk_reader_utils import detect_vtk_xml_type, make_vtk_reader


EXCLUDED_ARRAYS = {
    "vtkGhostType",
    "vtkValidPointMask",
    "Normals",
    "TextureCoordinates",
    "outputOffsets",
    "outputMonotonyffsets",
    "ttkOffsetScalarField",
}

SUPPORTED_VTK_SUFFIXES = {".vtu", ".vti", ".vtr", ".vts", ".vtp", ".vtk", ".pvtu", ".pvti", ".pvtr", ".pvts", ".pvtp"}


@dataclass(frozen=True)
class ProcessedSpec:
    key: str
    label: str
    scalar: str
    expected: int
    dimension: int


@dataclass(frozen=True)
class VidalSpec:
    key: str
    label: str
    filename: str
    expected: int
    dimension: int


@dataclass(frozen=True)
class FileInfo:
    path: Path
    arrays: Tuple[Tuple[str, int], ...]
    dimension: Optional[int]


@dataclass(frozen=True)
class GroupCandidate:
    files: Tuple[Path, ...]
    anchor: Path
    common_arrays: Tuple[str, ...]
    dimension: Optional[int]
    origin: str


@dataclass(frozen=True)
class DiscoveryResult:
    processed: Mapping[str, GroupCandidate]
    vidal: Mapping[str, FileInfo]
    all_files: Tuple[Path, ...]
    diagnostic: Path


ALIASES: Mapping[str, Tuple[str, ...]] = {
    "2006_earthquake_3D": ("earthquake", "seismic", "quake"),
    "2008_ionization_front_2D": ("ionization", "ionisation", "front", "2d"),
    "2008_ionization_front_3D": ("ionization", "ionisation", "front", "3d"),
    "2014_volcanic_eruptions_2D": ("volcan", "eruption", "so2"),
    "2016_viscous_fingering_3D": ("viscous", "finger", "fingering"),
    "2017_cloud_processes_2D": ("cloud", "process", "ccb"),
    "2018_asteroid_impact_3D_clustering": ("asteroid", "impact", "cluster", "tev"),
    "2018_asteroid_impact_3D_temporal_subsampling": (
        "asteroid",
        "impact",
        "temporal",
        "subsampling",
        "time",
    ),
    "2004_isabel_3D": ("isabel", "isabella", "velocity", "3d"),
    "starting_vortex": ("starting", "vortex"),
    "sea_surface_height": ("sea", "surface", "height", "ssh"),
    "vortex_street": ("vortex", "street"),
}


_ARRAY_CACHE: MutableMapping[Path, Tuple[Tuple[str, int], ...]] = {}
_DIMENSION_CACHE: MutableMapping[Path, Optional[int]] = {}


def natural_key(value: Any) -> List[Any]:
    return [
        int(part) if part.isdigit() else part.casefold()
        for part in re.split(r"(\d+)", str(value))
    ]


def normalized(value: Any) -> str:
    text = str(value).casefold()
    return re.sub(r"[^a-z0-9]+", "", text)


def path_text(path: Path) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(path).casefold())


def dimension_from_path(path: Path) -> Optional[int]:
    text = path_text(path)
    tokens = set(text.split())
    if "2d" in tokens or "2" in tokens and "dimension" in tokens:
        return 2
    if "3d" in tokens or "3" in tokens and "dimension" in tokens:
        return 3
    return None


def walk_files_following_links(root: Path) -> List[Path]:
    
    root = root.expanduser().absolute()
    result: List[Path] = []
    visited_dirs: Set[Tuple[int, int]] = set()
    seen_files: Set[Tuple[int, int]] = set()

    for current, dirs, files in os.walk(root, followlinks=True):
        current_path = Path(current)
        try:
            stat = current_path.stat()
            identity = (int(stat.st_dev), int(stat.st_ino))
        except OSError:
            dirs[:] = []
            continue
        if identity in visited_dirs:
            dirs[:] = []
            continue
        visited_dirs.add(identity)

        dirs[:] = sorted(dirs, key=natural_key)
        for filename in sorted(files, key=natural_key):
            path = current_path / filename
            try:
                if not path.is_file():
                    continue
                stat = path.stat()
                identity_file = (int(stat.st_dev), int(stat.st_ino))
            except OSError:
                continue
            if identity_file in seen_files:
                continue
            seen_files.add(identity_file)
            result.append(path)
    return result


def _reader_for_information(path: Path):
    
    
    return make_vtk_reader(path, vtk)


def _xml_point_array_information(path: Path) -> Tuple[Tuple[str, int], ...]:
    
    if detect_vtk_xml_type(path) is None:
        return tuple()
    try:
        with path.open("rb") as handle:
            data = handle.read(16 * 1024 * 1024)
    except OSError:
        return tuple()
    marker = data.find(b"<AppendedData")
    if marker >= 0:
        data = data[:marker]
    text = data.decode("utf-8", errors="ignore")
    sections = re.findall(
        r"<(?:P?PointData)\b[^>]*>(.*?)</(?:P?PointData)>",
        text,
        flags=re.IGNORECASE | re.DOTALL,
    )
    arrays: Dict[str, int] = {}
    for section in sections:
        tags = re.findall(r"<(?:P?DataArray)\b[^>]*>", section, flags=re.IGNORECASE)
        for tag in tags:
            name_match = re.search(
                r"\bName\s*=\s*[\"']([^\"']+)[\"']",
                tag,
                flags=re.IGNORECASE,
            )
            if name_match is None:
                continue
            component_match = re.search(
                r"\bNumberOfComponents\s*=\s*[\"'](\d+)[\"']",
                tag,
                flags=re.IGNORECASE,
            )
            name = unescape(name_match.group(1))
            components = int(component_match.group(1)) if component_match else 1
            arrays[name] = max(components, arrays.get(name, 0))
    return tuple(sorted(arrays.items(), key=lambda item: natural_key(item[0])))


def array_information(path: Path) -> Tuple[Tuple[str, int], ...]:
    cached = _ARRAY_CACHE.get(path)
    if cached is not None:
        return cached

    xml_arrays = _xml_point_array_information(path)
    if xml_arrays:
        _ARRAY_CACHE[path] = xml_arrays
        return xml_arrays

    reader = _reader_for_information(path)
    reader.UpdateInformation()
    arrays: List[Tuple[str, int]] = []

    
    if hasattr(reader, "GetNumberOfPointArrays"):
        count = int(reader.GetNumberOfPointArrays())
        for index in range(count):
            name = reader.GetPointArrayName(index)
            if not name:
                continue
            components = 1
            try:
                info = reader.GetPointDataArraySelection()
                del info
            except Exception:
                pass
            
            arrays.append((str(name), components))

    
    
    
    
    

    deduplicated: Dict[str, int] = {}
    for name, components in arrays:
        deduplicated[name] = max(int(components), deduplicated.get(name, 0))
    value = tuple(sorted(deduplicated.items(), key=lambda item: natural_key(item[0])))
    _ARRAY_CACHE[path] = value
    return value




def full_array_information(path: Path) -> Tuple[Tuple[str, int], ...]:
    
    reader = _reader_for_information(path)
    reader.Update()
    output = reader.GetOutputDataObject(0)
    arrays: List[Tuple[str, int]] = []
    if output is not None and output.GetPointData() is not None:
        point_data = output.GetPointData()
        for index in range(point_data.GetNumberOfArrays()):
            array = point_data.GetAbstractArray(index)
            name = array.GetName() if array is not None else None
            if name:
                arrays.append((str(name), int(array.GetNumberOfComponents())))
    return tuple(sorted(arrays, key=lambda item: natural_key(item[0])))


def dataset_dimension(path: Path) -> Optional[int]:
    cached = _DIMENSION_CACHE.get(path)
    if path in _DIMENSION_CACHE:
        return cached

    dimension: Optional[int] = None
    try:
        reader = _reader_for_information(path)
        reader.UpdateInformation()
        suffix = path.suffix.casefold()
        if suffix in {".vti", ".pvti"}:
            info = reader.GetOutputInformation(0)
            extent_key = vtk.vtkStreamingDemandDrivenPipeline.WHOLE_EXTENT()
            extent = info.Get(extent_key)
            if extent is not None and len(extent) == 6:
                dimension = sum(
                    1 for axis in range(3) if int(extent[2 * axis + 1]) > int(extent[2 * axis])
                )
        if dimension is None:
            reader.Update()
            output = reader.GetOutputDataObject(0)
            if output is not None and hasattr(output, "GetNumberOfCells"):
                maximum = 0
                number = int(output.GetNumberOfCells())
                
                for index in range(min(number, 512)):
                    cell = output.GetCell(index)
                    if cell is not None:
                        maximum = max(maximum, int(cell.GetCellDimension()))
                if maximum > 0:
                    dimension = maximum
    except Exception:
        dimension = None

    _DIMENSION_CACHE[path] = dimension
    return dimension


def usable_scalar_names(path: Path) -> Tuple[str, ...]:
    names = []
    for name, components in array_information(path):
        if components != 1:
            continue
        if name in EXCLUDED_ARRAYS or name.startswith("vtk"):
            continue
        names.append(name)
    return tuple(sorted(names, key=natural_key))


def _common_scalar_names(files: Sequence[Path]) -> Tuple[str, ...]:
    if not files:
        return tuple()
    indices = sorted({0, len(files) // 2, len(files) - 1})
    sets = [set(usable_scalar_names(files[index])) for index in indices]
    if not sets:
        return tuple()
    common = set.intersection(*sets)
    return tuple(sorted(common, key=natural_key))


def _group_key(files: Iterable[Path]) -> Tuple[str, ...]:
    return tuple(sorted((os.path.realpath(path) for path in files), key=natural_key))


def _filename_prefix(path: Path) -> str:
    stem = path.stem.casefold()
    stem = re.sub(r"(?:[_\-.](?:pd|st|jt|field|data))$", "", stem)
    stem = re.sub(
        r"(?:[_\-.]?(?:run|member|sample|step|time|timestep|realization|realisation)?\d+)+$",
        "",
        stem,
    )
    stem = re.sub(r"[_\-.]+$", "", stem)
    return stem


def discover_group_candidates(vtu_files: Sequence[Path], expected_counts: Set[int], root: Path) -> List[GroupCandidate]:
    possible: Dict[Tuple[str, ...], Tuple[Tuple[Path, ...], Path, str]] = {}

    def register(files: Iterable[Path], anchor: Path, origin: str) -> None:
        unique = tuple(sorted(set(files), key=lambda path: natural_key(path.name)))
        if len(unique) not in expected_counts:
            return
        key = _group_key(unique)
        current = possible.get(key)
        
        if current is None:
            possible[key] = (unique, anchor, origin)
            return
        old_files, old_anchor, old_origin = current
        old_score = (40 if "processed" in path_text(old_anchor) else 0) - len(old_anchor.parts)
        new_score = (40 if "processed" in path_text(anchor) else 0) - len(anchor.parts)
        if new_score > old_score:
            possible[key] = (unique, anchor, origin)

    by_parent: Dict[Path, List[Path]] = {}
    for path in vtu_files:
        by_parent.setdefault(path.parent, []).append(path)
    for parent, files in by_parent.items():
        register(files, parent, "same directory")

        by_prefix: Dict[str, List[Path]] = {}
        for path in files:
            prefix = _filename_prefix(path)
            if prefix:
                by_prefix.setdefault(prefix, []).append(path)
        for prefix, prefix_files in by_prefix.items():
            register(prefix_files, parent, f"prefix {prefix!r}")

    
    ancestor_groups: Dict[Path, List[Path]] = {}
    root_abs = root.absolute()
    for path in vtu_files:
        parent = path.parent
        for _ in range(5):
            if parent == root_abs.parent:
                break
            ancestor_groups.setdefault(parent, []).append(path)
            if parent == root_abs or parent.parent == parent:
                break
            parent = parent.parent
    for ancestor, files in ancestor_groups.items():
        register(files, ancestor, "same ancestor")

    candidates: List[GroupCandidate] = []
    for files, anchor, origin in possible.values():
        common = _common_scalar_names(files)
        dimension = dimension_from_path(anchor)
        if dimension is None and files:
            dimension = dataset_dimension(files[0])
        candidates.append(
            GroupCandidate(
                files=files,
                anchor=anchor,
                common_arrays=common,
                dimension=dimension,
                origin=origin,
            )
        )
    candidates.sort(key=lambda group: (len(group.files), natural_key(str(group.anchor))))
    return candidates


def _alias_score(key: str, value: Path) -> int:
    text = path_text(value)
    score = 0
    for alias in ALIASES.get(key, tuple()):
        token = alias.casefold()
        if token in text:
            score += 28 if len(token) >= 4 else 14
    return score


def _dimension_from_key(key: str) -> Optional[int]:
    if "_2D" in key or key.endswith("2D"):
        return 2
    if "_3D" in key or key.endswith("3D"):
        return 3
    return None


def score_processed(spec: ProcessedSpec, candidate: GroupCandidate) -> int:
    if len(candidate.files) != spec.expected:
        return -10_000

    score = 100
    wanted = normalized(spec.scalar)
    arrays = {normalized(name): name for name in candidate.common_arrays}
    if wanted in arrays:
        score += 650
    else:
        close = [name for key, name in arrays.items() if wanted in key or key in wanted]
        if close:
            score += 390
        else:
            return -10_000

    text = path_text(candidate.anchor)
    if "processed" in text or "curated" in text:
        score += 90
    if "nonprocessed" in normalized(text) or "raw" in text.split():
        score -= 180
    if "merge tree" in text or "mergetree" in normalized(text):
        score -= 300
    if "persistence diagram" in text or "persistencediagram" in normalized(text):
        score -= 300

    score += _alias_score(spec.key, candidate.anchor)
    wanted_dimension = spec.dimension or _dimension_from_key(spec.key)
    if candidate.dimension is not None and wanted_dimension is not None:
        score += 170 if candidate.dimension == wanted_dimension else -450
    dim_token = f"{wanted_dimension}d" if wanted_dimension else ""
    if dim_token and dim_token in text:
        score += 85
    return score


def member_layout(path: Path) -> Tuple[Optional[int], str]:
    arrays = [
        (name, components)
        for name, components in array_information(path)
        if name not in EXCLUDED_ARRAYS and not name.startswith("vtk")
    ]
    scalars = [name for name, components in arrays if components == 1]
    if len(scalars) in {12, 45, 48}:
        return len(scalars), "scalar arrays"
    for name, components in arrays:
        if components in {12, 45, 48}:
            return components, f"components of {name}"

    
    
    try:
        arrays = [
            (name, components)
            for name, components in full_array_information(path)
            if name not in EXCLUDED_ARRAYS and not name.startswith("vtk")
        ]
        scalars = [name for name, components in arrays if components == 1]
        if len(scalars) in {12, 45, 48}:
            return len(scalars), "scalar arrays"
        for name, components in arrays:
            if components in {12, 45, 48}:
                return components, f"components of {name}"
    except Exception:
        pass
    return None, "unknown structure"


def score_vidal(spec: VidalSpec, info: FileInfo) -> int:
    count, _ = member_layout(info.path)
    if count != spec.expected:
        return -10_000
    score = 700
    if normalized(info.path.name) == normalized(spec.filename):
        score += 700
    score += _alias_score(spec.key, info.path)
    text = path_text(info.path)
    if info.dimension is not None:
        score += 220 if info.dimension == spec.dimension else -500
    if f"{spec.dimension}d" in text:
        score += 90
    if "processed" in text:
        score -= 80
    return score


def _assign_unique(
    specs: Sequence[Any],
    candidates: Sequence[Any],
    score_function,
    threshold: int,
) -> Dict[str, Any]:
    viable: Dict[str, List[Tuple[int, int]]] = {}
    for spec in specs:
        rows = []
        for index, candidate in enumerate(candidates):
            score = int(score_function(spec, candidate))
            if score >= threshold:
                rows.append((score, index))
        rows.sort(reverse=True)
        viable[spec.key] = rows

    
    
    ordered_specs = sorted(
        specs,
        key=lambda spec: (
            len(viable[spec.key]),
            -(
                viable[spec.key][0][0] - viable[spec.key][1][0]
                if len(viable[spec.key]) > 1
                else (viable[spec.key][0][0] if viable[spec.key] else -10_000)
            ),
            spec.key,
        ),
    )

    assignment: Dict[str, Any] = {}
    used: Set[int] = set()
    pending = list(ordered_specs)
    progress = True
    while pending and progress:
        progress = False
        next_pending = []
        for spec in pending:
            options = [(score, index) for score, index in viable[spec.key] if index not in used]
            if not options:
                next_pending.append(spec)
                continue
            
            
            if len(options) == 1 or len(pending) == 1:
                _, index = options[0]
                assignment[spec.key] = candidates[index]
                used.add(index)
                progress = True
            else:
                next_pending.append(spec)
        pending = next_pending

    
    for spec in pending:
        options = [(score, index) for score, index in viable[spec.key] if index not in used]
        if options:
            _, index = options[0]
            assignment[spec.key] = candidates[index]
            used.add(index)
    return assignment


def _top_level_entries(root: Path) -> List[str]:
    try:
        return [entry.name for entry in sorted(root.iterdir(), key=lambda p: natural_key(p.name))]
    except OSError:
        return []


def discover_inputs(
    data_root: Path,
    processed_specs: Sequence[ProcessedSpec],
    vidal_specs: Sequence[VidalSpec],
    diagnostic_path: Path,
) -> DiscoveryResult:
    data_root = data_root.expanduser().absolute()
    all_files = walk_files_following_links(data_root)
    vtk_files = tuple(
        path for path in all_files if path.suffix.casefold() in SUPPORTED_VTK_SUFFIXES
    )
    vtu_files = tuple(path for path in vtk_files if path.suffix.casefold() in {".vtu", ".pvtu"})

    expected_counts = {spec.expected for spec in processed_specs}
    group_candidates = discover_group_candidates(vtu_files, expected_counts, data_root)
    processed = _assign_unique(
        processed_specs, group_candidates, score_processed, threshold=500
    )

    
    
    
    used_processed_files = {
        os.path.realpath(path)
        for candidate in processed.values()
        for path in candidate.files
    }
    vidal_candidates: List[FileInfo] = []
    
    
    image_candidates = [
        path for path in vtk_files if path.suffix.casefold() in {".vti", ".pvti"}
    ]
    for path in image_candidates:
        if os.path.realpath(path) in used_processed_files:
            continue
        count, _ = member_layout(path)
        if count not in {spec.expected for spec in vidal_specs}:
            continue
        vidal_candidates.append(
            FileInfo(path=path, arrays=array_information(path), dimension=dataset_dimension(path))
        )
    vidal = _assign_unique(vidal_specs, vidal_candidates, score_vidal, threshold=500)

    diagnostic_path.parent.mkdir(parents=True, exist_ok=True)
    lines: List[str] = [
        "AUTOMATIC DIAGNOSTIC OF THE WASSERSTEINMERGETREESDATA ARCHIVE",
        "",
        f"Analyzed root: {data_root}",
        f"Total number of files: {len(all_files)}",
        f"Recognized VTK files: {len(vtk_files)}",
        f"VTU/PVTU files: {len(vtu_files)}",
        "Top-level entries: " + (", ".join(_top_level_entries(data_root)) or "<none>"),
        "",
        "CANDIDATE FIELD GROUPS",
    ]
    if not group_candidates:
        lines.append("<none>")
    for index, candidate in enumerate(group_candidates):
        lines.append(
            f"[{index}] {candidate.anchor} | {len(candidate.files)} files | "
            f"dimension={candidate.dimension} | arrays={', '.join(candidate.common_arrays) or '<none>'} | "
            f"origin={candidate.origin}"
        )

    lines.extend(["", "PREPROCESSED COLLECTION ASSIGNMENTS"])
    for spec in processed_specs:
        candidate = processed.get(spec.key)
        if candidate is None:
            lines.append(f"- {spec.label}: NOT FOUND")
        else:
            lines.append(
                f"- {spec.label}: {candidate.anchor} ({len(candidate.files)} files, "
                f"field {spec.scalar}, dimension {candidate.dimension})"
            )

    lines.extend(["", "CANDIDATE MULTI-MEMBER FILES"])
    if not vidal_candidates:
        lines.append("<none>")
    for info in vidal_candidates:
        count, layout = member_layout(info.path)
        lines.append(
            f"- {info.path} | members={count} ({layout}) | dimension={info.dimension} | "
            f"arrays={', '.join(name for name, _ in info.arrays)}"
        )

    lines.extend(["", "VIDAL COLLECTION ASSIGNMENTS"])
    for spec in vidal_specs:
        info = vidal.get(spec.key)
        lines.append(
            f"- {spec.label}: {info.path if info is not None else 'NOT FOUND'}"
        )

    missing_processed = [spec.label for spec in processed_specs if spec.key not in processed]
    missing_vidal = [spec.label for spec in vidal_specs if spec.key not in vidal]
    if missing_processed or missing_vidal:
        lines.extend(["", "COLLECTIONS STILL MISSING"])
        for label in missing_processed + missing_vidal:
            lines.append(f"- {label}")

    diagnostic_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    return DiscoveryResult(
        processed=dict(processed),
        vidal=dict(vidal),
        all_files=vtk_files,
        diagnostic=diagnostic_path,
    )
