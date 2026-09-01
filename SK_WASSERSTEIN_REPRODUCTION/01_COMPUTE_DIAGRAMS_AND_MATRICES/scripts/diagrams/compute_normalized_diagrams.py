#!/usr/bin/env python3


from __future__ import annotations

import argparse
import csv
import gc
import math
import os
import re
import shutil
import sys
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple
from xml.sax.saxutils import escape


def load_ttk_plugin() -> None:
    plugin = os.environ.get("TTK_PLUGIN", "").strip()
    if not plugin:
        return
    path = Path(plugin).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"TTK plugin not found: {path}")
    import paraview.simple as pvs

    try:
        pvs.LoadPlugin(str(path), remote=False, ns=pvs.__dict__)
    except TypeError:
        pvs.LoadPlugin(str(path), remote=False)


load_ttk_plugin()

import paraview.simple as pvs
from paraview import servermanager

try:
    import vtk
except ImportError:
    from vtkmodules import all as vtk  

from discover_data import (
    DiscoveryResult,
    GroupCandidate,
    FileInfo,
    ProcessedSpec,
    VidalSpec,
    discover_inputs,
    normalized as normalize_identifier,
)
from vtk_reader_utils import (
    format_description as vtk_format_description,
    open_paraview_reader,
)


REQUIRED_CELL_ARRAYS = (
    "PairIdentifier",
    "PairType",
    "Persistence",
    "Birth",
    "IsFinite",
)
REQUIRED_POINT_ARRAYS = ("ttkVertexScalarField", "CriticalType")


@dataclass(frozen=True)
class ProcessedConfig:
    directory: str
    label: str
    scalar: str
    expected: int


@dataclass(frozen=True)
class VidalConfig:
    directory: str
    label: str
    filename: str
    dimension: str
    expected: int
    clusters: Tuple[int, ...]


PROCESSED_CONFIGS: Tuple[ProcessedConfig, ...] = (
    ProcessedConfig("2006_earthquake_3D", "Earthquake 3D", "VectorMag", 12),
    ProcessedConfig("2008_ionization_front_2D", "Ionization Front 2D", "density", 16),
    ProcessedConfig("2008_ionization_front_3D", "Ionization Front 3D", "density", 16),
    ProcessedConfig("2014_volcanic_eruptions_2D", "Volcanic Eruptions 2D", "SO2", 12),
    ProcessedConfig("2016_viscous_fingering_3D", "Viscous Fingering 3D", "SplatterValues", 15),
    ProcessedConfig("2017_cloud_processes_2D", "Cloud Processes 2D", "ccb", 12),
    ProcessedConfig("2018_asteroid_impact_3D_clustering", "Asteroid Impact 3D Clustering", "tev", 7),
    ProcessedConfig("2018_asteroid_impact_3D_temporal_subsampling", "Asteroid Impact 3D Temporal Subsampling", "scalar", 20),
)

VIDAL_CONFIGS: Tuple[VidalConfig, ...] = (
    VidalConfig(
        "2004_isabel_3D",
        "Isabel 3D",
        "isabella_velocity_goodEnsemble.vti",
        "3D",
        12,
        tuple([0] * 4 + [1] * 4 + [2] * 4),
    ),
    VidalConfig(
        "starting_vortex",
        "Starting Vortex 2D",
        "startingVortexGoodEnsemble.vti",
        "2D",
        12,
        tuple([0] * 6 + [1] * 6),
    ),
    VidalConfig(
        "sea_surface_height",
        "Sea Surface Height 2D",
        "seaSurfaceHeightGoodEnsemble.vti",
        "2D",
        48,
        tuple([0] * 12 + [1] * 12 + [2] * 12 + [3] * 12),
    ),
    VidalConfig(
        "vortex_street",
        "Vortex Street 2D",
        "vortexStreetGoodEnsemble2.vti",
        "2D",
        45,
        tuple([0] * 9 + [1] * 9 + [2] * 9 + [3] * 9 + [4] * 9),
    ),
)

EXCLUDED_ARRAYS = {
    "vtkGhostType",
    "vtkValidPointMask",
    "Normals",
    "TextureCoordinates",
}


def natural_key(value: Any) -> List[Any]:
    return [
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", str(value))
    ]


def delete_proxy(*proxies: Any) -> None:
    for proxy in proxies:
        if proxy is None:
            continue
        try:
            pvs.Delete(proxy)
        except Exception:
            pass


def new_ttk_pd(input_proxy: Any) -> Any:
    constructor = getattr(pvs, "TTKPersistenceDiagram", None)
    if constructor is not None:
        return constructor(Input=input_proxy)
    sm_proxy = servermanager.ProxyManager().NewProxy(
        "filters", "ttkPersistenceDiagram"
    )
    if sm_proxy is None:
        raise RuntimeError("The ttkPersistenceDiagram filter is not registered.")
    proxy = servermanager._getPyProxy(sm_proxy)
    proxy.Input = input_proxy
    return proxy


def set_scalar_field(proxy: Any, scalar: str) -> None:
    attempts: List[str] = []
    for property_name in ("ScalarField", "ScalarFieldNew"):
        for value in (["POINTS", scalar], scalar):
            try:
                setattr(proxy, property_name, value)
                return
            except Exception as exc:
                attempts.append(f"{property_name}={value!r}: {exc}")
    raise RuntimeError(
        f"Could not select field {scalar!r}. " + " | ".join(attempts)
    )


def configure_persistence_diagram(diagram: Any, scalar: str) -> None:
    set_scalar_field(diagram, scalar)
    settings = (
        ("BackEnd", 2),                  
        ("DMSDimensions", 0),            
        ("ShowInsideDomain", 0),         
        ("ClearDGCache", 1),             
        ("ForceInputOffsetScalarField", 0),
        ("IgnoreBoundary", 0),
    )
    for name, value in settings:
        try:
            setattr(diagram, name, value)
        except Exception:
            pass


def point_array_info(proxy: Any) -> List[Tuple[str, int]]:
    
    proxy.UpdatePipeline()
    info = proxy.GetPointDataInformation()
    arrays: List[Tuple[str, int]] = []
    for index in range(info.GetNumberOfArrays()):
        array = info.GetArray(index)
        if array is None:
            continue
        name = array.GetName()
        if name:
            arrays.append((str(name), int(array.GetNumberOfComponents())))
    return arrays


def open_data_reader(path: Path) -> Any:
    
    return open_paraview_reader(path, pvs)


def resolve_scalar_name(wanted: str, available: Sequence[str]) -> str:
    if wanted in available:
        return wanted
    wanted_normalized = normalize_identifier(wanted)
    exact = [name for name in available if normalize_identifier(name) == wanted_normalized]
    if len(exact) == 1:
        return exact[0]
    close = [
        name
        for name in available
        if wanted_normalized in normalize_identifier(name)
        or normalize_identifier(name) in wanted_normalized
    ]
    if len(close) == 1:
        return close[0]
    raise RuntimeError(
        f"Field {wanted!r} could not be identified unambiguously. "
        f"Available fields: {', '.join(available) or '<none>'}"
    )


def deep_fetch_unstructured(proxy: Any) -> Any:
    proxy.UpdatePipeline()
    fetched = servermanager.Fetch(proxy)
    grid = vtk.vtkUnstructuredGrid.SafeDownCast(fetched)
    if grid is None:
        class_name = fetched.GetClassName() if fetched is not None else "None"
        raise RuntimeError(
            "TTK PersistenceDiagram did not produce a vtkUnstructuredGrid "
            f"(received output: {class_name})."
        )
    result = vtk.vtkUnstructuredGrid()
    result.DeepCopy(grid)
    return result


def replace_field_value(field_data: Any, name: str, value: Any) -> None:
    if field_data.GetAbstractArray(name) is not None:
        field_data.RemoveArray(name)
    if isinstance(value, bool):
        array = vtk.vtkIntArray()
        array.SetNumberOfTuples(1)
        array.SetValue(0, int(value))
    elif isinstance(value, int):
        array = vtk.vtkIntArray()
        array.SetNumberOfTuples(1)
        array.SetValue(0, value)
    elif isinstance(value, float):
        array = vtk.vtkDoubleArray()
        array.SetNumberOfTuples(1)
        array.SetValue(0, value)
    else:
        array = vtk.vtkStringArray()
        array.SetNumberOfTuples(1)
        array.SetValue(0, str(value))
    array.SetName(name)
    field_data.AddArray(array)


def add_metadata(grid: Any, metadata: Mapping[str, Any]) -> None:
    field_data = grid.GetFieldData()
    for name, value in metadata.items():
        replace_field_value(field_data, str(name), value)


def validate_diagram(grid: Any, require_unit_square: bool = False) -> int:
    cell_data = grid.GetCellData()
    point_data = grid.GetPointData()
    missing_cells = [
        name for name in REQUIRED_CELL_ARRAYS if cell_data.GetArray(name) is None
    ]
    missing_points = [
        name for name in REQUIRED_POINT_ARRAYS if point_data.GetArray(name) is None
    ]
    if missing_cells or missing_points:
        raise RuntimeError(
            "Incomplete TTK diagram. "
            f"missing Cell Data={missing_cells}, missing Point Data={missing_points}"
        )

    pair_ids = cell_data.GetArray("PairIdentifier")
    births = cell_data.GetArray("Birth")
    persistence = cell_data.GetArray("Persistence")
    count = 0
    tolerance = 1.0e-9
    for index in range(pair_ids.GetNumberOfTuples()):
        if int(pair_ids.GetTuple1(index)) == -1:
            continue
        birth = float(births.GetTuple1(index))
        pers = float(persistence.GetTuple1(index))
        death = birth + pers
        if not (math.isfinite(birth) and math.isfinite(pers) and math.isfinite(death)):
            raise RuntimeError("A diagram contains a non-finite value.")
        if pers < -tolerance or death < birth - tolerance:
            raise RuntimeError("A diagram contains a negative persistence value.")
        if require_unit_square and (
            birth < -tolerance
            or death < -tolerance
            or birth > 1.0 + tolerance
            or death > 1.0 + tolerance
        ):
            raise RuntimeError(
                f"Point outside [0,1]^2: birth={birth}, mort={death}."
            )
        count += 1
    if count < 1:
        raise RuntimeError("The diagram contains no real persistence pair.")
    return count


def write_grid(grid: Any, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(grid)
    try:
        writer.SetDataModeToAppended()
        writer.EncodeAppendedDataOff()
        writer.SetCompressorTypeToZLib()
    except Exception:
        try:
            writer.SetDataModeToBinary()
        except Exception:
            pass
    if writer.Write() != 1:
        raise IOError(f"Could not write: {path}")


def read_grid(path: Path) -> Any:
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    output = vtk.vtkUnstructuredGrid.SafeDownCast(reader.GetOutput())
    if output is None:
        raise RuntimeError(f"Unreadable VTU file: {path}")
    result = vtk.vtkUnstructuredGrid()
    result.DeepCopy(output)
    return result


def raw_bounds(paths: Sequence[Path]) -> Tuple[float, float, int]:
    minimum = math.inf
    maximum = -math.inf
    total_pairs = 0
    for path in paths:
        grid = read_grid(path)
        validate_diagram(grid)
        cells = grid.GetCellData()
        pair_ids = cells.GetArray("PairIdentifier")
        births = cells.GetArray("Birth")
        persistence = cells.GetArray("Persistence")
        for index in range(pair_ids.GetNumberOfTuples()):
            if int(pair_ids.GetTuple1(index)) == -1:
                continue
            birth = float(births.GetTuple1(index))
            death = birth + float(persistence.GetTuple1(index))
            minimum = min(minimum, birth, death)
            maximum = max(maximum, birth, death)
            total_pairs += 1
    if not (math.isfinite(minimum) and math.isfinite(maximum)):
        raise RuntimeError("The global scalar range is non-finite.")
    if maximum <= minimum:
        raise RuntimeError(
            "The collection has zero scalar range and cannot be normalized."
        )
    return minimum, maximum, total_pairs


def clamp_unit(value: float) -> float:
    return min(1.0, max(0.0, value))


def normalize_grid(grid: Any, minimum: float, maximum: float) -> Any:
    result = vtk.vtkUnstructuredGrid()
    result.DeepCopy(grid)
    span = maximum - minimum
    cells = result.GetCellData()
    pair_ids = cells.GetArray("PairIdentifier")
    raw_births = cells.GetArray("Birth")
    raw_persistence = cells.GetArray("Persistence")

    births = vtk.vtkDoubleArray()
    births.SetName("Birth")
    births.SetNumberOfTuples(pair_ids.GetNumberOfTuples())
    persistence = vtk.vtkDoubleArray()
    persistence.SetName("Persistence")
    persistence.SetNumberOfTuples(pair_ids.GetNumberOfTuples())

    max_local_persistence = 0.0
    diagonal_indices: List[int] = []
    points = result.GetPoints()

    for index in range(pair_ids.GetNumberOfTuples()):
        pair_id = int(pair_ids.GetTuple1(index))
        if pair_id == -1:
            diagonal_indices.append(index)
            continue
        raw_birth = float(raw_births.GetTuple1(index))
        raw_death = raw_birth + float(raw_persistence.GetTuple1(index))
        birth = clamp_unit((raw_birth - minimum) / span)
        death = clamp_unit((raw_death - minimum) / span)
        if birth > death:
            birth, death = death, birth
        pers = max(0.0, death - birth)
        births.SetTuple1(index, birth)
        persistence.SetTuple1(index, pers)
        max_local_persistence = max(max_local_persistence, pers)

        
        
        if points is not None:
            cell = result.GetCell(index)
            if cell is not None and cell.GetNumberOfPoints() >= 2:
                p0 = cell.GetPointId(0)
                p1 = cell.GetPointId(1)
                points.SetPoint(p0, birth, birth, 0.0)
                points.SetPoint(p1, birth, death, 0.0)

    for index in diagonal_indices:
        births.SetTuple1(index, 0.0)
        persistence.SetTuple1(index, 2.0 * max_local_persistence)

    cells.RemoveArray("Birth")
    cells.RemoveArray("Persistence")
    cells.AddArray(births)
    cells.AddArray(persistence)
    if points is not None:
        points.Modified()

    add_metadata(
        result,
        {
            "NormalizedForSKOT": 1,
            "NormalizationScope": "one common range for the whole dataset",
            "NormalizationFormula": "(x-min)/(max-min)",
            "NormalizationMinimumRaw": float(minimum),
            "NormalizationMaximumRaw": float(maximum),
            "NormalizationRangeRaw": float(span),
        },
    )
    validate_diagram(result, require_unit_square=True)
    return result


def extract_field_data(grid: Any) -> Dict[str, Any]:
    result: Dict[str, Any] = {}
    field = grid.GetFieldData()
    for index in range(field.GetNumberOfArrays()):
        array = field.GetAbstractArray(index)
        name = array.GetName()
        if not name or array.GetNumberOfTuples() < 1:
            continue
        try:
            if array.IsA("vtkStringArray"):
                value: Any = array.GetValue(0)
            elif array.IsA("vtkDataArray"):
                components = array.GetNumberOfComponents()
                if components == 1:
                    value = array.GetTuple1(0)
                    if float(value).is_integer():
                        value = int(value)
                else:
                    value = ";".join(str(v) for v in array.GetTuple(0))
            else:
                value = str(array.GetVariantValue(0))
            result[str(name)] = value
        except Exception:
            continue
    return result


def write_csv(rows: Sequence[Mapping[str, Any]], path: Path) -> None:
    if not rows:
        return
    fieldnames: List[str] = []
    preferred = ["FILE", "DiagramIndex", "ClusterID", "TimeStep", "RunName"]
    for name in preferred:
        if any(name in row for row in rows):
            fieldnames.append(name)
    for row in rows:
        for name in row:
            if name not in fieldnames:
                fieldnames.append(name)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(dict(row))


def write_vtm(paths: Sequence[Path], destination: Path) -> None:
    lines = [
        '<?xml version="1.0"?>',
        '<VTKFile type="vtkMultiBlockDataSet" version="1.0" byte_order="LittleEndian">',
        "  <vtkMultiBlockDataSet>",
    ]
    for index, path in enumerate(paths):
        relative = os.path.relpath(path, destination.parent).replace(os.sep, "/")
        lines.append(
            f'    <DataSet index="{index}" name="{escape(path.stem)}" '
            f'file="{escape(relative)}"/>'
        )
    lines.extend(["  </vtkMultiBlockDataSet>", "</VTKFile>", ""])
    destination.write_text("\n".join(lines), encoding="utf-8")


def find_dataset_dir(root: Path, name: str) -> Path:
    matches = [path for path in root.rglob(name) if path.is_dir()]
    if not matches:
        raise FileNotFoundError(f"Directory {name!r} is missing from the official archive.")
    matches.sort(
        key=lambda path: (
            0 if (path / "processed").is_dir() else 1,
            len(path.parts),
            str(path),
        )
    )
    return matches[0]


def find_file(root: Path, filename: str) -> Path:
    matches = [path for path in root.rglob(filename) if path.is_file()]
    if not matches:
        raise FileNotFoundError(f"File {filename!r} is missing from the official archive.")
    matches.sort(key=lambda path: (len(path.parts), str(path)))
    return matches[0]


def compute_grid(input_proxy: Any, scalar: str, metadata: Mapping[str, Any]) -> Any:
    diagram = new_ttk_pd(input_proxy)
    configure_persistence_diagram(diagram, scalar)
    try:
        grid = deep_fetch_unstructured(diagram)
    finally:
        delete_proxy(diagram)
    add_metadata(grid, metadata)
    pair_count = validate_diagram(grid)
    replace_field_value(grid.GetFieldData(), "NumberOfPairs", pair_count)
    return grid


def raw_file_is_valid(path: Path) -> bool:
    if not path.is_file():
        return False
    try:
        validate_diagram(read_grid(path))
        return True
    except Exception:
        try:
            path.unlink()
        except Exception:
            pass
        return False


def completed_dataset(output_root: Path, directory: str, expected: int) -> Optional[Path]:
    
    folder = output_root / directory
    marker = folder / "DONE.ok"
    manifest = folder / "diagrams.vtm"
    diagrams = sorted(folder.glob("*.vtu"), key=lambda path: natural_key(path.name))
    if not (marker.is_file() and manifest.is_file() and len(diagrams) == expected):
        return None
    if any(path.stat().st_size <= 0 for path in diagrams):
        return None
    text = manifest.read_text(encoding="utf-8", errors="replace")
    if any(path.name not in text for path in diagrams):
        return None
    return manifest


def finalize_dataset(
    dataset_name: str,
    label: str,
    raw_paths: Sequence[Path],
    output_root: Path,
) -> Tuple[int, Path]:
    output_dir = output_root / dataset_name
    partial_dir = output_root / f".{dataset_name}.in_progress"
    if partial_dir.exists():
        shutil.rmtree(partial_dir)
    partial_dir.mkdir(parents=True, exist_ok=True)

    minimum, maximum, total_pairs = raw_bounds(raw_paths)
    final_paths: List[Path] = []
    rows: List[Dict[str, Any]] = []

    print(
        f"  Common normalization: min={minimum:.17g}, max={maximum:.17g}, "
        f"{len(raw_paths)} diagrams"
    )
    for index, raw_path in enumerate(raw_paths):
        grid = normalize_grid(read_grid(raw_path), minimum, maximum)
        add_metadata(grid, {"DiagramIndex": index, "Dataset": label})
        filename = raw_path.name.replace("_raw.vtu", ".vtu")
        final_path = partial_dir / filename
        write_grid(grid, final_path)
        row = {"FILE": filename}
        row.update(extract_field_data(grid))
        rows.append(row)
        final_paths.append(final_path)
        
        
        
        
        try:
            raw_path.unlink()
        except FileNotFoundError:
            pass
        del grid
        gc.collect()

    write_vtm(final_paths, partial_dir / "diagrams.vtm")
    write_csv(rows, partial_dir / "data.csv")
    write_csv(
        [
            {
                "Dataset": label,
                "Formula": "(x-min)/(max-min)",
                "MinimumRaw": minimum,
                "MaximumRaw": maximum,
                "RangeRaw": maximum - minimum,
                "NumberOfDiagrams": len(final_paths),
                "NumberOfPairs": total_pairs,
            }
        ],
        partial_dir / "normalization_SKOT.csv",
    )
    (partial_dir / "DONE.ok").write_text(
        "Collection computed and normalized in [0,1]^2.\n",
        encoding="utf-8",
    )

    if output_dir.exists():
        shutil.rmtree(output_dir)
    partial_dir.rename(output_dir)
    return len(final_paths), output_dir / "diagrams.vtm"


def process_processed(
    config: ProcessedConfig,
    source: GroupCandidate,
    output_root: Path,
    work_root: Path,
) -> Tuple[int, Path]:
    manifest = completed_dataset(output_root, config.directory, config.expected)
    if manifest is not None:
        print(f"\n=== {config.label}: already complete; preserving the result ===")
        return config.expected, manifest

    files = list(source.files)
    if len(files) != config.expected:
        raise RuntimeError(
            f"{config.label}: {config.expected} expected fields, {len(files)} found "
            f"in {source.anchor}."
        )
    scalar = resolve_scalar_name(config.scalar, source.common_arrays)

    raw_dir = work_root / "raw" / config.directory
    raw_dir.mkdir(parents=True, exist_ok=True)
    raw_paths: List[Path] = []

    print(f"\n=== {config.label}: {len(files)} diagrams ===")
    print(f"  Data detected automatically in: {source.anchor}")
    print(f"  Scalar field used: {scalar}")
    for index, input_file in enumerate(files):
        safe_stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", input_file.stem)
        raw_path = raw_dir / f"diagram_{index:03d}_{safe_stem}_PD_raw.vtu"
        raw_paths.append(raw_path)
        if raw_file_is_valid(raw_path):
            print(f"  [{index + 1}/{len(files)}] already computed: {input_file.name}")
            continue

        print(f"  [{index + 1}/{len(files)}] computing: {input_file.name}")
        if index == 0:
            print(f"  Actual VTK format detected: {vtk_format_description(input_file)}")
        reader = open_data_reader(input_file)
        try:
            reader.PointArrayStatus = [scalar]
        except Exception:
            pass
        arrays = [name for name, _ in point_array_info(reader)]
        actual_scalar = resolve_scalar_name(scalar, arrays)
        grid = compute_grid(
            reader,
            actual_scalar,
            {
                "DiagramIndex": index,
                "Dataset": config.label,
                "ScalarField": actual_scalar,
                "SourceFile": input_file.name,
                "SourcePathInArchive": str(input_file),
                "PersistenceDiagramBackend": "Discrete Morse Sandwich",
                "PersistenceDimensions": "all",
            },
        )
        write_grid(grid, raw_path)
        del grid
        delete_proxy(reader)
        gc.collect()

    count, manifest = finalize_dataset(
        config.directory, config.label, raw_paths, output_root
    )
    shutil.rmtree(raw_dir, ignore_errors=True)
    return count, manifest


def discover_vidal_members(reader: Any, expected: int, label: str):
    arrays = [
        (name, components)
        for name, components in point_array_info(reader)
        if name not in EXCLUDED_ARRAYS and not name.startswith("vtk")
    ]
    scalars = sorted(
        [name for name, components in arrays if components == 1], key=natural_key
    )
    if len(scalars) == expected:
        return "arrays", scalars

    multi = [(name, components) for name, components in arrays if components == expected]
    if len(multi) == 1:
        return "components", (multi[0][0], tuple(range(expected)))

    
    if len(scalars) > expected:
        likely = [
            name
            for name in scalars
            if name.lower() not in {"mask", "validpointmask", "ghost"}
        ]
        if len(likely) >= expected:
            print(
                f"  WARNING: {len(scalars)} scalar fields detected; "
                f"the first {expected} in natural order will be used."
            )
            return "arrays", likely[:expected]

    description = ", ".join(f"{name} ({components} comp.)" for name, components in arrays)
    raise RuntimeError(
        f"{label}: could not identify {expected} members in the VTI. "
        f"Detected arrays: {description or '<none>'}"
    )


def process_vidal(
    config: VidalConfig,
    source: FileInfo,
    output_root: Path,
    work_root: Path,
) -> Tuple[int, Path]:
    manifest = completed_dataset(output_root, config.directory, config.expected)
    if manifest is not None:
        print(f"\n=== {config.label}: already complete; preserving the result ===")
        return config.expected, manifest

    input_file = source.path
    reader = open_data_reader(input_file)
    mode, members = discover_vidal_members(reader, config.expected, config.label)
    raw_dir = work_root / "raw" / config.directory
    raw_dir.mkdir(parents=True, exist_ok=True)
    raw_paths: List[Path] = []

    print(f"\n=== {config.label}: {config.expected} diagrams ===")
    print(f"  Data detected automatically in: {input_file}")
    for index in range(config.expected):
        raw_path = raw_dir / f"diagram_{index:03d}_PD_raw.vtu"
        raw_paths.append(raw_path)
        if raw_file_is_valid(raw_path):
            print(f"  [{index + 1}/{config.expected}] already computed")
            continue

        calculator = None
        if mode == "arrays":
            scalar = str(members[index])
            try:
                reader.PointArrayStatus = [scalar]
            except Exception:
                pass
            domain = reader
            member_name = scalar
        else:
            array_name = str(members[0])
            component = int(members[1][index])
            try:
                reader.PointArrayStatus = [array_name]
            except Exception:
                pass
            scalar = f"EnsembleMember_{index:03d}"
            calculator = pvs.PythonCalculator(Input=reader)
            calculator.ArrayName = scalar
            calculator.Expression = (
                f"inputs[0].PointData[{array_name!r}][:, {component}]"
            )
            domain = calculator
            member_name = f"{array_name}[{component}]"

        print(f"  [{index + 1}/{config.expected}] computing: {member_name}")
        grid = compute_grid(
            domain,
            scalar,
            {
                "DiagramIndex": index,
                "Dataset": config.label,
                "ClusterID": int(config.clusters[index]),
                "Dim.": config.dimension,
                "Generation": "Simulation",
                "MemberArray": member_name,
                "SourceFile": input_file.name,
                "ScalarField": scalar,
                "PersistenceDiagramBackend": "Discrete Morse Sandwich",
                "PersistenceDimensions": "all",
            },
        )
        write_grid(grid, raw_path)
        del grid
        delete_proxy(calculator)
        gc.collect()

    delete_proxy(reader)
    gc.collect()
    count, manifest = finalize_dataset(
        config.directory, config.label, raw_paths, output_root
    )
    shutil.rmtree(raw_dir, ignore_errors=True)
    return count, manifest


def write_top_level_index(
    output_root: Path,
    manifests: Sequence[Tuple[str, Path, int]],
    failures: Sequence[Tuple[str, str]],
) -> None:
    lines = [
        "NORMALIZED PERSISTENCE-DIAGRAM COLLECTIONS FOR SK-OT",
        "",
        "Each diagrams.vtm file loads all .vtu diagrams in its dataset.",
        "The same normalization is applied to every diagram in a dataset:",
        "    x_normalized = (x - global_minimum) / (global_maximum - global_minimum)",
        "",
    ]
    for label, manifest, count in manifests:
        relative = manifest.relative_to(output_root)
        lines.append(f"- {label}: {count} diagrams -> {relative}")
    if failures:
        lines.extend(["", "FAILURES:"])
        for label, error in failures:
            lines.append(f"- {label}: {error}")
    (output_root / "COLLECTION_INDEX.txt").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def build_discovery(data_root: Path, output_root: Path) -> DiscoveryResult:
    processed_specs = [
        ProcessedSpec(
            key=config.directory,
            label=config.label,
            scalar=config.scalar,
            expected=config.expected,
            dimension=2 if "_2D" in config.directory else 3,
        )
        for config in PROCESSED_CONFIGS
    ]
    vidal_specs = [
        VidalSpec(
            key=config.directory,
            label=config.label,
            filename=config.filename,
            expected=config.expected,
            dimension=2 if config.dimension == "2D" else 3,
        )
        for config in VIDAL_CONFIGS
    ]
    diagnostic = output_root / "ARCHIVE_STRUCTURE_DIAGNOSTIC.txt"
    print("Automatically analyzing the actual archive structure...")
    result = discover_inputs(
        data_root=data_root,
        processed_specs=processed_specs,
        vidal_specs=vidal_specs,
        diagnostic_path=diagnostic,
    )
    print(f"  Analyzed structure: {len(result.all_files)} recognized VTK files.")
    print(f"  Diagnostic saved to: {diagnostic}")
    return result


def main(argv: Sequence[str] = sys.argv[1:]) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--work-root", required=True, type=Path)
    args = parser.parse_args(argv)

    data_root = args.data_root.expanduser().resolve()
    output_root = args.output_root.expanduser().resolve()
    work_root = args.work_root.expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    work_root.mkdir(parents=True, exist_ok=True)

    discovery = build_discovery(data_root, output_root)
    missing = [
        config.label for config in PROCESSED_CONFIGS
        if config.directory not in discovery.processed
    ] + [
        config.label for config in VIDAL_CONFIGS
        if config.directory not in discovery.vidal
    ]
    if missing:
        print("\nERROR: some collections could not be recognized automatically:", file=sys.stderr)
        for label in missing:
            print(f"  - {label}", file=sys.stderr)
        print(
            f"Details of the detected structure are in {discovery.diagnostic}",
            file=sys.stderr,
        )
        return 1

    manifests: List[Tuple[str, Path, int]] = []
    failures: List[Tuple[str, str]] = []

    for config in PROCESSED_CONFIGS:
        try:
            count, manifest = process_processed(
                config, discovery.processed[config.directory], output_root, work_root
            )
            manifests.append((config.label, manifest, count))
        except Exception as exc:
            failures.append((config.label, str(exc)))
            print(f"\nERROR for {config.label}: {exc}", file=sys.stderr)
            traceback.print_exc()

    for config in VIDAL_CONFIGS:
        try:
            count, manifest = process_vidal(
                config, discovery.vidal[config.directory], output_root, work_root
            )
            manifests.append((config.label, manifest, count))
        except Exception as exc:
            failures.append((config.label, str(exc)))
            print(f"\nERROR for {config.label}: {exc}", file=sys.stderr)
            traceback.print_exc()

    write_top_level_index(output_root, manifests, failures)
    expected_collections = len(PROCESSED_CONFIGS) + len(VIDAL_CONFIGS)
    print(
        f"\nSummary: {len(manifests)}/{expected_collections} completed collections."
    )
    if failures:
        print(
            "Some collections failed. Rerunning exactly the same script "
            "will reuse completed collections.",
            file=sys.stderr,
        )
        return 1

    (output_root / "ALL_DONE.txt").write_text(
        "All persistence diagrams were computed and normalized for SK-OT.\n"
        "No distance matrix was computed.\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
