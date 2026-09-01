#!/usr/bin/env python3


from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Sequence

from prepare_archives import (
    ArchiveInfo,
    archive_kind,
    estimated_uncompressed_size,
    extract_one,
    find_archives,
    human_size,
    safe_target_for,
)

VTK_SUFFIXES = {".vtu", ".pvtu", ".vti", ".pvti", ".vtr", ".pvtr", ".vts", ".pvts", ".vtp", ".pvtp", ".vtk"}


@dataclass(frozen=True)
class Dataset:
    key: str
    label: str
    expected: int


MAIN_DATASETS: tuple[Dataset, ...] = (
    Dataset("2006_earthquake_3D", "Earthquake 3D", 12),
    Dataset("2008_ionization_front_2D", "Ionization Front 2D", 16),
    Dataset("2008_ionization_front_3D", "Ionization Front 3D", 16),
    Dataset("2014_volcanic_eruptions_2D", "Volcanic Eruptions 2D", 12),
    Dataset("2016_viscous_fingering_3D", "Viscous Fingering 3D", 15),
    Dataset("2017_cloud_processes_2D", "Cloud Processes 2D", 12),
    Dataset("2018_asteroid_impact_3D_clustering", "Asteroid Impact 3D Clustering", 7),
    Dataset(
        "2018_asteroid_impact_3D_temporal_subsampling",
        "Asteroid Impact 3D Temporal Subsampling",
        20,
    ),
)

VIDAL_DATASETS: tuple[Dataset, ...] = (
    Dataset("2004_isabel_3D", "Isabel 3D", 12),
    Dataset("starting_vortex", "Starting Vortex 2D", 12),
    Dataset("sea_surface_height", "Sea Surface Height 2D", 48),
    Dataset("vortex_street", "Vortex Street 2D", 45),
)


def natural_key(value: str):
    import re

    return [int(part) if part.isdigit() else part.casefold() for part in re.split(r"(\d+)", value)]


def normalized(value: object) -> str:
    import re

    return re.sub(r"[^a-z0-9]+", "", str(value).casefold())


def iter_files(root: Path) -> Iterable[Path]:
    if not root.exists():
        return
    for current, dirs, files in os.walk(root, followlinks=False):
        dirs[:] = sorted(dirs, key=natural_key)
        base = Path(current)
        for name in sorted(files, key=natural_key):
            path = base / name
            try:
                if path.is_file() and not path.is_symlink():
                    yield path
            except OSError:
                continue


def contains_vtk(root: Path) -> bool:
    return any(path.suffix.casefold() in VTK_SUFFIXES for path in iter_files(root))


def directory_size(root: Path) -> int:
    total = 0
    for path in iter_files(root):
        try:
            total += path.stat().st_size
        except OSError:
            pass
    return total


def disk_free(path: Path) -> int:
    probe = path
    while not probe.exists() and probe.parent != probe:
        probe = probe.parent
    return shutil.disk_usage(probe).free


def result_complete(output_root: Path, dataset: Dataset) -> bool:
    folder = output_root / dataset.key
    marker = folder / "DONE.ok"
    manifest = folder / "diagrams.vtm"
    diagrams = sorted(folder.glob("*.vtu"), key=lambda p: natural_key(p.name))
    if not (marker.is_file() and manifest.is_file() and len(diagrams) == dataset.expected):
        return False
    if any(not path.is_file() or path.stat().st_size <= 0 for path in diagrams):
        return False
    try:
        text = manifest.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return False
    return all(path.name in text for path in diagrams)


def likely_source_directories(data_root: Path, dataset: Dataset) -> list[Path]:
    key_norm = normalized(dataset.key)
    candidates: list[Path] = []
    for current, dirs, _files in os.walk(data_root, topdown=True, followlinks=False):
        base = Path(current)
        for dirname in list(dirs):
            path = base / dirname
            name_norm = normalized(dirname)
            if not (
                name_norm.startswith(key_norm)
                or key_norm in name_norm
                or normalized(str(path.relative_to(data_root))).endswith(key_norm)
            ):
                continue
            if contains_vtk(path):
                candidates.append(path)
    
    
    
    candidates = list(dict.fromkeys(path.resolve() for path in candidates))
    candidates.sort(
        key=lambda path: (
            0 if (path / ".WMT_EXTRACTION_COMPLETE").is_file() else 1,
            len(path.parts),
            directory_size(path),
            str(path),
        )
    )
    
    
    
    roots: list[Path] = []
    for candidate in candidates:
        if any(candidate == root or candidate.is_relative_to(root) for root in roots):
            continue
        roots.append(candidate)
    return roots


def archive_for_dataset(data_root: Path, dataset: Dataset) -> Optional[ArchiveInfo]:
    key_norm = normalized(dataset.key)
    matches = []
    for info in find_archives(data_root):
        name_norm = normalized(info.path.name)
        path_norm = normalized(str(info.path.relative_to(data_root)))
        if key_norm in name_norm or key_norm in path_norm:
            matches.append(info)
    if not matches:
        return None
    matches.sort(key=lambda info: (len(info.path.parts), str(info.path)))
    return matches[0]


def source_for_dataset(data_root: Path, dataset: Dataset) -> Optional[Path]:
    directories = likely_source_directories(data_root, dataset)
    if directories:
        return directories[0]
    archive = archive_for_dataset(data_root, dataset)
    if archive is not None:
        target = safe_target_for(archive.path)
        if (target / ".WMT_EXTRACTION_COMPLETE").is_file() and contains_vtk(target):
            return target
    return None


def safe_remove_source(path: Path, data_root: Path, dataset: Dataset) -> int:
    try:
        resolved = path.resolve()
        root = data_root.resolve()
        resolved.relative_to(root)
    except (OSError, ValueError):
        raise RuntimeError(f"Refusing to delete a path outside the cache: {path}")

    if not (
        normalized(dataset.key) in normalized(resolved.name)
        or (resolved / ".WMT_EXTRACTION_COMPLETE").is_file()
    ):
        raise RuntimeError(f"Refusing to delete an unidentified directory: {resolved}")

    size = directory_size(resolved)
    shutil.rmtree(resolved, ignore_errors=False)
    return size


def remove_archive_if_present(data_root: Path, dataset: Dataset) -> int:
    info = archive_for_dataset(data_root, dataset)
    if info is None:
        return 0
    size = info.size
    info.path.unlink()
    return size


def remove_empty_parents(path: Path, stop: Path) -> None:
    current = path.parent
    stop = stop.resolve()
    while current.exists():
        try:
            current.resolve().relative_to(stop)
        except ValueError:
            break
        if current.resolve() == stop:
            break
        try:
            current.rmdir()
        except OSError:
            break
        current = current.parent


def invoke_collection(
    dataset: Dataset,
    source: Path,
    output_root: Path,
    work_root: Path,
    pvpython: Path,
    plugin: str,
    compute_script: Path,
) -> None:
    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"
    if plugin:
        env["TTK_PLUGIN"] = plugin
    else:
        env.pop("TTK_PLUGIN", None)
    command = [
        str(pvpython),
        str(compute_script),
        "--key",
        dataset.key,
        "--data-root",
        str(source),
        "--output-root",
        str(output_root),
        "--work-root",
        str(work_root),
    ]
    completed = subprocess.run(command, env=env, check=False)
    if completed.returncode != 0:
        raise RuntimeError(
            f"TTK computation for {dataset.label} stopped (code {completed.returncode})."
        )
    if not result_complete(output_root, dataset):
        raise RuntimeError(
            f"Computation for {dataset.label} ended without a complete final collection."
        )


def cleanup_completed_main(data_root: Path, output_root: Path, work_root: Path) -> None:
    for dataset in MAIN_DATASETS:
        if not result_complete(output_root, dataset):
            continue
        for source in likely_source_directories(data_root, dataset):
            if not source.exists():
                continue
            freed = safe_remove_source(source, data_root, dataset)
            print(f"  Resume cleanup : {dataset.label}, {human_size(freed)} freed.")
            remove_empty_parents(source, data_root)
        archived = remove_archive_if_present(data_root, dataset)
        if archived:
            print(f"  Archive no longer needed and deleted: {human_size(archived)} freed.")
        shutil.rmtree(work_root / "raw" / dataset.key, ignore_errors=True)


def choose_next_main(data_root: Path, output_root: Path):
    incomplete = [dataset for dataset in MAIN_DATASETS if not result_complete(output_root, dataset)]
    sources = []
    for dataset in incomplete:
        source = source_for_dataset(data_root, dataset)
        if source is not None:
            sources.append((directory_size(source), dataset, source))
    if sources:
        
        
        sources.sort(key=lambda row: (row[0], row[1].key))
        size, dataset, source = sources[0]
        return "source", dataset, source, size

    archives = []
    for dataset in incomplete:
        info = archive_for_dataset(data_root, dataset)
        if info is None:
            continue
        estimate = estimated_uncompressed_size(info)
        sort_size = estimate if estimate is not None else max(info.size * 4, info.size)
        temporal_penalty = 1 if "temporal_subsampling" in dataset.key else 0
        archives.append((temporal_penalty, sort_size, dataset, info))
    if archives:
        
        
        archives.sort(key=lambda row: (row[0], row[1], row[2].key))
        _, estimate, dataset, info = archives[0]
        return "archive", dataset, info, estimate
    return None


def process_main(
    data_root: Path,
    output_root: Path,
    work_root: Path,
    pvpython: Path,
    plugin: str,
    compute_script: Path,
    outer_archive: Optional[Path],
    report_path: Path,
) -> None:
    output_root.mkdir(parents=True, exist_ok=True)
    work_root.mkdir(parents=True, exist_ok=True)
    data_root.mkdir(parents=True, exist_ok=True)

    
    
    
    has_main_material = any(
        archive_for_dataset(data_root, dataset) is not None
        or source_for_dataset(data_root, dataset) is not None
        for dataset in MAIN_DATASETS
    )
    if outer_archive is not None and outer_archive.is_file() and has_main_material:
        size = outer_archive.stat().st_size
        outer_archive.unlink()
        print(f"Outer archive deleted: {human_size(size)} freed.")

    cleanup_completed_main(data_root, output_root, work_root)
    report_lines = [
        "LOW-DISK PROCESSING — MAIN COLLECTIONS",
        f"Data root: {data_root}",
        f"Results: {output_root}",
        "",
    ]

    while not all(result_complete(output_root, dataset) for dataset in MAIN_DATASETS):
        choice = choose_next_main(data_root, output_root)
        if choice is None:
            missing = [
                dataset.label
                for dataset in MAIN_DATASETS
                if not result_complete(output_root, dataset)
            ]
            raise RuntimeError(
                "No archive or extracted data was found for: "
                + ", ".join(missing)
            )

        kind, dataset, item, size_or_estimate = choice
        free_before = disk_free(data_root)
        print("\n" + "=" * 76)
        print(f"COLLECTION: {dataset.label}")
        print(f"Free space before: {human_size(free_before)}")

        if kind == "archive":
            info: ArchiveInfo = item
            print(
                f"Opening only its archive: {info.path.name} "
                f"({human_size(info.size)} compressed; "
                f"approximately {human_size(size_or_estimate)} extracted)."
            )
            source, count = extract_one(info)
            print(f"Archive opened: {count} file(s) in {source}")
        else:
            source = item
            print(
                "Data already opened by the previous run: "
                f"{source} ({human_size(size_or_estimate)})."
            )

        invoke_collection(
            dataset,
            source,
            output_root,
            work_root,
            pvpython,
            plugin,
            compute_script,
        )

        
        freed_source = safe_remove_source(source, data_root, dataset)
        remove_empty_parents(source, data_root)
        freed_archive = remove_archive_if_present(data_root, dataset)
        shutil.rmtree(work_root / "raw" / dataset.key, ignore_errors=True)
        free_after = disk_free(data_root)
        print(f"Collection validated: {dataset.expected} normalized diagrams.")
        print(
            f"Scalar fields deleted: {human_size(freed_source + freed_archive)} freed."
        )
        print(f"Free space after: {human_size(free_after)}")
        report_lines.append(
            f"- {dataset.label}: complete; {dataset.expected} diagrams; "
            f"{human_size(freed_source + freed_archive)} cleaned; "
            f"libre={human_size(free_after)}"
        )
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    marker = data_root / ".main_collections_complete"
    marker.write_text(
        "The eight main collections were computed and their fields deleted.\n",
        encoding="utf-8",
    )
    print("\nThe eight main collections are complete.")
    print(f"Current free space: {human_size(disk_free(data_root))}")


def find_vidal_file(data_root: Path, dataset: Dataset) -> Optional[Path]:
    exact_names = {
        "2004_isabel_3D": "isabella_velocity_goodEnsemble.vti",
        "starting_vortex": "startingVortexGoodEnsemble.vti",
        "sea_surface_height": "seaSurfaceHeightGoodEnsemble.vti",
        "vortex_street": "vortexStreetGoodEnsemble2.vti",
    }
    wanted = exact_names[dataset.key].casefold()
    matches = [path for path in iter_files(data_root) if path.name.casefold() == wanted]
    matches.sort(key=lambda path: (len(path.parts), str(path)))
    return matches[0] if matches else None


def process_vidal_phase(
    data_root: Path,
    output_root: Path,
    work_root: Path,
    pvpython: Path,
    plugin: str,
    compute_script: Path,
    report_path: Path,
) -> None:
    report_lines = [
        "LOW-DISK PROCESSING — VIDAL COLLECTIONS",
        f"Data root: {data_root}",
        "",
    ]
    
    
    for dataset in VIDAL_DATASETS:
        if result_complete(output_root, dataset):
            print(f"\n=== {dataset.label}: already complete ===")
            report_lines.append(f"- {dataset.label}: already complete")
            continue
        source_file = find_vidal_file(data_root, dataset)
        if source_file is None:
            raise RuntimeError(f"Vidal source file not found for {dataset.label}.")
        print("\n" + "=" * 76)
        print(f"COLLECTION: {dataset.label}")
        print(f"Source : {source_file}")
        invoke_collection(
            dataset,
            data_root,
            output_root,
            work_root,
            pvpython,
            plugin,
            compute_script,
        )
        report_lines.append(f"- {dataset.label}: complete; {dataset.expected} diagrams")
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    print("\nThe four Vidal collections are complete.")


def main(argv: Sequence[str] = sys.argv[1:]) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=("main", "vidal"), required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument("--pvpython", type=Path, required=True)
    parser.add_argument("--plugin", default="")
    parser.add_argument("--compute-script", type=Path, required=True)
    parser.add_argument("--outer-archive", type=Path)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args(argv)

    data_root = args.data_root.expanduser().resolve()
    output_root = args.output_root.expanduser().resolve()
    work_root = args.work_root.expanduser().resolve()
    pvpython = args.pvpython.expanduser().resolve()
    compute_script = args.compute_script.expanduser().resolve()
    report = args.report.expanduser().resolve()
    outer = args.outer_archive.expanduser().resolve() if args.outer_archive else None

    try:
        if args.phase == "main":
            process_main(
                data_root,
                output_root,
                work_root,
                pvpython,
                args.plugin,
                compute_script,
                outer,
                report,
            )
        else:
            process_vidal_phase(
                data_root,
                output_root,
                work_root,
                pvpython,
                args.plugin,
                compute_script,
                report,
            )
        return 0
    except Exception as exc:
        print(f"\nLOW-DISK MODE ERROR: {exc}", file=sys.stderr)
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
