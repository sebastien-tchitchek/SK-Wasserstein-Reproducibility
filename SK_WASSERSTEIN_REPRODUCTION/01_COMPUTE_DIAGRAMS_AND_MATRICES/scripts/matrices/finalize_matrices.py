#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import statistics
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from collections_config import (
    COLLECTIONS,
    DEFAULT_LEVELS,
    EXPECTED_PAIRWISE_COMPARISONS,
    EXPECTED_PERSISTENCE_PAIRS,
    EXPECTED_TOTAL_DIAGRAMS,
    EXPECTED_TOTAL_MATRICES,
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def method_folder(root: Path, key: str, method: str, level: int | None) -> Path:
    if method == "SKOT":
        return root / key / "SKOT" / f"L{int(level):02d}"
    if method == "SK_W2DeltaSk":
        return root / key / "SK_W2DeltaSk" / f"L{int(level):02d}"
    return root / key / "W2"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def validate_matrix(path: Path, expected: int) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(path)
    matrix = np.load(path, allow_pickle=False)
    if matrix.shape != (expected, expected):
        raise RuntimeError(f"{path}: shape {matrix.shape}, expected {(expected, expected)}")
    if not np.all(np.isfinite(matrix)):
        raise RuntimeError(f"{path}: NaN or infinity")
    scale = max(1.0, float(np.max(np.abs(matrix))))
    asym = float(np.max(np.abs(matrix - matrix.T)))
    diag = float(np.max(np.abs(np.diag(matrix))))
    minimum = float(np.min(matrix))
    if asym > 1e-9 * scale:
        raise RuntimeError(f"{path}: asymmetry {asym}")
    if diag > 1e-9 * scale:
        raise RuntimeError(f"{path}: diagonal {diag}")
    if minimum < -1e-10 * scale:
        raise RuntimeError(f"{path}: negative minimum {minimum}")
    return {
        "shape": f"{expected}x{expected}",
        "minimum": minimum,
        "maximum": float(np.max(matrix)),
        "max_asymmetry": asym,
        "max_absolute_diagonal": diag,
    }


def read_info(folder: Path) -> dict[str, Any]:
    path = folder / "info.json"
    if not path.is_file():
        raise FileNotFoundError(path)
    return json.loads(path.read_text(encoding="utf-8"))


def timing_values(info: dict[str, Any]) -> list[float]:
    raw = info.get("timing_runs_seconds")
    if isinstance(raw, list) and raw:
        values = [float(v) for v in raw]
    else:
        values = [float(info["compute_seconds"])]
    if not all(math.isfinite(v) and v >= 0 for v in values):
        raise RuntimeError("Invalid timing values")
    return values


def iter_archive_files(root: Path, archive: Path) -> Iterable[Path]:
    for path in sorted(root.rglob("*")):
        if path.is_file() and path.resolve() != archive.resolve():
            yield path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--timing-runs", type=int, default=1)
    parser.add_argument("--delta-lim", type=float, default=0.01)
    parser.add_argument("--pvpython", default="")
    parser.add_argument("--plugin", default="")
    parser.add_argument("--diagrams-root", default="")
    args = parser.parse_args()

    root = args.output_root.expanduser().resolve()
    archive = args.archive.expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(root)

    manifest_rows: list[dict[str, Any]] = []
    timing_rows: list[dict[str, Any]] = []
    total_diagrams = 0
    total_pairs = 0
    matrix_count = 0

    for cfg in COLLECTIONS:
        key = str(cfg["key"])
        expected = int(cfg["expected"])
        collection_dir = root / key
        samples_path = collection_dir / "samples.csv"
        if not samples_path.is_file():
            raise FileNotFoundError(samples_path)
        with samples_path.open(newline="", encoding="utf-8") as handle:
            samples = list(csv.DictReader(handle))
        if len(samples) != expected:
            raise RuntimeError(f"{key}: {len(samples)} samples, expected {expected}")
        counts = [int(float(row["number_of_pairs"])) for row in samples]
        total_diagrams += expected
        total_pairs += sum(counts)

        tasks: list[tuple[str, int | None]] = (
            [("SKOT", L) for L in DEFAULT_LEVELS]
            + [("SK_W2DeltaSk", L) for L in DEFAULT_LEVELS]
            + [("W2", None)]
        )
        for method, level in tasks:
            folder = method_folder(root, key, method, level)
            if not (folder / "DONE.ok").is_file():
                raise FileNotFoundError(folder / "DONE.ok")
            diag = validate_matrix(folder / "matrix_distance.npy", expected)
            info = read_info(folder)
            if info.get("method") != method:
                raise RuntimeError(f"{folder}: inconsistent method")
            expected_level = int(level) if level is not None else None
            if info.get("level") != expected_level:
                raise RuntimeError(f"{folder}: inconsistent level")
            if method == "W2":
                got_delta = float(info.get("parameters", {}).get("delta_lim", math.nan))
                if not math.isclose(got_delta, args.delta_lim, rel_tol=0, abs_tol=1e-15):
                    raise RuntimeError(f"{folder}: DeltaLim={got_delta}")
            if method == "SKOT":
                cost_path = folder / "matrix_SKOT_cost.npy"
                cost = np.load(cost_path, allow_pickle=False)
                distance = np.load(folder / "matrix_distance.npy", allow_pickle=False)
                if not np.allclose(cost, distance * distance, rtol=1e-12, atol=1e-14):
                    raise RuntimeError(f"{cost_path}: does not match d_SK^2")

            runs = timing_values(info)
            expected_runs = (
                max(1, args.timing_runs)
                if method == "W2" or (level == 30 and method in {"SKOT", "SK_W2DeltaSk"})
                else 1
            )
            if len(runs) < expected_runs:
                raise RuntimeError(
                    f"{folder}: {len(runs)}/{expected_runs} timing measurements"
                )
            mean_time = float(statistics.fmean(runs))
            median_time = float(statistics.median(runs))
            stdev_time = float(statistics.stdev(runs)) if len(runs) >= 2 else 0.0
            matrix_count += 1
            common = {
                "dataset_key": key,
                "dataset_label": cfg["label"],
                "paper_label": cfg["paper_label"],
                "number_of_diagrams": expected,
                "method": method,
                "L": "" if level is None else int(level),
                "relative_path": str((folder / "matrix_distance.npy").relative_to(root)),
                **diag,
                "timing_runs": len(runs),
                "compute_seconds_mean": mean_time,
                "compute_seconds_median": median_time,
                "compute_seconds_stdev": stdev_time,
            }
            manifest_rows.append(common)
            for run_index, seconds in enumerate(runs, start=1):
                timing_rows.append(
                    {
                        "dataset_key": key,
                        "dataset_label": cfg["label"],
                        "method": method,
                        "L": "" if level is None else int(level),
                        "run": run_index,
                        "compute_seconds": seconds,
                        "included_in_paper_timing": int(
                            method == "W2"
                            or (level == 30 and method in {"SKOT", "SK_W2DeltaSk"})
                        ),
                    }
                )

    if total_diagrams != EXPECTED_TOTAL_DIAGRAMS:
        raise RuntimeError(f"Total diagrams: {total_diagrams}/{EXPECTED_TOTAL_DIAGRAMS}")
    if total_pairs != EXPECTED_PERSISTENCE_PAIRS:
        raise RuntimeError(f"Total pairs: {total_pairs}/{EXPECTED_PERSISTENCE_PAIRS}")
    if matrix_count != EXPECTED_TOTAL_MATRICES:
        raise RuntimeError(f"Total matrices: {matrix_count}/{EXPECTED_TOTAL_MATRICES}")

    write_csv(root / "GLOBAL_MANIFEST.csv", manifest_rows)
    write_csv(root / "COMPUTE_TIME_RUNS.csv", timing_rows)
    parameters = {
        "benchmark": "SKOT 12 collections",
        "created_at_utc": utc_now(),
        "number_of_collections": len(COLLECTIONS),
        "number_of_diagrams": total_diagrams,
        "number_of_persistence_pairs": total_pairs,
        "number_of_pairwise_comparisons_per_method": EXPECTED_PAIRWISE_COMPARISONS,
        "levels": list(DEFAULT_LEVELS),
        "primary_level": 30,
        "W2": {
            "implementation": "TTK PersistenceDiagramDistanceMatrix auction approximation",
            "p": 2,
            "DeltaLim": args.delta_lim,
            "critical_pair_types": "all available types, computed separately by TTK and combined",
            "constraint": "complete diagrams",
        },
        "timing": {
            "runs_for_primary_methods": int(args.timing_runs),
            "aggregation_for_paper": (
                "single recorded run"
                if int(args.timing_runs) == 1
                else "arithmetic mean"
            ),
            "scope": "wall-clock time around filter UpdatePipeline only",
            "excluded": ["input loading", "Fetch", "output writing"],
        },
        "pvpython": args.pvpython,
        "ttk_plugin": args.plugin,
        "source_diagrams_root": args.diagrams_root,
    }
    (root / "EXPERIMENT_PARAMETERS.json").write_text(
        json.dumps(parameters, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    (root / "README.txt").write_text(
        """SKOT / W_GAMMA / W2 MATRICES — 12-COLLECTION BENCHMARK

This directory contains 108 matrices: d_SK,L and W_Gamma,L for
L=10,20,30,40, together with one complete W2 matrix for each of the
12 collections. Each samples.csv file defines the row and column order.
Primary runtimes are stored in timing_runs.json (one run by default).
""",
        encoding="utf-8",
    )

    
    checksum_lines: list[str] = []
    checksum_path = root / "SHA256SUMS.txt"
    for path in sorted(root.rglob("*")):
        if path.is_file() and path != checksum_path:
            checksum_lines.append(f"{sha256(path)}  {path.relative_to(root)}")
    checksum_path.write_text("\n".join(checksum_lines) + "\n", encoding="utf-8")

    archive.parent.mkdir(parents=True, exist_ok=True)
    if archive.exists():
        archive.unlink()
    with zipfile.ZipFile(
        archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9
    ) as zf:
        prefix = Path(root.name)
        for path in iter_archive_files(root, archive):
            zf.write(path, prefix / path.relative_to(root))

    print("VALIDATION_OK=1")
    print(f"COLLECTIONS={len(COLLECTIONS)}")
    print(f"DIAGRAMS={total_diagrams}")
    print(f"PERSISTENCE_PAIRS={total_pairs}")
    print(f"MATRICES={matrix_count}")
    print(f"ARCHIVE={archive}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
