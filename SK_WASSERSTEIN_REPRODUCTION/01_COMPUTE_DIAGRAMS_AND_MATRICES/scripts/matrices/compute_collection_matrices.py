#!/usr/bin/env pvpython

from __future__ import annotations

import argparse
import gc
import json
import math
import shutil
import statistics
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

import matrix_core as core

try:
    import numpy as np
except ImportError as exc:  
    raise RuntimeError("NumPy is required by the selected pvpython.") from exc


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def final_dir_for(output: Path, method: str, level: int | None) -> Path:
    if method == "SKOT":
        return output / "SKOT" / f"L{int(level):02d}"
    if method == "SK_W2DeltaSk":
        return output / "SK_W2DeltaSk" / f"L{int(level):02d}"
    if method == "W2":
        return output / "W2"
    raise ValueError(method)


def timing_target(method: str, level: int | None, requested: int) -> int:
    if method == "W2" or (level == 30 and method in {"SKOT", "SK_W2DeltaSk"}):
        return max(1, int(requested))
    return 1


def load_timing_runs(folder: Path, info: dict[str, Any]) -> list[float]:
    timing_path = folder / "timing_runs.json"
    values: list[float] = []
    if timing_path.is_file():
        try:
            payload = json.loads(timing_path.read_text(encoding="utf-8"))
            values = [float(v) for v in payload.get("compute_seconds", [])]
        except Exception:
            values = []
    if not values:
        raw = info.get("timing_runs_seconds", [])
        if isinstance(raw, list):
            values = [float(v) for v in raw if math.isfinite(float(v)) and float(v) >= 0]
    if not values and "compute_seconds" in info:
        value = float(info["compute_seconds"])
        if math.isfinite(value) and value >= 0:
            values = [value]
    return values


def write_timing_state(
    folder: Path,
    info: dict[str, Any],
    runs: list[float],
    target: int,
    method: str,
    level: int | None,
) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    payload = {
        "method": method,
        "level": level,
        "target_runs": int(target),
        "completed_runs": len(runs),
        "compute_seconds": runs,
        "updated_at_utc": utc_now(),
        "protocol": (
            "fresh filter proxy for each run; wall-clock time around "
            "UpdatePipeline(); input loading, Fetch(), and output writing excluded"
        ),
    }
    (folder / "timing_runs.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    info = dict(info)
    info["timing_runs_seconds"] = list(runs)
    info["timing_runs_target"] = int(target)
    info["timing_runs_completed"] = len(runs)
    info["timing_protocol"] = payload["protocol"]
    if runs:
        info["compute_seconds_mean"] = float(statistics.fmean(runs))
        info["compute_seconds_median"] = float(statistics.median(runs))
        info["compute_seconds_min"] = float(min(runs))
        info["compute_seconds_max"] = float(max(runs))
        info["compute_seconds_stdev"] = (
            float(statistics.stdev(runs)) if len(runs) >= 2 else 0.0
        )
    (folder / "info.json").write_text(
        json.dumps(info, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def ensure_timing_runs(
    reader: Any,
    folder: Path,
    method: str,
    level: int | None,
    delta_lim: float,
    target: int,
) -> dict[str, Any]:
    info_path = folder / "info.json"
    if not info_path.is_file():
        raise FileNotFoundError(info_path)
    info = json.loads(info_path.read_text(encoding="utf-8"))
    runs = load_timing_runs(folder, info)
    write_timing_state(folder, info, runs, target, method, level)

    while len(runs) < target:
        run_index = len(runs) + 1
        print(
            f"    timing {run_index}/{target} : {info.get('display_name', method)}",
            flush=True,
        )
        proxy = None
        try:
            proxy = core.pvs.TTKPersistenceDiagramDistanceMatrix(Input=reader)
            core.configure_filter(proxy, method, level, delta_lim)
            start = time.perf_counter()
            proxy.UpdatePipeline()
            elapsed = time.perf_counter() - start
            if not math.isfinite(elapsed) or elapsed < 0:
                raise RuntimeError(f"Invalid runtime: {elapsed}")
            runs.append(float(elapsed))
            info = json.loads(info_path.read_text(encoding="utf-8"))
            write_timing_state(folder, info, runs, target, method, level)
            print(f"      {elapsed:.6f} s", flush=True)
        finally:
            if proxy is not None:
                try:
                    core.pvs.Delete(proxy)
                except Exception:
                    pass
            gc.collect()

    return json.loads(info_path.read_text(encoding="utf-8"))


def copy_if_present(source: Path | None, destination: Path) -> None:
    if source is not None and source.is_file():
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--dataset-key", required=True)
    parser.add_argument("--dataset-label", required=True)
    parser.add_argument("--expected", type=int, required=True)
    parser.add_argument("--levels", type=int, nargs="+", required=True)
    parser.add_argument("--delta-lim", type=float, default=0.01)
    parser.add_argument("--timing-runs", type=int, default=1)
    parser.add_argument("--source-data-csv", type=Path, default=None)
    parser.add_argument("--normalization-csv", type=Path, default=None)
    return parser.parse_args(argv)


def main(argv: Sequence[str] = sys.argv[1:]) -> int:
    args = parse_args(argv)
    input_path = args.input.expanduser().resolve()
    output = args.output.expanduser().resolve()
    if not input_path.is_file():
        raise FileNotFoundError(input_path)
    levels = sorted(set(int(v) for v in args.levels))
    if levels != [10, 20, 30, 40]:
        raise ValueError(f"Expected levels: 10 20 30 40; received: {levels}")
    if args.timing_runs < 1:
        raise ValueError("--timing-runs must be >= 1")

    output.mkdir(parents=True, exist_ok=True)
    copy_if_present(args.source_data_csv, output / "source_data.csv")
    copy_if_present(args.normalization_csv, output / "normalization_SKOT.csv")
    (output / "source_manifest_path.txt").write_text(
        str(input_path) + "\n", encoding="utf-8"
    )

    print("\n" + "=" * 78)
    print(f"COLLECTION: {args.dataset_label}")
    print(f"Input      : {input_path}")
    print(f"Output     : {output}")
    print(f"Levels L   : {', '.join(map(str, levels))}")
    print(f"Timing runs: {args.timing_runs} for L=30 and W2")
    print("=" * 78, flush=True)

    start_read = time.perf_counter()
    reader = core.pvs.XMLMultiBlockDataReader(FileName=[str(input_path)])
    reader.UpdatePipeline()
    collection = core.servermanager.Fetch(reader)
    input_read_seconds = time.perf_counter() - start_read

    sample_rows, collection_summary = core.inspect_collection(
        collection, input_path, int(args.expected)
    )
    del collection
    gc.collect()
    labels = [
        str(row.get("vtm_file") or f"sample_{int(row['sample_index']):03d}")
        for row in sample_rows
    ]
    core.write_rows_csv(sample_rows, output / "samples.csv")
    collection_info = {
        "dataset_key": args.dataset_key,
        "dataset_label": args.dataset_label,
        "input_manifest": str(input_path),
        "expected_diagrams": int(args.expected),
        "levels": levels,
        "input_read_seconds": input_read_seconds,
        "collection_summary": collection_summary,
        "timing_runs_requested": int(args.timing_runs),
        "created_at_utc": utc_now(),
    }
    (output / "collection_info.json").write_text(
        json.dumps(collection_info, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(
        f"  {collection_summary['number_of_diagrams']} diagrams, "
        f"{collection_summary['number_of_pairs_total']} total pairs",
        flush=True,
    )

    method_infos: list[dict[str, Any]] = []
    tasks: list[tuple[str, int | None]] = (
        [("SKOT", L) for L in levels]
        + [("SK_W2DeltaSk", L) for L in levels]
        + [("W2", None)]
    )
    for method, level in tasks:
        info = core.run_method(
            reader,
            output,
            method,
            level,
            int(args.expected),
            labels,
            float(args.delta_lim),
        )
        folder = final_dir_for(output, method, level)
        target = timing_target(method, level, int(args.timing_runs))
        info = ensure_timing_runs(
            reader, folder, method, level, float(args.delta_lim), target
        )
        method_infos.append(info)

    (output / "METHODS_COMPLETE.json").write_text(
        json.dumps(method_infos, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    (output / "COLLECTION_COMPLETE.ok").write_text(
        f"All matrices for {args.dataset_label} are complete.\n",
        encoding="utf-8",
    )
    try:
        core.pvs.Delete(reader)
    except Exception:
        pass
    print(f"\nCOLLECTION_COMPLETE={args.dataset_key}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
