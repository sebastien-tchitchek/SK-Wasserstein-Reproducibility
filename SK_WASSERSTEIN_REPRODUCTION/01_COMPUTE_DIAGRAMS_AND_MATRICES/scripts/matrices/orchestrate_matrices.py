#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

from collections_config import COLLECTIONS, DEFAULT_LEVELS


def run_with_heartbeat(command: list[str], env: dict[str, str], label: str, heartbeat: int, timeout_minutes: float) -> None:
    print("Commande :", " ".join(command), flush=True)
    process = subprocess.Popen(command, env=env, start_new_session=True)
    start = time.monotonic()
    deadline = start + timeout_minutes * 60.0 if timeout_minutes > 0 else None
    try:
        while True:
            try:
                code = process.wait(timeout=max(1, heartbeat))
                if code != 0:
                    raise RuntimeError(f"{label} stopped with code {code}.")
                return
            except subprocess.TimeoutExpired:
                elapsed = time.monotonic() - start
                print(
                    f"  {label} still active — elapsed: {elapsed/60:.1f} min",
                    flush=True,
                )
                if deadline is not None and time.monotonic() >= deadline:
                    raise TimeoutError(
                        f"{label} exceeded the limit of {timeout_minutes:.1f} minutes."
                    )
    except BaseException:
        if process.poll() is None:
            try:
                os.killpg(process.pid, signal.SIGTERM)
                process.wait(timeout=20)
            except Exception:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except Exception:
                    pass
        raise


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pvpython", type=Path, required=True)
    parser.add_argument("--plugin", default="")
    parser.add_argument("--diagrams-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--compute-script", type=Path, required=True)
    parser.add_argument("--timing-runs", type=int, default=1)
    parser.add_argument("--delta-lim", type=float, default=0.01)
    parser.add_argument("--heartbeat-seconds", type=int, default=120)
    parser.add_argument("--collection-timeout-minutes", type=float, default=600.0)
    args = parser.parse_args()

    pvpython = args.pvpython.expanduser().resolve()
    diagrams_root = args.diagrams_root.expanduser().resolve()
    output_root = args.output_root.expanduser().resolve()
    compute_script = args.compute_script.expanduser().resolve()
    if not pvpython.is_file():
        raise FileNotFoundError(pvpython)
    if not diagrams_root.is_dir():
        raise FileNotFoundError(diagrams_root)
    if not compute_script.is_file():
        raise FileNotFoundError(compute_script)
    output_root.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"
    if args.plugin:
        env["TTK_PLUGIN"] = str(Path(args.plugin).expanduser().resolve())
    else:
        env.pop("TTK_PLUGIN", None)

    for index, cfg in enumerate(COLLECTIONS, start=1):
        key = str(cfg["key"])
        source = diagrams_root / key
        manifest = source / "diagrams.vtm"
        if not manifest.is_file():
            raise FileNotFoundError(manifest)
        print("\n" + "#" * 80, flush=True)
        print(f"COLLECTION {index}/{len(COLLECTIONS)} : {cfg['label']}", flush=True)
        print("#" * 80, flush=True)
        command = [
            str(pvpython),
            str(compute_script),
            "--input", str(manifest),
            "--output", str(output_root / key),
            "--dataset-key", key,
            "--dataset-label", str(cfg["label"]),
            "--expected", str(cfg["expected"]),
            "--levels", *[str(v) for v in DEFAULT_LEVELS],
            "--delta-lim", str(args.delta_lim),
            "--timing-runs", str(args.timing_runs),
        ]
        source_data = source / "data.csv"
        normalization = source / "normalization_SKOT.csv"
        if source_data.is_file():
            command += ["--source-data-csv", str(source_data)]
        if normalization.is_file():
            command += ["--normalization-csv", str(normalization)]
        run_with_heartbeat(
            command,
            env,
            str(cfg["label"]),
            max(10, int(args.heartbeat_seconds)),
            float(args.collection_timeout_minutes),
        )

    print("\nALL 12 COLLECTIONS ARE COMPLETE.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
