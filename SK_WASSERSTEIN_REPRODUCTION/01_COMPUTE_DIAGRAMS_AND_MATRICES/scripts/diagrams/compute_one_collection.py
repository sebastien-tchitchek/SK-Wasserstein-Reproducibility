#!/usr/bin/env python3


from __future__ import annotations

import argparse
import sys
import traceback
from pathlib import Path
from typing import Sequence

from compute_normalized_diagrams import (
    PROCESSED_CONFIGS,
    VIDAL_CONFIGS,
    completed_dataset,
    process_processed,
    process_vidal,
)
from discover_data import ProcessedSpec, VidalSpec, discover_inputs


def main(argv: Sequence[str] = sys.argv[1:]) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--key", required=True)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--work-root", required=True, type=Path)
    args = parser.parse_args(argv)

    data_root = args.data_root.expanduser().resolve()
    output_root = args.output_root.expanduser().resolve()
    work_root = args.work_root.expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    work_root.mkdir(parents=True, exist_ok=True)

    processed = {config.directory: config for config in PROCESSED_CONFIGS}
    vidal = {config.directory: config for config in VIDAL_CONFIGS}
    if args.key not in processed and args.key not in vidal:
        print(f"Unknown collection: {args.key}", file=sys.stderr)
        return 2

    config = processed.get(args.key) or vidal[args.key]
    manifest = completed_dataset(output_root, config.directory, config.expected)
    if manifest is not None:
        print(f"COLLECTION_ALREADY_COMPLETE={config.directory}")
        print(f"MANIFEST={manifest}")
        return 0

    diagnostic_dir = output_root / "_DIAGNOSTICS"
    diagnostic_dir.mkdir(parents=True, exist_ok=True)
    diagnostic = diagnostic_dir / f"{args.key}.txt"

    try:
        if args.key in processed:
            pconfig = processed[args.key]
            spec = ProcessedSpec(
                key=pconfig.directory,
                label=pconfig.label,
                scalar=pconfig.scalar,
                expected=pconfig.expected,
                dimension=2 if "_2D" in pconfig.directory else 3,
            )
            discovery = discover_inputs(
                data_root=data_root,
                processed_specs=[spec],
                vidal_specs=[],
                diagnostic_path=diagnostic,
            )
            source = discovery.processed.get(args.key)
            if source is None:
                raise RuntimeError(
                    f"The {pconfig.expected} fields of {pconfig.label} were not "
                    f"recognized in {data_root}. See {diagnostic}."
                )
            count, manifest = process_processed(
                pconfig, source, output_root, work_root
            )
        else:
            vconfig = vidal[args.key]
            spec = VidalSpec(
                key=vconfig.directory,
                label=vconfig.label,
                filename=vconfig.filename,
                expected=vconfig.expected,
                dimension=2 if vconfig.dimension == "2D" else 3,
            )
            discovery = discover_inputs(
                data_root=data_root,
                processed_specs=[],
                vidal_specs=[spec],
                diagnostic_path=diagnostic,
            )
            source = discovery.vidal.get(args.key)
            if source is None:
                raise RuntimeError(
                    f"The multi-member file for {vconfig.label} was not recognized "
                    f"in {data_root}. See {diagnostic}."
                )
            count, manifest = process_vidal(
                vconfig, source, output_root, work_root
            )

        print(f"COLLECTION_COMPLETE={config.directory}")
        print(f"DIAGRAM_COUNT={count}")
        print(f"MANIFEST={manifest}")
        return 0
    except Exception as exc:
        print(f"COLLECTION_ERROR={config.label}: {exc}", file=sys.stderr)
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
