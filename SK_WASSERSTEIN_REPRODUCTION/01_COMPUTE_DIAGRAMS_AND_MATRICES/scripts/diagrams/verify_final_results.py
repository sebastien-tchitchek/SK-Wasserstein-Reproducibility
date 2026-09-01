#!/usr/bin/env python3


from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

CONFIGS = (
    ("2006_earthquake_3D", "Earthquake 3D", 12),
    ("2008_ionization_front_2D", "Ionization Front 2D", 16),
    ("2008_ionization_front_3D", "Ionization Front 3D", 16),
    ("2014_volcanic_eruptions_2D", "Volcanic Eruptions 2D", 12),
    ("2016_viscous_fingering_3D", "Viscous Fingering 3D", 15),
    ("2017_cloud_processes_2D", "Cloud Processes 2D", 12),
    ("2018_asteroid_impact_3D_clustering", "Asteroid Impact 3D Clustering", 7),
    ("2018_asteroid_impact_3D_temporal_subsampling", "Asteroid Impact 3D Temporal Subsampling", 20),
    ("2004_isabel_3D", "Isabel 3D", 12),
    ("starting_vortex", "Starting Vortex 2D", 12),
    ("sea_surface_height", "Sea Surface Height 2D", 48),
    ("vortex_street", "Vortex Street 2D", 45),
)


def natural_key(value: str):
    return [int(part) if part.isdigit() else part.casefold() for part in re.split(r"(\d+)", value)]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    root = args.output_root.expanduser().resolve()
    failures: list[str] = []
    lines = [
        "NORMALIZED PERSISTENCE-DIAGRAM COLLECTIONS FOR SK-OT",
        "",
        "Each diagrams.vtm file loads the normalized .vtu diagrams of its dataset.",
        "The same affine transformation is applied to every diagram in a dataset.",
        "No distance matrix was computed.",
        "",
    ]
    total = 0
    for key, label, expected in CONFIGS:
        folder = root / key
        marker = folder / "DONE.ok"
        manifest = folder / "diagrams.vtm"
        diagrams = sorted(folder.glob("*.vtu"), key=lambda p: natural_key(p.name))
        ok = marker.is_file() and manifest.is_file() and len(diagrams) == expected
        if ok:
            text = manifest.read_text(encoding="utf-8", errors="replace")
            ok = all(path.stat().st_size > 0 and path.name in text for path in diagrams)
        if not ok:
            failures.append(f"{label}: incomplete result ({len(diagrams)}/{expected})")
            continue
        total += len(diagrams)
        lines.append(f"- {label}: {expected} diagrams -> {key}/diagrams.vtm")

    if failures or total != 227:
        print("Results are still incomplete:", file=sys.stderr)
        for failure in failures:
            print(f"  - {failure}", file=sys.stderr)
        print(f"Validated total: {total}/227", file=sys.stderr)
        return 1

    (root / "COLLECTION_INDEX.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (root / "ALL_DONE.txt").write_text(
        "All 227 persistence diagrams were computed and normalized for SK-OT.\n"
        "No distance matrix was computed.\n",
        encoding="utf-8",
    )
    print("FINAL_VERIFICATION_OK=1")
    print("COLLECTION_COUNT=12")
    print("DIAGRAM_COUNT=227")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
