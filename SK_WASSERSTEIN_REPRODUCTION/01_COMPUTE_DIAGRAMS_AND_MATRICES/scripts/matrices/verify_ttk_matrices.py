#!/usr/bin/env pvpython


from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

from ttk_bootstrap import ensure_ttk_plugin_loaded

ensure_ttk_plugin_loaded()

from paraview import servermanager
import paraview.simple as pvs


def prop(proxy, name: str):
    value = proxy.SMProxy.GetProperty(name)
    if value is None:
        raise RuntimeError(
            f"Missing custom property: {name!r}. "
            "The wrong TTK plugin is probably loaded."
        )
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    args = parser.parse_args()
    input_path = args.input.expanduser().resolve()
    if not input_path.is_file():
        raise FileNotFoundError(input_path)

    constructor = getattr(pvs, "TTKPersistenceDiagramDistanceMatrix", None)
    if constructor is None:
        raise RuntimeError("TTKPersistenceDiagramDistanceMatrix is not registered.")

    reader = pvs.XMLMultiBlockDataReader(FileName=[str(input_path)])
    reader.UpdatePipeline()
    collection = servermanager.Fetch(reader)
    if collection is None or collection.GetClassName() != "vtkMultiBlockDataSet":
        got = collection.GetClassName() if collection is not None else "None"
        raise RuntimeError(f"diagrams.vtm did not produce a vtkMultiBlockDataSet ({got}).")
    n = int(collection.GetNumberOfBlocks())
    if n < 2:
        raise RuntimeError(f"The collection contains only {n} diagram(s).")

    matrix_filter = constructor(Input=reader)
    required = (
        "Critical pairs",
        "n",
        "DeltaLim",
        "AntiAlpha",
        "Lambda",
        "Constraint",
        "HilbertInt",
        "ChoiceHilbertDistance",
        "L",
    )
    for name in required:
        prop(matrix_filter, name)

    
    prop(matrix_filter, "HilbertInt").SetElement(0, 1)
    prop(matrix_filter, "ChoiceHilbertDistance").SetElement(0, 2)
    prop(matrix_filter, "L").SetElement(0, 10)
    matrix_filter.SMProxy.UpdateVTKObjects()
    prop(matrix_filter, "ChoiceHilbertDistance").SetElement(0, 4)
    matrix_filter.SMProxy.UpdateVTKObjects()

    pvs.Delete(matrix_filter)
    pvs.Delete(reader)
    print("TTK_MATRIX_FILTER_OK=1")
    print("TTK_SKOT_CHOICE_2_OK=1")
    print("TTK_SK_W2DELTASK_CHOICE_4_OK=1")
    print(f"TEST_DIAGRAM_COUNT={n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
