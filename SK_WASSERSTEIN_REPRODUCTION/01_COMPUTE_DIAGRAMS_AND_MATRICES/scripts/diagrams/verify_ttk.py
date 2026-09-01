#!/usr/bin/env python3


from __future__ import annotations

import argparse
import os
import sys
import tempfile
from pathlib import Path


def load_plugin() -> None:
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


def new_ttk_pd(input_proxy=None):
    import paraview.simple as pvs
    from paraview import servermanager

    constructor = getattr(pvs, "TTKPersistenceDiagram", None)
    if constructor is not None:
        if input_proxy is None:
            return constructor()
        return constructor(Input=input_proxy)

    sm_proxy = servermanager.ProxyManager().NewProxy(
        "filters", "ttkPersistenceDiagram"
    )
    if sm_proxy is None:
        raise RuntimeError("The ttkPersistenceDiagram proxy is not registered.")
    py_proxy = servermanager._getPyProxy(sm_proxy)
    if input_proxy is not None:
        py_proxy.Input = input_proxy
    return py_proxy


def set_scalar(proxy, name: str) -> None:
    errors = []
    for prop in ("ScalarField", "ScalarFieldNew"):
        for value in (["POINTS", name], name):
            try:
                setattr(proxy, prop, value)
                return
            except Exception as exc:  
                errors.append(f"{prop}={value!r}: {exc}")
    raise RuntimeError(
        "Could not select the scalar field. " + " | ".join(errors)
    )


def basic_check() -> bool:
    from paraview import servermanager

    manager = servermanager.ProxyManager()
    diagram = manager.NewProxy("filters", "ttkPersistenceDiagram")
    if diagram is None:
        raise RuntimeError("TTK PersistenceDiagram is unavailable.")

    matrix = manager.NewProxy("filters", "ttkPersistenceDiagramDistanceMatrix")
    custom_skot = bool(
        matrix is not None and matrix.GetProperty("ChoiceHilbertDistance") is not None
    )
    return custom_skot


def full_check() -> None:
    import paraview.simple as pvs
    from paraview import servermanager

    source = pvs.Wavelet()

    
    
    source.UpdatePipeline()
    point_info = source.GetPointDataInformation()
    source_arrays = {}
    for index in range(point_info.GetNumberOfArrays()):
        array = point_info.GetArray(index)
        if array is not None and array.GetName():
            source_arrays[str(array.GetName())] = int(array.GetNumberOfComponents())
    if source_arrays.get("RTData") != 1:
        raise RuntimeError(
            "The ParaView Point Data API did not return RTData correctly. "
            f"Detected arrays: {source_arrays}"
        )

    
    
    
    
    
    try:
        import vtk
    except ImportError:  
        from vtkmodules import all as vtk  
    from vtk_reader_utils import open_paraview_reader, reader_spec

    required_xml_readers = (
        "XMLImageDataReader",
        "XMLRectilinearGridReader",
        "XMLStructuredGridReader",
        "XMLPolyDataReader",
        "XMLUnstructuredGridReader",
    )
    missing_readers = [name for name in required_xml_readers if getattr(pvs, name, None) is None]
    if missing_readers:
        raise RuntimeError(
            "Missing ParaView 5.13 XML readers: " + ", ".join(missing_readers)
        )

    with tempfile.TemporaryDirectory(prefix="skw_reader_test_") as temporary:
        disguised = Path(temporary) / "image_payload_with_vtu_suffix.vtu"
        image = vtk.vtkImageData()
        image.SetDimensions(2, 2, 2)
        values = vtk.vtkFloatArray()
        values.SetName("ReaderSentinel")
        values.SetNumberOfTuples(8)
        for index in range(8):
            values.SetValue(index, float(index))
        image.GetPointData().AddArray(values)
        writer = vtk.vtkXMLImageDataWriter()
        writer.SetFileName(str(disguised))
        writer.SetInputData(image)
        if writer.Write() != 1:
            raise RuntimeError("Could not create the VTK reader test file.")

        spec = reader_spec(disguised)
        if spec is None or spec.paraview_constructor != "XMLImageDataReader":
            raise RuntimeError(
                "Actual VTK type detection failed for a .vtu file "
                "containing ImageData."
            )
        disguised_reader = open_paraview_reader(disguised, pvs)
        disguised_reader.UpdatePipeline()
        disguised_info = disguised_reader.GetPointDataInformation()
        disguised_arrays = {
            str(disguised_info.GetArray(i).GetName())
            for i in range(disguised_info.GetNumberOfArrays())
            if disguised_info.GetArray(i) is not None
            and disguised_info.GetArray(i).GetName()
        }
        if "ReaderSentinel" not in disguised_arrays:
            raise RuntimeError(
                "The reader selected from the XML header did not read "
                f"ReaderSentinel. Detected arrays: {sorted(disguised_arrays)}"
            )
        try:
            pvs.Delete(disguised_reader)
        except Exception:
            pass

    diagram = new_ttk_pd(source)
    set_scalar(diagram, "RTData")
    for prop, value in (
        ("BackEnd", 2),
        ("DMSDimensions", 0),
        ("ShowInsideDomain", 0),
        ("ClearDGCache", 1),
        ("ForceInputOffsetScalarField", 0),
    ):
        try:
            setattr(diagram, prop, value)
        except Exception:
            pass

    diagram.UpdatePipeline()
    data = servermanager.Fetch(diagram)
    if data is None or not data.IsA("vtkUnstructuredGrid"):
        raise RuntimeError("The TTK filter did not produce a vtkUnstructuredGrid.")
    cell_data = data.GetCellData()
    point_data = data.GetPointData()
    required_cells = ("PairIdentifier", "PairType", "Persistence", "Birth")
    required_points = ("ttkVertexScalarField", "CriticalType")
    missing = [name for name in required_cells if cell_data.GetArray(name) is None]
    missing += [name for name in required_points if point_data.GetArray(name) is None]
    if missing:
        raise RuntimeError("Missing TTK arrays: " + ", ".join(missing))

    try:
        pvs.Delete(diagram)
        pvs.Delete(source)
    except Exception:
        pass


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--full", action="store_true")
    args = parser.parse_args()

    try:
        load_plugin()
        custom_skot = basic_check()
        if args.full:
            full_check()
    except Exception as exc:
        print(f"ERROR={exc}")
        return 1

    print("TTK_PD_OK=1")
    print(f"TTK_SKOT_CUSTOM={int(custom_skot)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
