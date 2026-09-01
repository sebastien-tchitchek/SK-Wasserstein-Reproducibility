#!/usr/bin/env python3


from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional


@dataclass(frozen=True)
class ReaderSpec:
    dataset_type: str
    paraview_constructor: str
    vtk_constructor: str
    detected_from_header: bool


def _normalized_type(value: str) -> str:
    key = re.sub(r"[^a-z0-9]+", "", value.casefold())
    if key.startswith("vtk"):
        key = key[3:]
    return key




_XML_TYPE_SPECS = {
    "imagedata": ("XMLImageDataReader", "vtkXMLImageDataReader"),
    "pimagedata": ("XMLPartitionedImageDataReader", "vtkXMLPImageDataReader"),
    "rectilineargrid": ("XMLRectilinearGridReader", "vtkXMLRectilinearGridReader"),
    "prectilineargrid": (
        "XMLPartitionedRectilinearGridReader",
        "vtkXMLPRectilinearGridReader",
    ),
    "structuredgrid": ("XMLStructuredGridReader", "vtkXMLStructuredGridReader"),
    "pstructuredgrid": (
        "XMLPartitionedStructuredGridReader",
        "vtkXMLPStructuredGridReader",
    ),
    "polydata": ("XMLPolyDataReader", "vtkXMLPolyDataReader"),
    "ppolydata": ("XMLPartitionedPolydataReader", "vtkXMLPPolyDataReader"),
    "unstructuredgrid": (
        "XMLUnstructuredGridReader",
        "vtkXMLUnstructuredGridReader",
    ),
    "punstructuredgrid": (
        "XMLPartitionedUnstructuredGridReader",
        "vtkXMLPUnstructuredGridReader",
    ),
    "multiblockdataset": ("XMLMultiBlockDataReader", "vtkXMLMultiBlockDataReader"),
}

_SUFFIX_SPECS = {
    ".vti": ("ImageData", "XMLImageDataReader", "vtkXMLImageDataReader"),
    ".pvti": (
        "PImageData",
        "XMLPartitionedImageDataReader",
        "vtkXMLPImageDataReader",
    ),
    ".vtr": (
        "RectilinearGrid",
        "XMLRectilinearGridReader",
        "vtkXMLRectilinearGridReader",
    ),
    ".pvtr": (
        "PRectilinearGrid",
        "XMLPartitionedRectilinearGridReader",
        "vtkXMLPRectilinearGridReader",
    ),
    ".vts": ("StructuredGrid", "XMLStructuredGridReader", "vtkXMLStructuredGridReader"),
    ".pvts": (
        "PStructuredGrid",
        "XMLPartitionedStructuredGridReader",
        "vtkXMLPStructuredGridReader",
    ),
    ".vtp": ("PolyData", "XMLPolyDataReader", "vtkXMLPolyDataReader"),
    ".pvtp": (
        "PPolyData",
        "XMLPartitionedPolydataReader",
        "vtkXMLPPolyDataReader",
    ),
    ".vtu": (
        "UnstructuredGrid",
        "XMLUnstructuredGridReader",
        "vtkXMLUnstructuredGridReader",
    ),
    ".pvtu": (
        "PUnstructuredGrid",
        "XMLPartitionedUnstructuredGridReader",
        "vtkXMLPUnstructuredGridReader",
    ),
    ".vtm": ("vtkMultiBlockDataSet", "XMLMultiBlockDataReader", "vtkXMLMultiBlockDataReader"),
    ".vtk": ("LegacyVTK", "LegacyVTKReader", "vtkDataSetReader"),
}

_VTKFILE_TYPE_RE = re.compile(
    rb"<\s*VTKFile\b[^>]*\btype\s*=\s*[\"']([^\"']+)[\"']",
    flags=re.IGNORECASE | re.DOTALL,
)


def read_file_prefix(path: Path, maximum: int = 1024 * 1024) -> bytes:
    try:
        with Path(path).open("rb") as stream:
            return stream.read(maximum)
    except OSError:
        return b""


def detect_vtk_xml_type(path: Path) -> Optional[str]:
    
    match = _VTKFILE_TYPE_RE.search(read_file_prefix(Path(path)))
    if match is None:
        return None
    return match.group(1).decode("ascii", errors="replace").strip()


def is_legacy_vtk(path: Path) -> bool:
    prefix = read_file_prefix(Path(path), maximum=512).lstrip().casefold()
    return prefix.startswith(b"# vtk datafile version")


def reader_spec(path: Path) -> Optional[ReaderSpec]:
    path = Path(path)
    xml_type = detect_vtk_xml_type(path)
    if xml_type:
        names = _XML_TYPE_SPECS.get(_normalized_type(xml_type))
        if names is not None:
            return ReaderSpec(xml_type, names[0], names[1], True)
    if is_legacy_vtk(path):
        return ReaderSpec("LegacyVTK", "LegacyVTKReader", "vtkDataSetReader", True)
    fallback = _SUFFIX_SPECS.get(path.suffix.casefold())
    if fallback is None:
        return None
    return ReaderSpec(fallback[0], fallback[1], fallback[2], False)


def open_paraview_reader(path: Path, pvs: Any) -> Any:
    
    path = Path(path).expanduser().resolve()
    spec = reader_spec(path)
    errors: list[str] = []

    if spec is not None:
        constructor = getattr(pvs, spec.paraview_constructor, None)
        if constructor is None:
            errors.append(
                f"proxy {spec.paraview_constructor} missing for type {spec.dataset_type}"
            )
        else:
            
            
            attempts = (
                ("FileName", [str(path)]),
                ("FileName", str(path)),
                ("FileNames", [str(path)]),
                ("FileNames", str(path)),
            )
            for keyword, value in attempts:
                try:
                    return constructor(**{keyword: value})
                except Exception as exc:  
                    errors.append(f"{spec.paraview_constructor}({keyword}) : {exc}")

        
        
        if spec.detected_from_header:
            raise RuntimeError(
                f"Could not open {path.name} as {spec.dataset_type} with "
                f"{spec.paraview_constructor}. " + " | ".join(errors)
            )

    opener = getattr(pvs, "OpenDataFile", None)
    if opener is None:
        details = " | ".join(errors) if errors else "no specification detected"
        raise RuntimeError(f"No ParaView reader is available for {path}: {details}")
    reader = opener(str(path))
    if reader is None:
        details = " | ".join(errors) if errors else ""
        raise RuntimeError(f"ParaView cannot open {path}. {details}")
    return reader


def make_vtk_reader(path: Path, vtk_module: Any) -> Any:
    
    path = Path(path).expanduser().resolve()
    spec = reader_spec(path)
    if spec is None:
        raise ValueError(f"Unsupported VTK format: {path}")
    constructor = getattr(vtk_module, spec.vtk_constructor, None)
    if constructor is None:
        
        if spec.dataset_type != "LegacyVTK":
            constructor = getattr(vtk_module, "vtkXMLGenericDataObjectReader", None)
    if constructor is None:
        raise RuntimeError(
            f"VTK reader {spec.vtk_constructor} unavailable for {path.name}"
        )
    reader = constructor()
    reader.SetFileName(str(path))
    return reader


def format_description(path: Path) -> str:
    spec = reader_spec(Path(path))
    if spec is None:
        return "inconnu"
    origin = "XML header" if spec.detected_from_header else "extension"
    return f"{spec.dataset_type} ({origin}, reader {spec.paraview_constructor})"
