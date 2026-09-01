#!/usr/bin/env pvpython


from __future__ import annotations

import argparse
import csv
import gc
import json
import math
import os
import re
import shutil
import sys
import time
import traceback
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from ttk_bootstrap import ensure_ttk_plugin_loaded

ensure_ttk_plugin_loaded()

from paraview import servermanager
import paraview.simple as pvs

try:
    import numpy as np
except ImportError as exc:
    raise RuntimeError("NumPy is required and should be included with ParaView.") from exc


DISTANCE_COLUMN_RE = re.compile(r"Diagram(\d+)$")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def set_xml_property(proxy: Any, name: str, value: Any) -> None:
    prop = proxy.SMProxy.GetProperty(name)
    if prop is None:
        raise RuntimeError(
            f"Missing ParaView/TTK property: {name!r}. "
            "Check that the custom TTK plugin is loaded."
        )
    if isinstance(value, (list, tuple)):
        prop.SetNumberOfElements(len(value))
        for index, item in enumerate(value):
            prop.SetElement(index, item)
    else:
        prop.SetElement(0, value)
    proxy.SMProxy.UpdateVTKObjects()


def variant_to_text(value: Any) -> str:
    try:
        return value.ToString()
    except Exception:
        return str(value)


def first_field_values(field_data: Any) -> dict[str, str]:
    result: dict[str, str] = {}
    if field_data is None:
        return result
    for idx in range(field_data.GetNumberOfArrays()):
        array = field_data.GetAbstractArray(idx)
        if array is None or not array.GetName() or array.GetNumberOfTuples() < 1:
            continue
        try:
            result[str(array.GetName())] = variant_to_text(array.GetVariantValue(0))
        except Exception:
            result[str(array.GetName())] = ""
    return result


def parse_vtm_files(path: Path) -> list[str]:
    try:
        root = ET.parse(path).getroot()
    except Exception:
        return []
    datasets: list[tuple[int, str]] = []
    fallback_index = 0
    for element in root.iter():
        if element.tag.split("}")[-1] != "DataSet":
            continue
        filename = element.attrib.get("file", "")
        if not filename:
            continue
        try:
            index = int(element.attrib.get("index", fallback_index))
        except ValueError:
            index = fallback_index
        datasets.append((index, filename))
        fallback_index += 1
    datasets.sort(key=lambda item: item[0])
    return [filename for _, filename in datasets]


def percentile(values: Sequence[float], q: float) -> float:
    if not values:
        return float("nan")
    return float(np.percentile(np.asarray(values, dtype=float), q))


def inspect_collection(collection: Any, input_path: Path, expected: int) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if collection is None or collection.GetClassName() != "vtkMultiBlockDataSet":
        got = collection.GetClassName() if collection is not None else "None"
        raise TypeError(f"Expected input: vtkMultiBlockDataSet; received: {got}.")
    n = int(collection.GetNumberOfBlocks())
    if n != expected:
        raise RuntimeError(f"Unexpected number of diagrams: {n}/{expected}.")

    vtm_files = parse_vtm_files(input_path)
    if vtm_files and len(vtm_files) != n:
        raise RuntimeError(
            f"The VTM manifest references {len(vtm_files)} files for {n} blocks."
        )

    rows: list[dict[str, Any]] = []
    global_min = math.inf
    global_max = -math.inf
    all_counts: list[int] = []
    total_pairs = 0
    normalized_count = 0
    normalization_ranges: list[tuple[float, float]] = []

    for block_index in range(n):
        grid = collection.GetBlock(block_index)
        if grid is None or grid.GetClassName() != "vtkUnstructuredGrid":
            got = grid.GetClassName() if grid is not None else "None"
            raise TypeError(f"Invalid block {block_index}: {got}.")
        cells = grid.GetCellData()
        pair_ids = cells.GetArray("PairIdentifier")
        births = cells.GetArray("Birth")
        pers = cells.GetArray("Persistence")
        pair_types = cells.GetArray("PairType")
        if pair_ids is None or births is None or pers is None:
            raise RuntimeError(
                f"Block {block_index}: PairIdentifier, Birth, or Persistence is missing."
            )
        n_cells = int(pair_ids.GetNumberOfTuples())
        if births.GetNumberOfTuples() != n_cells or pers.GetNumberOfTuples() != n_cells:
            raise RuntimeError(f"Block {block_index}: inconsistent cell arrays.")

        real_count = 0
        type_counts: dict[str, int] = {"pair_type_-1": 0, "pair_type_0": 0, "pair_type_1": 0, "pair_type_2": 0, "pair_type_other": 0}
        persistences: list[float] = []
        local_min = math.inf
        local_max = -math.inf
        for pair_index in range(n_cells):
            if int(round(float(pair_ids.GetTuple1(pair_index)))) == -1:
                continue
            birth = float(births.GetTuple1(pair_index))
            persistence = float(pers.GetTuple1(pair_index))
            death = birth + persistence
            if not (math.isfinite(birth) and math.isfinite(persistence) and math.isfinite(death)):
                raise RuntimeError(f"Block {block_index}, pair {pair_index}: non-finite value.")
            if persistence < -1.0e-12:
                raise RuntimeError(f"Block {block_index}, pair {pair_index}: negative persistence.")
            local_min = min(local_min, birth, death)
            local_max = max(local_max, birth, death)
            persistences.append(max(0.0, persistence))
            real_count += 1
            if pair_types is None:
                type_counts["pair_type_other"] += 1
            else:
                value = int(round(float(pair_types.GetTuple1(pair_index))))
                key = f"pair_type_{value}" if value in (-1, 0, 1, 2) else "pair_type_other"
                type_counts[key] += 1

        if real_count < 1:
            raise RuntimeError(f"Block {block_index} contains no real persistence pair.")
        global_min = min(global_min, local_min)
        global_max = max(global_max, local_max)
        all_counts.append(real_count)
        total_pairs += real_count

        metadata = first_field_values(grid.GetFieldData())
        normalized_marker = metadata.get(
            "NormalizedForSKOT", metadata.get("SKNormalized", "0")
        )
        try:
            is_normalized = int(round(float(normalized_marker))) == 1
        except ValueError:
            is_normalized = False
        if is_normalized:
            normalized_count += 1
            raw_min_text = metadata.get(
                "NormalizationMinimumRaw", metadata.get("SKNormalizationMinRaw")
            )
            raw_max_text = metadata.get(
                "NormalizationMaximumRaw", metadata.get("SKNormalizationMaxRaw")
            )
            try:
                if raw_min_text is not None and raw_max_text is not None:
                    normalization_ranges.append(
                        (float(raw_min_text), float(raw_max_text))
                    )
            except ValueError:
                pass

        row: dict[str, Any] = {
            "sample_index": block_index,
            "vtm_file": vtm_files[block_index] if vtm_files else "",
            "number_of_pairs": real_count,
            "birth_death_min": local_min,
            "birth_death_max": local_max,
            "persistence_min": min(persistences),
            "persistence_median": percentile(persistences, 50),
            "persistence_max": max(persistences),
        }
        row.update(type_counts)
        for key, value in metadata.items():
            row[f"metadata_{key}"] = value
        rows.append(row)

    tolerance = 1.0e-10
    if global_min < -tolerance or global_max > 1.0 + tolerance:
        raise RuntimeError(
            f"The diagrams are not normalized in [0,1]^2: "
            f"detected range [{global_min:.17g}, {global_max:.17g}]."
        )
    if 0 < normalized_count < n:
        raise RuntimeError("The collection mixes normalized and non-normalized diagrams.")
    common_range_confirmed = False
    if normalized_count == n and len(normalization_ranges) == n:
        ref_min, ref_max = normalization_ranges[0]
        tol = 1.0e-12 * max(1.0, abs(ref_min), abs(ref_max))
        common_range_confirmed = all(
            abs(a - ref_min) <= tol and abs(b - ref_max) <= tol
            for a, b in normalization_ranges[1:]
        )
        if not common_range_confirmed:
            raise RuntimeError("The diagrams do not all use the same normalization.")

    summary = {
        "number_of_diagrams": n,
        "number_of_pairs_total": total_pairs,
        "number_of_pairs_min": min(all_counts),
        "number_of_pairs_median": float(np.median(np.asarray(all_counts, dtype=float))),
        "number_of_pairs_max": max(all_counts),
        "birth_death_min": global_min,
        "birth_death_max": global_max,
        "normalized_metadata_count": normalized_count,
        "common_normalization_range_confirmed": common_range_confirmed,
    }
    return rows, summary


def write_rows_csv(rows: Sequence[Mapping[str, Any]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(str(key))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(dict(row))


def table_to_matrix(table: Any) -> tuple[np.ndarray, list[str], list[str]]:
    if table is None or table.GetClassName() != "vtkTable":
        got = table.GetClassName() if table is not None else "None"
        raise RuntimeError(f"The filter did not produce vtkTable ({got}).")
    names = [table.GetColumnName(i) for i in range(table.GetNumberOfColumns())]
    distance_columns = sorted(
        [name for name in names if DISTANCE_COLUMN_RE.fullmatch(name or "")],
        key=lambda name: int(DISTANCE_COLUMN_RE.fullmatch(name).group(1)),
    )
    n_rows = int(table.GetNumberOfRows())
    if n_rows != len(distance_columns):
        raise RuntimeError(
            f"Non-square matrix: {n_rows} rows, {len(distance_columns)} columns."
        )
    matrix = np.empty((n_rows, n_rows), dtype=np.float64)
    for col_index, name in enumerate(distance_columns):
        column = table.GetColumnByName(name)
        for row_index in range(n_rows):
            matrix[row_index, col_index] = float(
                column.GetVariantValue(row_index).ToDouble()
            )
    metadata_columns = [name for name in names if name not in distance_columns]
    return matrix, distance_columns, metadata_columns


def table_metadata_rows(table: Any, columns: Sequence[str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for row_index in range(table.GetNumberOfRows()):
        row = {"sample_index": str(row_index)}
        for name in columns:
            column = table.GetColumnByName(name)
            row[name] = variant_to_text(column.GetVariantValue(row_index))
        rows.append(row)
    return rows


def write_matrix_csv(matrix: np.ndarray, path: Path, labels: Sequence[str] | None = None) -> None:
    n = int(matrix.shape[0])
    if labels is None or len(labels) != n:
        labels = [str(i) for i in range(n)]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Sample"] + list(labels))
        for index, row in enumerate(matrix):
            writer.writerow([labels[index]] + [f"{float(value):.17g}" for value in row])


def max_triangle_violation(matrix: np.ndarray) -> float:
    n = matrix.shape[0]
    maximum = 0.0
    for j in range(n):
        violation = matrix - (matrix[:, j][:, None] + matrix[j, :][None, :])
        candidate = float(np.max(violation))
        if candidate > maximum:
            maximum = candidate
    return maximum


def matrix_diagnostics(matrix: np.ndarray) -> dict[str, Any]:
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise RuntimeError("The matrix is not square.")
    if not np.all(np.isfinite(matrix)):
        raise RuntimeError("The matrix contains NaN or infinity.")
    n = matrix.shape[0]
    scale = max(1.0, float(np.max(np.abs(matrix))))
    max_asym = float(np.max(np.abs(matrix - matrix.T))) if n else 0.0
    max_diag = float(np.max(np.abs(np.diag(matrix)))) if n else 0.0
    min_value = float(np.min(matrix)) if n else 0.0
    if max_asym > 1.0e-9 * scale:
        raise RuntimeError(f"Matrix is too asymmetric: {max_asym:.6g}.")
    if max_diag > 1.0e-9 * scale:
        raise RuntimeError(f"Nonzero diagonal: {max_diag:.6g}.")
    if min_value < -1.0e-10 * scale:
        raise RuntimeError(f"Negative distance: {min_value:.6g}.")

    
    jmat = np.eye(n) - np.ones((n, n), dtype=float) / float(n)
    gram = -0.5 * jmat.dot(matrix * matrix).dot(jmat)
    gram = 0.5 * (gram + gram.T)
    eigenvalues = np.linalg.eigvalsh(gram)
    abs_sum = float(np.sum(np.abs(eigenvalues)))
    negative_mass = float(np.sum(np.abs(eigenvalues[eigenvalues < 0.0])))
    return {
        "shape": [int(n), int(n)],
        "minimum": min_value,
        "maximum": float(np.max(matrix)) if n else 0.0,
        "max_asymmetry": max_asym,
        "max_absolute_diagonal": max_diag,
        "max_triangle_violation": max_triangle_violation(matrix),
        "gram_min_eigenvalue": float(eigenvalues[0]) if n else 0.0,
        "gram_negative_eigenvalue_mass_ratio": negative_mass / abs_sum if abs_sum > 0 else 0.0,
    }


def method_complete(
    folder: Path,
    expected: int,
    need_skot_cost: bool,
    method: str,
    level: int | None,
    delta_lim: float,
) -> bool:
    try:
        if not (folder / "DONE.ok").is_file():
            return False
        matrix_path = folder / "matrix_distance.npy"
        info_path = folder / "info.json"
        if not (
            matrix_path.is_file()
            and (folder / "matrix_distance.csv").is_file()
            and (folder / "matrix_distance_labeled.csv").is_file()
            and info_path.is_file()
        ):
            return False
        info = json.loads(info_path.read_text(encoding="utf-8"))
        if info.get("method") != method:
            return False
        if info.get("level") != (int(level) if level is not None else None):
            return False
        parameters = info.get("parameters", {})
        if method == "W2" and abs(
            float(parameters.get("delta_lim", float("nan"))) - float(delta_lim)
        ) > 1.0e-15:
            return False
        if parameters.get("critical_pairs_xml_value") != -1:
            return False
        if parameters.get("constraint_xml_value") != 0:
            return False
        matrix = np.load(matrix_path, allow_pickle=False)
        if matrix.shape != (expected, expected) or not np.all(np.isfinite(matrix)):
            return False
        if need_skot_cost:
            if not (
                (folder / "matrix_SKOT_cost.npy").is_file()
                and (folder / "matrix_SKOT_cost.csv").is_file()
                and (folder / "matrix_SKOT_cost_labeled.csv").is_file()
            ):
                return False
            cost = np.load(folder / "matrix_SKOT_cost.npy", allow_pickle=False)
            if cost.shape != (expected, expected) or not np.all(np.isfinite(cost)):
                return False
        return True
    except Exception:
        return False


def configure_filter(
    proxy: Any,
    method: str,
    level: int | None,
    delta_lim: float,
) -> dict[str, Any]:
    
    set_xml_property(proxy, "Critical pairs", -1)
    set_xml_property(proxy, "Constraint", 0)
    set_xml_property(proxy, "MaxNumberOfPairs", 20)
    set_xml_property(proxy, "n", "2")
    set_xml_property(proxy, "DeltaLim", float(delta_lim))
    set_xml_property(proxy, "AntiAlpha", 0.0)  
    set_xml_property(proxy, "Lambda", 1.0)

    parameters: dict[str, Any] = {
        "critical_pairs": "all",
        "critical_pairs_xml_value": -1,
        "constraint": "full_diagrams",
        "constraint_xml_value": 0,
        "wasserstein_p": 2,
        "delta_lim": float(delta_lim),
        "anti_alpha": 0.0,
        "internal_alpha": 1.0,
        "lambda": 1.0,
    }
    if method == "SKOT":
        if level is None:
            raise ValueError("Missing L for SKOT")
        set_xml_property(proxy, "HilbertInt", 1)
        set_xml_property(proxy, "ChoiceHilbertDistance", 2)
        set_xml_property(proxy, "L", int(level))
        set_xml_property(proxy, "GLevel", 0)
        parameters.update(
            {
                "use_sk_distance": True,
                "choice_hilbert_distance": 2,
                "choice_name": "SK-OT",
                "L": int(level),
            }
        )
    elif method == "SK_W2DeltaSk":
        if level is None:
            raise ValueError("Missing L for SK_W2DeltaSk")
        set_xml_property(proxy, "HilbertInt", 1)
        set_xml_property(proxy, "ChoiceHilbertDistance", 4)
        set_xml_property(proxy, "L", int(level))
        set_xml_property(proxy, "GLevel", 0)
        parameters.update(
            {
                "use_sk_distance": True,
                "choice_hilbert_distance": 4,
                "choice_name": "SK-W2DeltaSk",
                "L": int(level),
            }
        )
    elif method == "W2":
        set_xml_property(proxy, "HilbertInt", 0)
        parameters.update(
            {
                "use_sk_distance": False,
                "choice_name": "TTK classical 2-Wasserstein",
            }
        )
    else:
        raise ValueError(method)
    return parameters


def run_method(
    reader: Any,
    collection_output: Path,
    method: str,
    level: int | None,
    expected: int,
    labels: Sequence[str],
    delta_lim: float,
) -> dict[str, Any]:
    if method == "SKOT":
        final_dir = collection_output / "SKOT" / f"L{int(level):02d}"
        display_name = f"sqrt(SKOT), L={level}"
        semantics = "TTK filter output d_SK,L = sqrt(SKOT_L)"
        need_skot_cost = True
    elif method == "SK_W2DeltaSk":
        final_dir = collection_output / "SK_W2DeltaSk" / f"L{int(level):02d}"
        display_name = f"SK-W2DeltaSk, L={level}"
        semantics = "square root of the SK-induced squared 2D matching cost"
        need_skot_cost = False
    else:
        final_dir = collection_output / "W2"
        display_name = "classical TTK W2"
        semantics = "TTK persistence-diagram 2-Wasserstein distance"
        need_skot_cost = False

    if method_complete(
        final_dir, expected, need_skot_cost, method, level, delta_lim
    ):
        info = json.loads((final_dir / "info.json").read_text(encoding="utf-8"))
        print(f"  already computed: {display_name}")
        return info

    if final_dir.exists():
        shutil.rmtree(final_dir)
    temp_dir = final_dir.parent / f".{final_dir.name}.in_progress"
    if temp_dir.exists():
        shutil.rmtree(temp_dir)
    temp_dir.mkdir(parents=True, exist_ok=True)

    matrix_filter = None
    try:
        print(f"\n  COMPUTATION: {display_name}")
        matrix_filter = pvs.TTKPersistenceDiagramDistanceMatrix(Input=reader)
        parameters = configure_filter(matrix_filter, method, level, delta_lim)

        start_compute = time.perf_counter()
        matrix_filter.UpdatePipeline()
        compute_seconds = time.perf_counter() - start_compute

        start_fetch = time.perf_counter()
        table = servermanager.Fetch(matrix_filter)
        fetch_seconds = time.perf_counter() - start_fetch

        matrix, _, metadata_columns = table_to_matrix(table)
        if matrix.shape != (expected, expected):
            raise RuntimeError(
                f"Matrix {matrix.shape}, whereas {(expected, expected)} was expected."
            )
        diagnostics = matrix_diagnostics(matrix)

        start_save = time.perf_counter()
        np.save(temp_dir / "matrix_distance.npy", matrix)
        write_matrix_csv(matrix, temp_dir / "matrix_distance.csv")
        write_matrix_csv(matrix, temp_dir / "matrix_distance_labeled.csv", labels)

        if method == "SKOT":
            skot_cost = matrix * matrix
            np.save(temp_dir / "matrix_SKOT_cost.npy", skot_cost)
            write_matrix_csv(skot_cost, temp_dir / "matrix_SKOT_cost.csv")
            write_matrix_csv(
                skot_cost, temp_dir / "matrix_SKOT_cost_labeled.csv", labels
            )

        ttk_metadata_path = collection_output / "samples_from_TTK_table.csv"
        if not ttk_metadata_path.is_file() and metadata_columns:
            write_rows_csv(table_metadata_rows(table, metadata_columns), ttk_metadata_path)

        vtt_saved = False
        try:
            pvs.SaveData(str(temp_dir / "matrix_table.vtt"), proxy=matrix_filter)
            vtt_saved = (temp_dir / "matrix_table.vtt").is_file()
        except Exception as exc:
            (temp_dir / "VTT_WARNING.txt").write_text(
                f"The optional VTT table could not be saved: {exc}\n",
                encoding="utf-8",
            )
        save_seconds = time.perf_counter() - start_save

        info = {
            "status": "complete",
            "method": method,
            "display_name": display_name,
            "matrix_semantics": semantics,
            "level": int(level) if level is not None else None,
            "computed_at_utc": utc_now(),
            "compute_seconds": compute_seconds,
            "fetch_seconds": fetch_seconds,
            "save_seconds": save_seconds,
            "total_seconds_excluding_input_read": compute_seconds + fetch_seconds + save_seconds,
            "vtt_saved": vtt_saved,
            "parameters": parameters,
            "diagnostics": diagnostics,
            "files": {
                "distance_npy": "matrix_distance.npy",
                "distance_csv": "matrix_distance.csv",
                "distance_labeled_csv": "matrix_distance_labeled.csv",
            },
        }
        if method == "SKOT":
            info["files"].update(
                {
                    "skot_cost_npy": "matrix_SKOT_cost.npy",
                    "skot_cost_csv": "matrix_SKOT_cost.csv",
                    "skot_cost_labeled_csv": "matrix_SKOT_cost_labeled.csv",
                }
            )
        (temp_dir / "info.json").write_text(
            json.dumps(info, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )
        (temp_dir / "DONE.ok").write_text(
            f"{display_name} computed and validated.\n", encoding="utf-8"
        )
        final_dir.parent.mkdir(parents=True, exist_ok=True)
        temp_dir.rename(final_dir)
        print(
            f"  completed in {compute_seconds:.3f} s (filter computation), "
            f"max={diagnostics['maximum']:.6g}"
        )
        return info
    except Exception as exc:
        (temp_dir / "ERROR.txt").write_text(
            f"{type(exc).__name__}: {exc}\n\n{traceback.format_exc()}",
            encoding="utf-8",
        )
        raise
    finally:
        if matrix_filter is not None:
            try:
                pvs.Delete(matrix_filter)
            except Exception:
                pass
        gc.collect()


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
    parser.add_argument("--source-data-csv", type=Path, default=None)
    parser.add_argument("--normalization-csv", type=Path, default=None)
    return parser.parse_args(argv)


def main(argv: Sequence[str] = sys.argv[1:]) -> int:
    args = parse_args(argv)
    input_path = args.input.expanduser().resolve()
    output = args.output.expanduser().resolve()
    if not input_path.is_file():
        raise FileNotFoundError(input_path)
    levels = sorted(set(int(value) for value in args.levels))
    if any(value < 0 or value > 50 for value in levels):
        raise ValueError("L levels must be between 0 and 50.")

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
    print("=" * 78)

    start_read = time.perf_counter()
    reader = pvs.XMLMultiBlockDataReader(FileName=[str(input_path)])
    reader.UpdatePipeline()
    collection = servermanager.Fetch(reader)
    input_read_seconds = time.perf_counter() - start_read

    sample_rows, collection_summary = inspect_collection(
        collection, input_path, int(args.expected)
    )
    
    
    
    del collection
    gc.collect()
    labels = [
        str(row.get("vtm_file") or f"sample_{int(row['sample_index']):03d}")
        for row in sample_rows
    ]
    write_rows_csv(sample_rows, output / "samples.csv")
    collection_info = {
        "dataset_key": args.dataset_key,
        "dataset_label": args.dataset_label,
        "input_manifest": str(input_path),
        "expected_diagrams": int(args.expected),
        "levels": levels,
        "input_read_seconds": input_read_seconds,
        "collection_summary": collection_summary,
        "created_at_utc": utc_now(),
    }
    (output / "collection_info.json").write_text(
        json.dumps(collection_info, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(
        f"  {collection_summary['number_of_diagrams']} diagrams, "
        f"{collection_summary['number_of_pairs_total']} total pairs, "
        f"range [{collection_summary['birth_death_min']:.6g}, "
        f"{collection_summary['birth_death_max']:.6g}]"
    )

    method_infos: list[dict[str, Any]] = []
    
    
    for level in levels:
        method_infos.append(
            run_method(
                reader,
                output,
                "SKOT",
                level,
                int(args.expected),
                labels,
                float(args.delta_lim),
            )
        )
    for level in levels:
        method_infos.append(
            run_method(
                reader,
                output,
                "SK_W2DeltaSk",
                level,
                int(args.expected),
                labels,
                float(args.delta_lim),
            )
        )
    method_infos.append(
        run_method(
            reader,
            output,
            "W2",
            None,
            int(args.expected),
            labels,
            float(args.delta_lim),
        )
    )

    (output / "METHODS_COMPLETE.json").write_text(
        json.dumps(method_infos, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    (output / "COLLECTION_COMPLETE.ok").write_text(
        f"All matrices for {args.dataset_label} are complete.\n",
        encoding="utf-8",
    )
    try:
        pvs.Delete(reader)
    except Exception:
        pass
    print(f"\nCOLLECTION_COMPLETE={args.dataset_key}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
