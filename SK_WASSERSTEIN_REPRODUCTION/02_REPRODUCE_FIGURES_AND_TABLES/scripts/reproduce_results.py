#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
from pathlib import Path
from typing import Any, Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import spearmanr
from sklearn.cluster import AgglomerativeClustering, KMeans, SpectralClustering
from sklearn.metrics import adjusted_rand_score

LEVELS = (10, 20, 30, 40)
PRIMARY_L = 30
REFERENCE_L = 40

COLLECTIONS = (
    ("2004_isabel_3D", "Isabel", 12, 3),
    ("2006_earthquake_3D", "Earthquake", 12, 3),
    ("2008_ionization_front_2D", "Ionization 2D", 16, 4),
    ("2008_ionization_front_3D", "Ionization 3D", 16, 4),
    ("2014_volcanic_eruptions_2D", "Volcanic", 12, 3),
    ("2016_viscous_fingering_3D", "Viscous Fingering", 15, 3),
    ("2017_cloud_processes_2D", "Cloud Processes", 12, 3),
    ("2018_asteroid_impact_3D_clustering", "Asteroid ensemble", 7, 2),
    ("2018_asteroid_impact_3D_temporal_subsampling", "Asteroid temporal", 20, 4),
    ("sea_surface_height", "Sea Surface Height", 48, 4),
    ("starting_vortex", "Starting Vortex", 12, 2),
    ("vortex_street", "Vortex Street", 45, 5),
)
LABELS = {key: label for key, label, _, _ in COLLECTIONS}
EXPECTED_DIAGRAMS = 227
EXPECTED_PERSISTENCE_PAIRS = 3_376_524
EXPECTED_COMPARISONS = 3_004


def upper(matrix: np.ndarray) -> np.ndarray:
    return matrix[np.triu_indices_from(matrix, k=1)]


def matrix_path(root: Path, key: str, method: str, level: int | None = None) -> Path:
    if method == "W2":
        return root / key / "W2" / "matrix_distance.npy"
    return root / key / method / f"L{int(level):02d}" / "matrix_distance.npy"


def info_path(root: Path, key: str, method: str, level: int | None = None) -> Path:
    if method == "W2":
        return root / key / "W2" / "info.json"
    return root / key / method / f"L{int(level):02d}" / "info.json"


def load_matrix(root: Path, key: str, method: str, level: int | None = None) -> np.ndarray:
    return np.load(matrix_path(root, key, method, level), allow_pickle=False)


def load_info(root: Path, key: str, method: str, level: int | None = None) -> dict[str, Any]:
    return json.loads(info_path(root, key, method, level).read_text(encoding="utf-8"))


def compute_time(info: dict[str, Any]) -> float:
    if "compute_seconds_mean" in info:
        return float(info["compute_seconds_mean"])
    runs = info.get("timing_runs_seconds")
    if isinstance(runs, list) and runs:
        return float(np.mean(np.asarray(runs, dtype=float)))
    return float(info["compute_seconds"])


def validate_matrix(matrix: np.ndarray, n: int, name: str) -> None:
    if matrix.shape != (n, n):
        raise RuntimeError(f"{name}: shape {matrix.shape}, expected {(n, n)}")
    if not np.all(np.isfinite(matrix)):
        raise RuntimeError(f"{name}: non-finite values")
    scale = max(1.0, float(np.max(np.abs(matrix))))
    if float(np.max(np.abs(matrix - matrix.T))) > 1e-9 * scale:
        raise RuntimeError(f"{name}: non-symmetric matrix")
    if float(np.max(np.abs(np.diag(matrix)))) > 1e-9 * scale:
        raise RuntimeError(f"{name}: nonzero diagonal")
    if float(np.min(matrix)) < -1e-10 * scale:
        raise RuntimeError(f"{name}: negative distance")


def labels_from_samples(path: Path) -> np.ndarray:
    frame = pd.read_csv(path)
    candidates = [c for c in frame.columns if "clusterid" in c.lower()]
    if not candidates:
        raise RuntimeError(f"No ClusterID column in {path}")
    return pd.to_numeric(frame[candidates[0]], errors="raise").astype(int).to_numpy()


def agglomerative_precomputed(matrix: np.ndarray, k: int) -> np.ndarray:
    return AgglomerativeClustering(
        n_clusters=k, metric="precomputed", linkage="average"
    ).fit_predict(matrix)


def centered_gram(matrix: np.ndarray) -> np.ndarray:
    n = matrix.shape[0]
    centering = np.eye(n) - np.ones((n, n), dtype=float) / n
    gram = -0.5 * centering @ (matrix * matrix) @ centering
    return 0.5 * (gram + gram.T)


def hilbert_coordinates(matrix: np.ndarray, tol_factor: float = 1e-12) -> np.ndarray:
    values, vectors = np.linalg.eigh(centered_gram(matrix))
    order = np.argsort(values)[::-1]
    values = values[order]
    vectors = vectors[:, order]
    tolerance = tol_factor * max(1.0, float(np.max(np.abs(values))))
    keep = values > tolerance
    return vectors[:, keep] * np.sqrt(values[keep])


def nearest_overlap(a: np.ndarray, b: np.ndarray, k: int = 3) -> float:
    scores: list[float] = []
    for i in range(a.shape[0]):
        aa = a[i].copy()
        bb = b[i].copy()
        aa[i] = np.inf
        bb[i] = np.inf
        ia = np.argsort(aa)[:k]
        ib = np.argsort(bb)[:k]
        scores.append(len(set(ia).intersection(ib)) / k)
    return float(np.mean(scores))


def median_sigma(matrix: np.ndarray) -> float:
    values = upper(matrix)
    values = values[values > 0]
    if not len(values):
        raise RuntimeError("Cannot define sigma: all distances are zero")
    return float(np.median(values))


def gaussian(matrix: np.ndarray, sigma: float) -> np.ndarray:
    return np.exp(-(matrix * matrix) / (2.0 * sigma * sigma))


def spectral_partition(affinity: np.ndarray, k: int) -> np.ndarray:
    return SpectralClustering(
        n_clusters=k,
        affinity="precomputed",
        assign_labels="kmeans",
        n_init=100,
        random_state=0,
    ).fit_predict(affinity)


def true_boundaries(labels: np.ndarray) -> list[int]:
    return [i for i in range(1, len(labels)) if labels[i] != labels[i - 1]]


def segment_kernel(kernel: np.ndarray, n_segments: int) -> tuple[np.ndarray, list[int]]:
    n = kernel.shape[0]
    diag_prefix = np.r_[0.0, np.cumsum(np.diag(kernel))]
    summed = np.pad(kernel, ((1, 0), (1, 0))).cumsum(0).cumsum(1)

    def block_sum(a: int, b: int) -> float:
        return float(summed[b, b] - summed[a, b] - summed[b, a] + summed[a, a])

    def cost(a: int, b: int) -> float:
        m = b - a
        return float(diag_prefix[b] - diag_prefix[a] - block_sum(a, b) / m)

    dp = np.full((n_segments + 1, n + 1), np.inf)
    previous = np.full((n_segments + 1, n + 1), -1, dtype=int)
    dp[0, 0] = 0.0
    for q in range(1, n_segments + 1):
        for b in range(q, n + 1):
            for a in range(q - 1, b):
                candidate = dp[q - 1, a] + cost(a, b)
                if candidate < dp[q, b]:
                    dp[q, b] = candidate
                    previous[q, b] = a
    boundaries: list[int] = []
    b = n
    for q in range(n_segments, 1, -1):
        b = int(previous[q, b])
        boundaries.append(b)
    boundaries.sort()
    prediction = np.zeros(n, dtype=int)
    starts = [0] + boundaries
    ends = boundaries + [n]
    for label, (a, b) in enumerate(zip(starts, ends)):
        prediction[a:b] = label
    return prediction, boundaries


def boundary_mae(reference: Sequence[int], predicted: Sequence[int]) -> float:
    if not reference:
        return 0.0
    return float(np.mean(np.abs(np.asarray(reference) - np.asarray(predicted))))


def tex_int(value: int | float) -> str:
    return f"{int(round(value)):,}".replace(",", r"\,")


def tex_float(value: float, digits: int = 3) -> str:
    return f"{float(value):.{digits}f}"


def tex_sci(value: float, digits: int = 2) -> str:
    if value == 0:
        return "0"
    exponent = int(math.floor(math.log10(abs(value))))
    mantissa = value / (10**exponent)
    return rf"{mantissa:.{digits}f}\times10^{{{exponent}}}"


def format_boundaries(values: Sequence[int]) -> str:
    return "[" + ", ".join(str(int(v)) for v in values) + "]"


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip() + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def bold_best(values: Sequence[float], index: int) -> str:
    formatted = f"{values[index]:.3f}"
    if math.isclose(values[index], max(values), rel_tol=0, abs_tol=5e-13):
        return rf"\textbf{{{formatted}}}"
    return formatted


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def locate_matrix_root(path: Path) -> Path:
    path = path.expanduser().resolve()
    if all((path / key).is_dir() for key in LABELS):
        return path
    candidates = [p for p in path.rglob("SKOT_WGAMMA_W2_MATRICES_12_COLLECTIONS") if p.is_dir()]
    if len(candidates) == 1:
        return candidates[0]
    raise RuntimeError(f"Could not identify the matrix root in {path}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrices", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--paper-root", type=Path, default=None)
    args = parser.parse_args()

    root = locate_matrix_root(args.matrices)
    output = args.output.expanduser().resolve()
    csv_dir = output / "csv"
    fig_dir = output / "figures"
    table_dir = output / "tables"
    for directory in (csv_dir, fig_dir, table_dir):
        directory.mkdir(parents=True, exist_ok=True)

    dataset_rows: list[dict[str, Any]] = []
    convergence_rows: list[dict[str, Any]] = []
    w2_rows: list[dict[str, Any]] = []
    clustering_rows: list[dict[str, Any]] = []
    gaussian_rows: list[dict[str, Any]] = []
    segmentation_rows: list[dict[str, Any]] = []
    all_dsk_pairs: list[np.ndarray] = []
    all_wgamma_pairs: list[np.ndarray] = []
    all_w2_pairs: list[np.ndarray] = []
    asteroid_coordinates: np.ndarray | None = None
    asteroid_kernel: np.ndarray | None = None
    asteroid_labels: np.ndarray | None = None
    asteroid_reference_bounds: list[int] = []
    asteroid_predicted_bounds: list[int] = []

    for key, label, expected_n, expected_k in COLLECTIONS:
        samples_path = root / key / "samples.csv"
        samples = pd.read_csv(samples_path)
        labels = labels_from_samples(samples_path)
        if len(samples) != expected_n or len(labels) != expected_n:
            raise RuntimeError(f"{key}: inconsistent sample count")
        unique_labels = np.unique(labels)
        if len(unique_labels) != expected_k:
            raise RuntimeError(f"{key}: {len(unique_labels)} groups, expected {expected_k}")
        sizes = pd.to_numeric(samples["number_of_pairs"], errors="raise").astype(int).to_numpy()
        dataset_rows.append(
            {
                "dataset_key": key,
                "dataset": label,
                "n": expected_n,
                "k": expected_k,
                "pairs_min": int(sizes.min()),
                "pairs_median": float(np.median(sizes)),
                "pairs_max": int(sizes.max()),
                "pairs_total": int(sizes.sum()),
            }
        )

        matrices: dict[tuple[str, int | None], np.ndarray] = {}
        for method in ("SKOT", "SK_W2DeltaSk"):
            for level in LEVELS:
                matrix = load_matrix(root, key, method, level)
                validate_matrix(matrix, expected_n, f"{key}/{method}/L{level}")
                matrices[(method, level)] = matrix
        w2 = load_matrix(root, key, "W2")
        validate_matrix(w2, expected_n, f"{key}/W2")
        matrices[("W2", None)] = w2

        for method in ("SKOT", "SK_W2DeltaSk"):
            reference = matrices[(method, REFERENCE_L)]
            ref_upper = upper(reference)
            for level in (10, 20, 30):
                matrix = matrices[(method, level)]
                denominator = float(np.linalg.norm(reference))
                convergence_rows.append(
                    {
                        "dataset_key": key,
                        "dataset": label,
                        "method": method,
                        "L": level,
                        "relative_frobenius_error": float(np.linalg.norm(matrix - reference) / denominator) if denominator else 0.0,
                        "spearman_vs_L40": float(spearmanr(upper(matrix), ref_upper).statistic),
                        "NN_at_3_vs_L40": nearest_overlap(matrix, reference, 3),
                    }
                )

        dsk = matrices[("SKOT", PRIMARY_L)]
        wgamma = matrices[("SK_W2DeltaSk", PRIMARY_L)]
        dsk_u = upper(dsk)
        wg_u = upper(wgamma)
        w2_u = upper(w2)
        all_dsk_pairs.append(dsk_u)
        all_wgamma_pairs.append(wg_u)
        all_w2_pairs.append(w2_u)

        pred_dsk_avg = agglomerative_precomputed(dsk, expected_k)
        pred_wg_avg = agglomerative_precomputed(wgamma, expected_k)
        pred_w2_avg = agglomerative_precomputed(w2, expected_k)

        coordinates = hilbert_coordinates(dsk)
        pred_hilbert = KMeans(
            n_clusters=expected_k, n_init=100, random_state=0
        ).fit_predict(coordinates)
        sk_sigma = median_sigma(dsk)
        sk_kernel = gaussian(dsk, sk_sigma)
        pred_sk_spectral = spectral_partition(sk_kernel, expected_k)
        w2_sigma = median_sigma(w2)
        w2_affinity = gaussian(w2, w2_sigma)
        pred_w2_spectral = spectral_partition(w2_affinity, expected_k)

        dsk_time = compute_time(load_info(root, key, "SKOT", PRIMARY_L))
        wg_time = compute_time(load_info(root, key, "SK_W2DeltaSk", PRIMARY_L))
        w2_time = compute_time(load_info(root, key, "W2"))
        mask = w2_u > 1e-15
        ratio_w2_bound = w2_u / np.maximum(math.sqrt(2.0) * dsk_u, 1e-300)
        ratio_wg_bound = wg_u / np.maximum(math.sqrt(2.0) * dsk_u, 1e-300)
        w2_rows.append(
            {
                "dataset_key": key,
                "dataset": label,
                "n": expected_n,
                "median_number_of_pairs": float(np.median(sizes)),
                "spearman_dSK_vs_W2": float(spearmanr(dsk_u, w2_u).statistic),
                "spearman_WGamma_vs_W2": float(spearmanr(wg_u, w2_u).statistic),
                "NN_at_3_dSK_vs_W2": nearest_overlap(dsk, w2, 3),
                "NN_at_3_WGamma_vs_W2": nearest_overlap(wgamma, w2, 3),
                "max_W2_over_sqrt2_dSK": float(np.max(ratio_w2_bound)),
                "max_WGamma_over_sqrt2_dSK": float(np.max(ratio_wg_bound)),
                "W2_above_WGamma_count": int(np.sum(w2_u > wg_u + 1e-9)),
                "WGamma_above_sqrt2_dSK_count": int(np.sum(wg_u > math.sqrt(2.0) * dsk_u + 1e-9)),
                "median_WGamma_over_W2": float(np.median(wg_u[mask] / w2_u[mask])),
                "ARI_partition_dSK_vs_W2": adjusted_rand_score(pred_dsk_avg, pred_w2_avg),
                "ARI_partition_WGamma_vs_W2": adjusted_rand_score(pred_wg_avg, pred_w2_avg),
                "ARI_reference_W2": adjusted_rand_score(labels, pred_w2_avg),
                "ARI_reference_dSK": adjusted_rand_score(labels, pred_dsk_avg),
                "ARI_reference_WGamma": adjusted_rand_score(labels, pred_wg_avg),
                "W2_filter_seconds_mean": w2_time,
                "dSK_filter_seconds_mean": dsk_time,
                "WGamma_filter_seconds_mean": wg_time,
                "speedup_dSK_over_W2": w2_time / dsk_time,
                "speedup_WGamma_over_W2": w2_time / wg_time,
            }
        )

        ari_average = adjusted_rand_score(labels, pred_dsk_avg)
        ari_hilbert = adjusted_rand_score(labels, pred_hilbert)
        ari_spectral = adjusted_rand_score(labels, pred_sk_spectral)
        ari_w2_spectral = adjusted_rand_score(labels, pred_w2_spectral)
        clustering_rows.append(
            {
                "dataset_key": key,
                "dataset": label,
                "n": expected_n,
                "k": expected_k,
                "ARI_average_linkage_dSK": ari_average,
                "ARI_Hilbert_kmeans": ari_hilbert,
                "ARI_Gaussian_SKOT_spectral": ari_spectral,
                "ARI_average_linkage_W2": adjusted_rand_score(labels, pred_w2_avg),
                "ARI_Gaussian_W2_spectral": ari_w2_spectral,
                "sigma_SKOT": sk_sigma,
                "sigma_W2": w2_sigma,
            }
        )
        gaussian_rows.append(
            {
                "dataset_key": key,
                "dataset": label,
                "ARI_partition_Gaussian_SKOT_vs_Gaussian_W2": adjusted_rand_score(pred_sk_spectral, pred_w2_spectral),
                "ARI_reference_Gaussian_SKOT": ari_spectral,
                "ARI_reference_Gaussian_W2": ari_w2_spectral,
            }
        )

        references = true_boundaries(labels)
        if len(references) == expected_k - 1:
            pred_segment, predicted = segment_kernel(sk_kernel, expected_k)
            segmentation_rows.append(
                {
                    "dataset_key": key,
                    "dataset": label,
                    "reference_boundaries": references,
                    "predicted_boundaries": predicted,
                    "boundary_MAE": boundary_mae(references, predicted),
                    "ARI": adjusted_rand_score(labels, pred_segment),
                }
            )

        if key == "2018_asteroid_impact_3D_temporal_subsampling":
            asteroid_coordinates = coordinates
            asteroid_kernel = sk_kernel
            asteroid_labels = labels
            asteroid_reference_bounds = references
            asteroid_predicted_bounds = segment_kernel(sk_kernel, expected_k)[1]

    dataset_df = pd.DataFrame(dataset_rows)
    convergence_df = pd.DataFrame(convergence_rows)
    w2_df = pd.DataFrame(w2_rows)
    clustering_df = pd.DataFrame(clustering_rows)
    gaussian_df = pd.DataFrame(gaussian_rows)
    segmentation_df = pd.DataFrame(segmentation_rows)

    total_diagrams = int(dataset_df.n.sum())
    total_pairs = int(dataset_df.pairs_total.sum())
    total_comparisons = int(sum(n * (n - 1) // 2 for n in dataset_df.n))
    if (total_diagrams, total_pairs, total_comparisons) != (
        EXPECTED_DIAGRAMS,
        EXPECTED_PERSISTENCE_PAIRS,
        EXPECTED_COMPARISONS,
    ):
        raise RuntimeError(
            f"Unexpected benchmark: {(total_diagrams, total_pairs, total_comparisons)}"
        )

    convergence_summary_rows: list[dict[str, Any]] = []
    for method in ("SKOT", "SK_W2DeltaSk"):
        for level in (10, 20, 30):
            subset = convergence_df[
                (convergence_df.method == method) & (convergence_df.L == level)
            ]
            convergence_summary_rows.append(
                {
                    "method": method,
                    "L": level,
                    "median_relative_error": float(subset.relative_frobenius_error.median()),
                    "maximum_relative_error": float(subset.relative_frobenius_error.max()),
                    "minimum_spearman": float(subset.spearman_vs_L40.min()),
                    "minimum_NN_at_3": float(subset.NN_at_3_vs_L40.min()),
                }
            )
    convergence_summary_df = pd.DataFrame(convergence_summary_rows)

    dataset_df.to_csv(csv_dir / "dataset_summary.csv", index=False)
    convergence_df.to_csv(csv_dir / "convergence_all.csv", index=False)
    convergence_summary_df.to_csv(csv_dir / "convergence_summary.csv", index=False)
    w2_df.to_csv(csv_dir / "w2_comparison_L30.csv", index=False)
    clustering_df.to_csv(csv_dir / "clustering_all_L30.csv", index=False)
    gaussian_df.to_csv(csv_dir / "w2_gaussian_affinity_partition_agreement.csv", index=False)
    
    seg_csv = segmentation_df.copy()
    if len(seg_csv):
        seg_csv["reference_boundaries"] = seg_csv.reference_boundaries.map(json.dumps)
        seg_csv["predicted_boundaries"] = seg_csv.predicted_boundaries.map(json.dumps)
    seg_csv.to_csv(csv_dir / "kernel_contiguous_segmentation_supplementary.csv", index=False)

    
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    for method, label in (("SKOT", r"$d_{\mathrm{SK},L}$"), ("SK_W2DeltaSk", r"$W_{\Gamma,L}$")):
        subset = convergence_summary_df[convergence_summary_df.method == method]
        ax.plot(subset.L, subset.median_relative_error, marker="o", label=label + " median")
        ax.plot(subset.L, subset.maximum_relative_error, marker="s", linestyle="--", label=label + " maximum")
    ax.set_yscale("log")
    ax.set_xlabel(r"SK refinement level $L$")
    ax.set_ylabel(r"Relative Frobenius error versus $L=40$")
    ax.set_xticks([10, 20, 30])
    ax.legend()
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(fig_dir / f"fig_convergence_L.{ext}", dpi=220 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)

    all_dsk = np.concatenate(all_dsk_pairs)
    all_wg = np.concatenate(all_wgamma_pairs)
    all_w2 = np.concatenate(all_w2_pairs)
    fig, axes = plt.subplots(1, 2, figsize=(11.3, 4.65))
    axes[0].scatter(all_dsk, all_w2, s=9, alpha=0.28, edgecolors="none")
    maximum = max(float(all_dsk.max()), float(all_w2.max()))
    line = np.linspace(0, maximum, 300)
    axes[0].plot(line, math.sqrt(2.0) * line, linestyle="--", linewidth=1.5)
    axes[0].set_xlabel(r"$d_{\mathrm{SK},30}$")
    axes[0].set_ylabel(r"$W_2$")
    axes[0].text(
        0.03,
        0.96,
        "12 collections, 3,004 diagram pairs\n"
        + rf"median collection-wise Spearman $\rho = {w2_df.spearman_dSK_vs_W2.median():.3f}$",
        transform=axes[0].transAxes,
        va="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85),
    )
    axes[1].scatter(all_wg, all_w2, s=9, alpha=0.28, edgecolors="none")
    maximum = max(float(all_wg.max()), float(all_w2.max()))
    line = np.linspace(0, maximum, 300)
    axes[1].plot(line, line, linestyle="--", linewidth=1.5)
    axes[1].set_xlabel(r"$W_{\Gamma,30}$")
    axes[1].set_ylabel(r"$W_2$")
    axes[1].text(
        0.03,
        0.96,
        "12 collections, 3,004 diagram pairs\n"
        + rf"median collection-wise Spearman $\rho = {w2_df.spearman_WGamma_vs_W2.median():.3f}$",
        transform=axes[1].transAxes,
        va="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85),
    )
    for ax in axes:
        ax.set_xlim(left=-0.02)
        ax.set_ylim(bottom=-0.02)
        ax.grid(alpha=0.15)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(fig_dir / f"fig_w2_scatter_12_collections.{ext}", dpi=220 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)

    ordered = w2_df.sort_values("speedup_dSK_over_W2").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(8.2, 6.2))
    y = np.arange(len(ordered))
    height = 0.36
    ax.barh(y - height / 2, ordered.speedup_dSK_over_W2, height=height, label=r"$d_{\mathrm{SK},30}$")
    ax.barh(y + height / 2, ordered.speedup_WGamma_over_W2, height=height, label=r"$W_{\Gamma,30}$")
    ax.set_xscale("log")
    ax.set_xlabel(r"Speedup over $W_2$")
    ax.set_yticks(y)
    ax.set_yticklabels(ordered.dataset)
    ax.axvline(w2_df.speedup_dSK_over_W2.median(), linestyle=":", linewidth=1.3, label=rf"$d_{{\mathrm{{SK}}}}$ median: {w2_df.speedup_dSK_over_W2.median():.0f}$\times$")
    ax.grid(axis="x", alpha=0.2)
    ax.legend(loc="lower right")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(fig_dir / f"fig_w2_speedup_12_collections.{ext}", dpi=220 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 5.8))
    y = np.arange(len(w2_df))
    height = 0.36
    ax.barh(y - height / 2, w2_df.ARI_partition_dSK_vs_W2, height=height, label=r"$d_{\mathrm{SK},30}$ versus $W_2$")
    ax.barh(y + height / 2, w2_df.ARI_partition_WGamma_vs_W2, height=height, label=r"$W_{\Gamma,30}$ versus $W_2$")
    ax.set_xlim(-0.05, 1.04)
    ax.set_xlabel("Adjusted Rand index between average-linkage partitions")
    ax.set_yticks(y)
    ax.set_yticklabels(w2_df.dataset)
    ax.axvline(1.0, linestyle=":", linewidth=1.2)
    ax.legend(loc="lower right")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(fig_dir / f"fig_w2_partition_agreement_12_collections.{ext}", dpi=220 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10.0, 5.0))
    x = np.arange(len(clustering_df))
    width = 0.25
    ax.bar(x - width, clustering_df.ARI_average_linkage_dSK, width, label="Average linkage on $d_{SK}$")
    ax.bar(x, clustering_df.ARI_Hilbert_kmeans, width, label="Hilbert $k$-means")
    ax.bar(x + width, clustering_df.ARI_Gaussian_SKOT_spectral, width, label="Gaussian SKOT spectral")
    ax.axhline(1.0, linestyle=":")
    ax.set_ylabel("ARI versus benchmark reference partition")
    ax.set_xticks(x)
    ax.set_xticklabels(clustering_df.dataset, rotation=35, ha="right")
    ax.legend()
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(fig_dir / f"fig_clustering_ari.{ext}", dpi=220 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)

    if asteroid_coordinates is None or asteroid_kernel is None or asteroid_labels is None:
        raise RuntimeError("Asteroid temporal case is missing")
    coordinates2 = np.zeros((len(asteroid_coordinates), 2), dtype=float)
    coordinates2[:, : min(2, asteroid_coordinates.shape[1])] = asteroid_coordinates[:, :2]
    fig, ax = plt.subplots(figsize=(7.2, 5.8))
    ax.plot(coordinates2[:, 0], coordinates2[:, 1], linestyle="--", alpha=0.65, label="Temporal order")
    for phase in np.unique(asteroid_labels):
        indices = np.where(asteroid_labels == phase)[0]
        ax.scatter(coordinates2[indices, 0], coordinates2[indices, 1], s=55, label=f"Phase {phase + 1}")
    for index in (0, 4, 5, 9, 10, 14, 15, 19):
        ax.annotate(str(index), (coordinates2[index, 0], coordinates2[index, 1]), xytext=(4, 5), textcoords="offset points", fontsize=8)
    ax.set_xlabel("Hilbert MDS coordinate 1")
    ax.set_ylabel("Hilbert MDS coordinate 2")
    ax.set_title("Asteroid Impact temporal sequence")
    ax.legend(fontsize=8)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(fig_dir / f"fig_asteroid_hilbert_mds.{ext}", dpi=220 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.6, 5.7))
    image = ax.imshow(asteroid_kernel, origin="upper", aspect="equal")
    fig.colorbar(image, ax=ax, label="Kernel similarity")
    first = True
    for boundary in asteroid_reference_bounds:
        ax.axvline(boundary - 0.5, linestyle="-", linewidth=1.4, label="Reference boundary" if first else None)
        ax.axhline(boundary - 0.5, linestyle="-", linewidth=1.4)
        first = False
    first = True
    for boundary in asteroid_predicted_bounds:
        ax.axvline(boundary - 0.5, linestyle="--", linewidth=1.4, label="Kernel segmentation" if first else None)
        ax.axhline(boundary - 0.5, linestyle="--", linewidth=1.4)
        first = False
    ax.set_xlabel("Diagram index")
    ax.set_ylabel("Diagram index")
    ax.set_title("Gaussian SKOT kernel")
    ax.legend(loc="lower left", fontsize=8)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(fig_dir / f"fig_asteroid_skot_kernel.{ext}", dpi=220 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)

    
    dataset_lines = []
    for row in dataset_rows:
        dataset_lines.append(
            f"{row['dataset']} & {row['n']} & {row['k']} & {tex_int(row['pairs_min'])} & {tex_int(row['pairs_median'])} & {tex_int(row['pairs_max'])} \\\\"
        )
    write_text(
        table_dir / "table_dataset_summary.tex",
        r"""\begin{table}[t]
\centering
\caption{
Scientific collections used in the evaluation.
Here, \(n\) denotes the number of persistence diagrams in the
collection, \(k\) the number of reference groups supplied by the
benchmark metadata, and \(|X|\) the number of off-diagonal
persistence pairs in a diagram \(X\). Diagram cardinalities aggregate
all retained critical-pair types and homological dimensions returned
by the Discrete Morse Sandwich backend.
}
\label{tab:dataset_summary}
\scriptsize
\begin{tabular}{lrrrrr}
\toprule
Dataset & \(n\) & \(k\) & Min. \(|X|\) & Median \(|X|\) & Max. \(|X|\) \\
\midrule
""" + "\n".join(dataset_lines) + r"""
\bottomrule
\end{tabular}
\end{table}""",
    )

    conv_lines = []
    for method, tex_method in (("SKOT", r"d_{\mathrm{SK},L}"), ("SK_W2DeltaSk", r"W_{\Gamma,L}")):
        subset = convergence_summary_df[convergence_summary_df.method == method]
        for _, row in subset.iterrows():
            max_err = tex_sci(row.maximum_relative_error) if row.maximum_relative_error < 1 else f"{row.maximum_relative_error:.2f}"
            conv_lines.append(
                rf"\({tex_method}\) & {int(row.L)} & \({tex_sci(row.median_relative_error)}\) & \({max_err}\) & {row.minimum_spearman:.6f} & {row.minimum_NN_at_3:.3f} \\"
            )
        if method == "SKOT":
            conv_lines.append(r"\midrule")
    write_text(
        table_dir / "table_convergence.tex",
        r"""\begin{table}[t]
\centering
\caption{
Convergence to the high-resolution \(L=40\) numerical reference over
the 12 collections. For each collection \(c\),
\(E_c(L)\) is the relative Frobenius error,
\(\rho_c(L)\) is the Spearman correlation between the strict
upper-triangular entries of the \(L\)- and \(L=40\)-level matrices,
and \(\operatorname{NN@3}_c(L)\) is the mean overlap of their
three-nearest-neighbour sets. The table reports the median and
maximum of \(E_c(L)\), and the minimum values of \(\rho_c(L)\) and
\(\operatorname{NN@3}_c(L)\), over the collections. Thus, the last
two columns report the worst collection.
}
\label{tab:selector_convergence}
\small
\begin{tabular}{lccccc}
\toprule
Method & \(L\) & Median \(E_c(L)\) & Max. \(E_c(L)\) & Min. \(\rho_c(L)\) & Min. \(\operatorname{NN@3}_c(L)\) \\
\midrule
""" + "\n".join(conv_lines) + r"""
\bottomrule
\end{tabular}
\end{table}""",
    )

    w2_lines = []
    for _, row in w2_df.iterrows():
        w2_lines.append(
            f"{row.dataset} & {int(row.n)} & {tex_int(row.median_number_of_pairs)} & {row.spearman_dSK_vs_W2:.3f} & {row.spearman_WGamma_vs_W2:.3f} & {row.NN_at_3_dSK_vs_W2:.3f} & {row.NN_at_3_WGamma_vs_W2:.3f} & {row.speedup_dSK_over_W2:.0f}$\\times$ & {row.speedup_WGamma_over_W2:.0f}$\\times$ \\\\"
        )
    w2_lines.append(
        f"Median & -- & -- & {w2_df.spearman_dSK_vs_W2.median():.3f} & {w2_df.spearman_WGamma_vs_W2.median():.3f} & {w2_df.NN_at_3_dSK_vs_W2.median():.3f} & {w2_df.NN_at_3_WGamma_vs_W2.median():.3f} & {w2_df.speedup_dSK_over_W2.median():.0f}$\\times$ & {w2_df.speedup_WGamma_over_W2.median():.0f}$\\times$ \\\\"
    )
    write_text(
        table_dir / "table_w2_comparison.tex",
        r"""\begin{table}[t]
\centering
\caption{Faithfulness and computational performance at \(L=30\) on
the 12-collection benchmark (3,004 distinct diagram pairs). \(\rho\)
is Spearman rank correlation, NN@3 is the mean overlap of the
three-nearest-neighbour sets, and speedups use the mean filter-compute
time recorded by the reproduction run(s).}
\label{tab:w2_comparison}
\scriptsize
\setlength{\tabcolsep}{2.6pt}
\resizebox{\linewidth}{!}{%
\begin{tabular}{lrrrrrrrr}
\toprule
Dataset & \(n\) & Med. \(|X|\) & \(\rho(d_{\mathrm{SK}},W_2)\) & \(\rho(W_\Gamma,W_2)\) & NN@3 \(d_{\mathrm{SK}}\) & NN@3 \(W_\Gamma\) & Speedup \(d_{\mathrm{SK}}\) & Speedup \(W_\Gamma\) \\
\midrule
""" + "\n".join(w2_lines[:-1]) + "\n\\midrule\n" + w2_lines[-1] + r"""
\bottomrule
\end{tabular}}
\end{table}""",
    )

    partition_lines = []
    for _, row in w2_df.iterrows():
        partition_lines.append(
            f"{row.dataset} & {row.ARI_partition_dSK_vs_W2:.3f} & {row.ARI_partition_WGamma_vs_W2:.3f} & {row.ARI_reference_W2:.3f} & {row.ARI_reference_dSK:.3f} & {row.ARI_reference_WGamma:.3f} \\\\"
        )
    partition_lines.append(
        f"Mean & {w2_df.ARI_partition_dSK_vs_W2.mean():.3f} & {w2_df.ARI_partition_WGamma_vs_W2.mean():.3f} & {w2_df.ARI_reference_W2.mean():.3f} & {w2_df.ARI_reference_dSK.mean():.3f} & {w2_df.ARI_reference_WGamma.mean():.3f} \\\\"
    )
    write_text(
        table_dir / "table_w2_partition_agreement_supplementary.tex",
        r"""\begin{table}[t]
\centering
\caption{
Agreement between average-linkage partitions on the complete
12-collection benchmark. The first two ARI columns compare the
partitions obtained from the SK constructions with the partition
obtained from \(W_2\). The last three columns compare each partition
with the reference groups supplied by the benchmark metadata. The
number of groups \(k\) is fixed from these metadata.
}
\label{tab:w2_partition_agreement}
\scriptsize
\begin{tabular}{lrrrrr}
\toprule
Dataset & ARI \(d_{\mathrm{SK}}/W_2\) & ARI \(W_\Gamma/W_2\) & ARI ref./\(W_2\) & ARI ref./\(d_{\mathrm{SK}}\) & ARI ref./\(W_\Gamma\) \\
\midrule
""" + "\n".join(partition_lines[:-1]) + "\n\\midrule\n" + partition_lines[-1] + r"""
\bottomrule
\end{tabular}
\end{table}""",
    )

    cluster_lines = []
    for _, row in clustering_df.iterrows():
        values = [row.ARI_average_linkage_dSK, row.ARI_Hilbert_kmeans, row.ARI_Gaussian_SKOT_spectral]
        cluster_lines.append(
            f"{row.dataset} & {int(row.n)} & {int(row.k)} & {bold_best(values, 0)} & {bold_best(values, 1)} & {bold_best(values, 2)} \\\\"
        )
    means = [clustering_df.ARI_average_linkage_dSK.mean(), clustering_df.ARI_Hilbert_kmeans.mean(), clustering_df.ARI_Gaussian_SKOT_spectral.mean()]
    cluster_lines.append(f"Mean & -- & -- & {bold_best(means,0)} & {bold_best(means,1)} & {bold_best(means,2)} \\\\ ")
    write_text(
        table_dir / "table_kernel_clustering.tex",
        r"""\begin{table}[t]
\centering
\caption{Adjusted Rand index (ARI) at \(L=30\), measured with respect to the
reference partition supplied by the dataset metadata. The number of
clusters \(k\) is fixed from these metadata, while the reference
labels are not used to compute the cluster assignments.
``Hilbert \(k\)-means'' uses the full-dimensional Euclidean
realization of \(d_{\mathrm{SK},30}\), and
``Gaussian SKOT spectral'' uses
\(k_\sigma(X,Y)=\exp[-d_{\mathrm{SK},30}(X,Y)^2/(2\sigma^2)]\)
with the median-distance heuristic. The ``Average linkage'' column
reproduces the ``ARI ref./\(d_{\mathrm{SK}}\)'' column of
\autoref{tab:w2_partition_agreement}; it is repeated here as the
metric-space baseline for comparison with Hilbert \(k\)-means and
Gaussian SKOT spectral clustering.}
\label{tab:kernel_clustering}
\scriptsize
\setlength{\tabcolsep}{4pt}
\begin{tabular}{lrrrrr}
\toprule
Dataset & \(n\) & \(k\) & Average linkage & Hilbert \(k\)-means & Gaussian SKOT spectral \\
\midrule
""" + "\n".join(cluster_lines[:-1]) + "\n\\midrule\n" + cluster_lines[-1] + r"""
\bottomrule
\end{tabular}
\end{table}""",
    )

    gaussian_lines = []
    for _, row in gaussian_df.iterrows():
        gaussian_lines.append(
            f"{row.dataset} & {row.ARI_partition_Gaussian_SKOT_vs_Gaussian_W2:.3f} & {row.ARI_reference_Gaussian_SKOT:.3f} & {row.ARI_reference_Gaussian_W2:.3f} \\\\"
        )
    gaussian_lines.append(
        f"Mean & {gaussian_df.ARI_partition_Gaussian_SKOT_vs_Gaussian_W2.mean():.3f} & {gaussian_df.ARI_reference_Gaussian_SKOT.mean():.3f} & {gaussian_df.ARI_reference_Gaussian_W2.mean():.3f} \\\\"
    )
    write_text(
        table_dir / "table_w2_gaussian_agreement_supplementary.tex",
        r"""\begin{table}[t]
\centering
\caption{Agreement between spectral partitions obtained from Gaussian
SKOT and analogous Gaussian \(W_2\) affinities. Each bandwidth is the
median positive distance of its own matrix.}
\label{tab:w2_gaussian_agreement_supp}
\scriptsize
\begin{tabular}{lrrr}
\toprule
Dataset & ARI between partitions & ARI ref./SKOT & ARI ref./\(W_2\) \\
\midrule
""" + "\n".join(gaussian_lines[:-1]) + "\n\\midrule\n" + gaussian_lines[-1] + r"""
\bottomrule
\end{tabular}
\end{table}""",
    )

    segmentation_lines = []
    for _, row in segmentation_df.iterrows():
        segmentation_lines.append(
            f"{row.dataset} & {format_boundaries(row.reference_boundaries)} & {format_boundaries(row.predicted_boundaries)} & {row.boundary_MAE:.2f} & {row.ARI:.3f} \\\\"
        )
    write_text(
        table_dir / "table_kernel_segmentation_supplementary.tex",
        r"""\begin{table}[t]
\centering
\caption{
Fixed-\(k\) contiguous segmentation using the Gaussian SKOT kernel at
\(L=30\). For each collection, \(k\) is the number of reference groups
supplied by the metadata, and the Gaussian bandwidth is the median of
the positive pairwise \(d_{\mathrm{SK},30}\) distances, without
label-dependent tuning. Only the 11 collections whose reference
labels form \(k\) contiguous blocks in the stored sample order are
included. A boundary \(b\) denotes the first index of the next
half-open segment; for example, boundaries \([4,8]\) define
\([0,4)\), \([4,8)\), and \([8,n)\). Boundary MAE is the mean absolute
difference, measured in sample indices, between corresponding
reference and predicted boundaries. ARI compares the complete
predicted segmentation with the metadata reference partition. The
stored order is explicitly temporal for some collections and may
represent another ordered parameterization for others.
}
\label{tab:kernel_segmentation_supp}
\scriptsize
\begin{tabular}{lrrrr}
\toprule
Dataset & Reference boundaries & Predicted boundaries & Boundary MAE & ARI vs. reference \\
\midrule
""" + "\n".join(segmentation_lines) + r"""
\bottomrule
\end{tabular}
\end{table}""",
    )

    
    exact_avg_dsk = int(np.isclose(w2_df.ARI_partition_dSK_vs_W2, 1.0).sum())
    exact_avg_wg = int(np.isclose(w2_df.ARI_partition_WGamma_vs_W2, 1.0).sum())
    exact_spectral = int(np.isclose(gaussian_df.ARI_partition_Gaussian_SKOT_vs_Gaussian_W2, 1.0).sum())
    exact_segmentations = int(np.isclose(segmentation_df.ARI, 1.0).sum())
    summary = {
        "benchmark": {
            "collections": 12,
            "diagrams": total_diagrams,
            "persistence_pairs": total_pairs,
            "pairwise_comparisons_per_method": total_comparisons,
        },
        "convergence": {
            "dSK_L30_median_relative_error": float(convergence_summary_df.query("method == 'SKOT' and L == 30").median_relative_error.iloc[0]),
            "dSK_L30_max_relative_error": float(convergence_summary_df.query("method == 'SKOT' and L == 30").maximum_relative_error.iloc[0]),
            "WGamma_L30_median_relative_error": float(convergence_summary_df.query("method == 'SK_W2DeltaSk' and L == 30").median_relative_error.iloc[0]),
            "WGamma_L30_max_relative_error": float(convergence_summary_df.query("method == 'SK_W2DeltaSk' and L == 30").maximum_relative_error.iloc[0]),
        },
        "faithfulness_W2": {
            "median_spearman_dSK": float(w2_df.spearman_dSK_vs_W2.median()),
            "median_spearman_WGamma": float(w2_df.spearman_WGamma_vs_W2.median()),
            "median_NN3_dSK": float(w2_df.NN_at_3_dSK_vs_W2.median()),
            "median_NN3_WGamma": float(w2_df.NN_at_3_WGamma_vs_W2.median()),
            "max_W2_over_sqrt2_dSK": float(w2_df.max_W2_over_sqrt2_dSK.max()),
            "max_WGamma_over_sqrt2_dSK": float(w2_df.max_WGamma_over_sqrt2_dSK.max()),
            "bound_violations": int(w2_df.W2_above_WGamma_count.sum() + w2_df.WGamma_above_sqrt2_dSK_count.sum()),
            "average_linkage_exact_dSK_vs_W2": exact_avg_dsk,
            "average_linkage_exact_WGamma_vs_W2": exact_avg_wg,
            "mean_ARI_partition_dSK_vs_W2": float(w2_df.ARI_partition_dSK_vs_W2.mean()),
            "mean_ARI_partition_WGamma_vs_W2": float(w2_df.ARI_partition_WGamma_vs_W2.mean()),
        },
        "timing": {
            "uses_mean_of_recorded_runs": True,
            "median_speedup_dSK": float(w2_df.speedup_dSK_over_W2.median()),
            "median_speedup_WGamma": float(w2_df.speedup_WGamma_over_W2.median()),
            "minimum_speedup_dSK": float(w2_df.speedup_dSK_over_W2.min()),
            "maximum_speedup_dSK": float(w2_df.speedup_dSK_over_W2.max()),
            "sum_W2_seconds": float(w2_df.W2_filter_seconds_mean.sum()),
            "sum_dSK_seconds": float(w2_df.dSK_filter_seconds_mean.sum()),
        },
        "clustering_reference": {
            "mean_ARI_average_linkage_dSK": float(clustering_df.ARI_average_linkage_dSK.mean()),
            "mean_ARI_Hilbert_kmeans": float(clustering_df.ARI_Hilbert_kmeans.mean()),
            "mean_ARI_Gaussian_SKOT_spectral": float(clustering_df.ARI_Gaussian_SKOT_spectral.mean()),
            "mean_ARI_average_linkage_W2": float(clustering_df.ARI_average_linkage_W2.mean()),
            "perfect_Gaussian_SKOT_spectral": int(np.isclose(clustering_df.ARI_Gaussian_SKOT_spectral, 1.0).sum()),
            "improved_vs_average_linkage": int((clustering_df.ARI_Gaussian_SKOT_spectral > clustering_df.ARI_average_linkage_dSK + 1e-12).sum()),
            "worsened_vs_average_linkage": int((clustering_df.ARI_Gaussian_SKOT_spectral < clustering_df.ARI_average_linkage_dSK - 1e-12).sum()),
        },
        "Gaussian_SKOT_vs_W2_spectral": {
            "exact_partitions": exact_spectral,
            "median_ARI_between_partitions": float(gaussian_df.ARI_partition_Gaussian_SKOT_vs_Gaussian_W2.median()),
            "mean_ARI_between_partitions": float(gaussian_df.ARI_partition_Gaussian_SKOT_vs_Gaussian_W2.mean()),
        },
        "contiguous_segmentation": {
            "eligible_collections": int(len(segmentation_df)),
            "exact_segmentations": exact_segmentations,
        },
    }
    write_text(output / "SUMMARY.json", json.dumps(summary, indent=2, ensure_ascii=False))

    expected_checks = {
        "median_spearman_dSK": (summary["faithfulness_W2"]["median_spearman_dSK"], 0.878828823, 2e-6),
        "median_spearman_WGamma": (summary["faithfulness_W2"]["median_spearman_WGamma"], 0.923729164, 2e-6),
        "mean_ARI_average": (summary["clustering_reference"]["mean_ARI_average_linkage_dSK"], 0.666882642, 2e-6),
        "mean_ARI_Hilbert": (summary["clustering_reference"]["mean_ARI_Hilbert_kmeans"], 0.756430422, 2e-6),
        "mean_ARI_spectral": (summary["clustering_reference"]["mean_ARI_Gaussian_SKOT_spectral"], 0.799863463, 2e-6),
        "mean_ARI_spectral_partition_agreement": (summary["Gaussian_SKOT_vs_W2_spectral"]["mean_ARI_between_partitions"], 0.936, 5e-4),
    }
    validation_lines = ["VALIDATION OF NON-TIMING RESULTS"]
    for name, (actual, expected, tolerance) in expected_checks.items():
        ok = abs(actual - expected) <= tolerance
        validation_lines.append(f"{name}: actual={actual:.12g}, expected={expected:.12g}, tolerance={tolerance:g}, ok={ok}")
        if not ok:
            raise RuntimeError(f"Validation failed for {name}")
    validation_lines += [
        f"bound_violations: {summary['faithfulness_W2']['bound_violations']} (expected 0)",
        f"average_linkage_exact_dSK_vs_W2: {exact_avg_dsk} (expected 8)",
        f"Gaussian_SKOT_vs_W2_exact: {exact_spectral} (expected 10)",
        f"exact_segmentations: {exact_segmentations} (expected 6)",
    ]
    if summary["faithfulness_W2"]["bound_violations"] != 0 or exact_avg_dsk != 8 or exact_spectral != 10 or exact_segmentations != 6:
        raise RuntimeError("Structural validation failed")
    write_text(output / "RESULT_VALIDATION.txt", "\n".join(validation_lines))

    key_results = f"""REPRODUCED KEY RESULTS

Benchmark: 12 collections, {total_diagrams} diagrams, {tex_int(total_pairs)} persistence pairs, {total_comparisons} comparisons per method.
Median Spearman with W2: d_SK={w2_df.spearman_dSK_vs_W2.median():.3f}, W_Gamma={w2_df.spearman_WGamma_vs_W2.median():.3f}.
Median NN@3 with W2: d_SK={w2_df.NN_at_3_dSK_vs_W2.median():.3f}, W_Gamma={w2_df.NN_at_3_WGamma_vs_W2.median():.3f}.
Violations of W2 <= W_Gamma <= sqrt(2)d_SK: 0.
Average linkage identical to W2: {exact_avg_dsk}/12 for d_SK.
Mean ARI versus the reference partition: average-linkage d_SK={clustering_df.ARI_average_linkage_dSK.mean():.3f}, Hilbert k-means={clustering_df.ARI_Hilbert_kmeans.mean():.3f}, Gaussian SKOT spectral={clustering_df.ARI_Gaussian_SKOT_spectral.mean():.3f}.
Identical SKOT/W2 spectral partitions: {exact_spectral}/12.
Exact contiguous segmentations: {exact_segmentations}/{len(segmentation_df)} eligible collections.

Speedups are recomputed from the runtimes stored in the info.json files; they therefore depend on the machine and may differ from the paper values.
"""
    write_text(output / "KEY_RESULTS.txt", key_results)

    
    if args.paper_root is not None:
        paper = args.paper_root.expanduser().resolve()
        paper_figures = paper / "figs" / "skot_experiments"
        paper_tables = paper / "tables"
        paper_figures.mkdir(parents=True, exist_ok=True)
        paper_tables.mkdir(parents=True, exist_ok=True)
        for path in fig_dir.iterdir():
            if path.is_file():
                shutil.copy2(path, paper_figures / path.name)
        for path in table_dir.iterdir():
            if path.is_file():
                shutil.copy2(path, paper_tables / path.name)
        write_text(paper / "REPRODUCTION_RESULTS.txt", key_results)

    checksums: list[str] = []
    checksum_path = output / "SHA256SUMS.txt"
    for path in sorted(output.rglob("*")):
        if path.is_file() and path != checksum_path:
            checksums.append(f"{sha256(path)}  {path.relative_to(output)}")
    write_text(checksum_path, "\n".join(checksums))

    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print("REPRODUCTION_OK=1")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
