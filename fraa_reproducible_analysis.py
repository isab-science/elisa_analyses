#!/usr/bin/env python3
"""
Compare paired ELISA curve exports for an antigen-coated plate versus a blocked
control plate.

Usage
-----
python fraa_reproducible_analysis.py --antigen-input curves_APOER2.xlsx --blocked-input curves_Blocked.xlsx --antigen-name ApoER2
python fraa_reproducible_analysis.py  # opens file pickers
"""
from __future__ import annotations

import argparse
import math
import re
import warnings
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

try:
    import tkinter as tk
    from tkinter import filedialog
except ImportError:
    tk = None
    filedialog = None

import matplotlib.pyplot as plt
import numpy as np
import openpyxl
import pandas as pd
from openpyxl.utils.cell import column_index_from_string
from scipy.optimize import curve_fit


def sanitize_slug(text: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")
    return slug or "antigen"


def clean_text(value: object) -> str:
    return str(value).replace("\xa0", " ").strip()


def normalize_header(value: object) -> str:
    text = clean_text(value).lower()
    return re.sub(r"[^a-z0-9]+", "", text)


def parse_numeric_value(value: object) -> float:
    if value is None:
        return np.nan
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return float(value)

    text = clean_text(value)
    if not text:
        return np.nan

    numbers: List[float] = []
    for part in text.split(";"):
        token = part.strip()
        if not token or token.upper() == "NA":
            continue
        try:
            numbers.append(float(token))
        except ValueError:
            continue

    if not numbers:
        return np.nan
    return float(np.mean(numbers))


def parse_column_range(spec: str) -> List[int]:
    token = spec.replace(" ", "").upper()
    if ":" in token:
        start, end = token.split(":", 1)
        start_idx = column_index_from_string(start)
        end_idx = column_index_from_string(end)
        if end_idx < start_idx:
            raise SystemExit(f"Invalid measurement column range: {spec}")
        return list(range(start_idx, end_idx + 1))
    return [column_index_from_string(part) for part in token.split(",") if part]


def choose_input_workbook(title: str) -> Path:
    if tk is None or filedialog is None:
        raise SystemExit("A file picker was requested, but tkinter is unavailable.")

    root = tk.Tk()
    root.withdraw()
    root.attributes("-topmost", True)
    try:
        selected = filedialog.askopenfilename(
            title=title,
            filetypes=[
                ("Excel workbooks", "*.xlsx *.xlsm *.xls"),
                ("All files", "*.*"),
            ],
        )
    finally:
        root.destroy()

    if not selected:
        raise SystemExit(f"No file selected for: {title}")
    return Path(selected)


def choose_antigen_name(path: Path) -> str:
    text = path.stem
    text = re.sub(r"^curves[_\-\s]*", "", text, flags=re.IGNORECASE)
    text = re.sub(r"[_\-]+", " ", text).strip()
    return text or "Antigen"


def build_sample_metadata(value: object) -> Dict[str, object] | None:
    if value is None:
        return None
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        patient_id = int(value)
        return {
            "sample_key": f"patient:{patient_id}",
            "sample_id": f"Patient {patient_id}",
            "patient_id": patient_id,
            "group": "Patient",
        }

    text = clean_text(value)
    if not text:
        return None
    if re.fullmatch(r"\d+", text):
        patient_id = int(text)
        return {
            "sample_key": f"patient:{patient_id}",
            "sample_id": f"Patient {patient_id}",
            "patient_id": patient_id,
            "group": "Patient",
        }
    if "control" in text.lower():
        key = re.sub(r"\s+", " ", text.lower())
        return {
            "sample_key": f"control:{key}",
            "sample_id": text,
            "patient_id": np.nan,
            "group": "Control",
        }
    return None


def find_header_column(worksheet: openpyxl.worksheet.worksheet.Worksheet, patterns: Sequence[str]) -> int:
    normalized_patterns = [re.sub(r"[^a-z0-9]+", "", p.lower()) for p in patterns]
    headers = [worksheet.cell(1, col).value for col in range(1, worksheet.max_column + 1)]
    for col_idx, header in enumerate(headers, start=1):
        normalized = normalize_header(header)
        if any(pattern in normalized for pattern in normalized_patterns):
            return col_idx
    raise SystemExit(f"Could not find a workbook column matching {patterns}.")


def dilution_factors_from_headers(
    worksheet: openpyxl.worksheet.worksheet.Worksheet, measurement_cols: Sequence[int]
) -> np.ndarray:
    dilution_factors: List[float] = []
    for col_idx in measurement_cols:
        header_value = parse_numeric_value(worksheet.cell(1, col_idx).value)
        if not np.isfinite(header_value) or header_value <= 0:
            raise SystemExit(
                f"Measurement column {col_idx} does not have a numeric dilution header in row 1."
            )
        dilution_factors.append(1.0 / header_value)
    return np.asarray(dilution_factors, dtype=float)


def ec50_from_log10_column(value: object) -> float:
    log10_ec50 = parse_numeric_value(value)
    if not np.isfinite(log10_ec50):
        return np.nan
    try:
        return float(10 ** abs(log10_ec50))
    except OverflowError:
        return np.nan


def interpolate_logx_at_y(x1: float, y1: float, x2: float, y2: float, target_y: float) -> float:
    if not np.isfinite(y1) or not np.isfinite(y2):
        return x1
    if abs(y2 - y1) < 1e-12:
        return (x1 + x2) / 2.0
    weight = (target_y - y1) / (y2 - y1)
    return x1 + weight * (x2 - x1)


def fallback_ec50_from_measurements(dilution_factors: np.ndarray, y_values: Iterable[float]) -> float:
    y = np.asarray(list(y_values), dtype=float)
    mask = np.isfinite(y)
    if mask.sum() < 2:
        return np.nan

    x = np.log10(dilution_factors[mask])
    y = y[mask]
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    y_monotone = np.minimum.accumulate(y)
    target = float(np.nanmin(y_monotone) + 0.5 * (np.nanmax(y_monotone) - np.nanmin(y_monotone)))

    for idx in range(len(x) - 1):
        high = y_monotone[idx]
        low = y_monotone[idx + 1]
        if high >= target >= low:
            return float(10 ** interpolate_logx_at_y(x[idx], high, x[idx + 1], low, target))

    if target > y_monotone[0]:
        return float(10 ** interpolate_logx_at_y(x[0], y_monotone[0], x[1], y_monotone[1], target))
    return float(10 ** interpolate_logx_at_y(x[-2], y_monotone[-2], x[-1], y_monotone[-1], target))


def load_condition_workbook(
    path: Path,
    measurement_cols: Sequence[int],
    condition_prefix: str,
) -> Tuple[pd.DataFrame, np.ndarray]:
    workbook = openpyxl.load_workbook(path, data_only=True)
    worksheet = workbook[workbook.sheetnames[0]]
    patient_col = find_header_column(worksheet, ["patient_id", "patient id"])
    log10_ec50_col = find_header_column(worksheet, ["log10(ec50)", "log10ec50"])
    dilution_factors = dilution_factors_from_headers(worksheet, measurement_cols)

    measurement_names = [f"{condition_prefix}_d{i}" for i in range(1, len(measurement_cols) + 1)]
    records: List[Dict[str, object]] = []
    for row_idx in range(2, worksheet.max_row + 1):
        sample_meta = build_sample_metadata(worksheet.cell(row_idx, patient_col).value)
        if sample_meta is None:
            continue

        measurements = [parse_numeric_value(worksheet.cell(row_idx, col_idx).value) for col_idx in measurement_cols]
        if np.isfinite(measurements).sum() < 2:
            continue

        ec50 = ec50_from_log10_column(worksheet.cell(row_idx, log10_ec50_col).value)
        if (
            not np.isfinite(ec50)
            or ec50 <= 0
            or ec50 > dilution_factors.max() * 1e3
            or ec50 < dilution_factors.min() * 1e-3
        ):
            ec50 = fallback_ec50_from_measurements(dilution_factors, measurements)

        record = dict(sample_meta)
        record[f"{condition_prefix}_ec50"] = ec50
        for col_name, value in zip(measurement_names, measurements):
            record[col_name] = value
        records.append(record)

    if not records:
        raise SystemExit(f"No usable rows were found in {path}.")

    df = pd.DataFrame(records)
    aggregate_cols = [f"{condition_prefix}_ec50", *measurement_names]
    aggregated = (
        df.groupby(["sample_key", "sample_id", "patient_id", "group"], dropna=False)[aggregate_cols]
        .mean()
        .reset_index()
    )
    return aggregated, dilution_factors


def merge_condition_data(antigen_df: pd.DataFrame, blocked_df: pd.DataFrame) -> pd.DataFrame:
    merged = antigen_df.merge(blocked_df, on="sample_key", how="inner", suffixes=("_antigen", "_blocked"))
    if merged.empty:
        raise SystemExit("No overlapping samples were found between the antigen and blocked workbooks.")

    for field in ["sample_id", "patient_id", "group"]:
        left = f"{field}_antigen"
        right = f"{field}_blocked"
        mismatch = merged[left].notna() & merged[right].notna() & (merged[left] != merged[right])
        if mismatch.any():
            raise SystemExit(f"Sample metadata mismatch detected for field '{field}'.")
        merged[field] = merged[left].combine_first(merged[right])
        merged = merged.drop(columns=[left, right])

    return merged


def fourpl(logx: np.ndarray, bottom: float, top: float, hill: float, logec50: float) -> np.ndarray:
    return bottom + (top - bottom) / (1.0 + 10.0 ** ((logx - logec50) * hill))


def dilution_values(row: pd.Series, prefix: str, n_dilutions: int) -> np.ndarray:
    return row[[f"{prefix}_d{i}" for i in range(1, n_dilutions + 1)]].astype(float).to_numpy()

def trapezoid_with_nans(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 2:
        return np.nan
    return float(np.trapezoid(y[mask], x[mask]))


def approximate_logec50_se(
    dilution_factors: np.ndarray,
    y_values: Iterable[float],
    ec50_reported: float,
) -> Dict[str, float]:
    y = np.asarray(list(y_values), dtype=float)
    mask = np.isfinite(y)
    if mask.sum() < 4 or not np.isfinite(ec50_reported) or ec50_reported <= 0:
        return {
            "bottom": np.nan,
            "top": np.nan,
            "hill": np.nan,
            "rmse": np.nan,
            "se_log10_ec50": np.nan,
        }

    logx = np.log10(dilution_factors[mask])
    y = y[mask]
    logec50 = float(np.log10(ec50_reported))
    ymin = float(np.min(y))
    ymax = float(np.max(y))
    yrange = max(ymax - ymin, 1e-9)

    def model(local_logx: np.ndarray, bottom: float, top: float, hill: float) -> np.ndarray:
        return fourpl(local_logx, bottom, top, hill, logec50)

    try:
        popt, _ = curve_fit(
            model,
            logx,
            y,
            p0=[ymin, ymax, 1.0],
            bounds=([ymin - yrange, ymin, 0.01], [ymax, ymax + yrange, 20.0]),
            maxfev=100000,
        )
    except Exception:
        return {
            "bottom": np.nan,
            "top": np.nan,
            "hill": np.nan,
            "rmse": np.nan,
            "se_log10_ec50": np.nan,
        }

    yhat = model(logx, *popt)
    residuals = y - yhat
    df_resid = len(y) - len(popt)
    rmse = float(np.sqrt(np.sum(residuals ** 2) / df_resid)) if df_resid > 0 else np.nan
    bottom, top, hill = [float(v) for v in popt]
    slope = abs((top - bottom) * hill * np.log(10.0) / 4.0)
    se_log10_ec50 = rmse / slope if np.isfinite(rmse) and slope > 0 else np.nan

    return {
        "bottom": bottom,
        "top": top,
        "hill": hill,
        "rmse": rmse,
        "se_log10_ec50": se_log10_ec50,
    }


def make_errorbars(center: float, se_log10: float) -> Tuple[float, float]:
    if not np.isfinite(center) or center <= 0 or not np.isfinite(se_log10) or se_log10 < 0:
        return 0.0, 0.0
    lower = center - 10 ** (np.log10(center) - se_log10)
    upper = 10 ** (np.log10(center) + se_log10) - center
    return float(lower), float(upper)


def compute_derived_metrics(df: pd.DataFrame, dilution_factors: np.ndarray) -> pd.DataFrame:
    out = df.copy()
    n_dilutions = len(dilution_factors)
    log_dilutions = np.log10(dilution_factors)

    antigen_max: List[float] = []
    blocked_max: List[float] = []
    antigen_auc: List[float] = []
    blocked_auc: List[float] = []
    net_auc: List[float] = []
    ratio_by_dilution: List[np.ndarray] = []
    antigen_se: List[float] = []
    blocked_se: List[float] = []
    antigen_rmse: List[float] = []
    blocked_rmse: List[float] = []
    xerr_low: List[float] = []
    xerr_high: List[float] = []
    yerr_low: List[float] = []
    yerr_high: List[float] = []

    for _, row in out.iterrows():
        antigen_y = dilution_values(row, "antigen", n_dilutions)
        blocked_y = dilution_values(row, "blocked", n_dilutions)

        antigen_max.append(float(np.nanmax(antigen_y)))
        blocked_max.append(float(np.nanmax(blocked_y)))
        antigen_auc.append(trapezoid_with_nans(log_dilutions, antigen_y))
        blocked_auc.append(trapezoid_with_nans(log_dilutions, blocked_y))
        net_auc.append(trapezoid_with_nans(log_dilutions, antigen_y - blocked_y))

        with np.errstate(divide="ignore", invalid="ignore"):
            ratios = antigen_y / blocked_y
        ratio_by_dilution.append(ratios)

        antigen_stats = approximate_logec50_se(dilution_factors, antigen_y, row["antigen_ec50"])
        blocked_stats = approximate_logec50_se(dilution_factors, blocked_y, row["blocked_ec50"])

        antigen_se.append(antigen_stats["se_log10_ec50"])
        blocked_se.append(blocked_stats["se_log10_ec50"])
        antigen_rmse.append(antigen_stats["rmse"])
        blocked_rmse.append(blocked_stats["rmse"])

        xl, xh = make_errorbars(row["blocked_ec50"], blocked_stats["se_log10_ec50"])
        yl, yh = make_errorbars(row["antigen_ec50"], antigen_stats["se_log10_ec50"])
        xerr_low.append(xl)
        xerr_high.append(xh)
        yerr_low.append(yl)
        yerr_high.append(yh)

    out["antigen_max"] = antigen_max
    out["blocked_max"] = blocked_max
    out["antigen_auc"] = antigen_auc
    out["blocked_auc"] = blocked_auc
    out["net_auc"] = net_auc
    out["antigen_se_log10_ec50"] = antigen_se
    out["blocked_se_log10_ec50"] = blocked_se
    out["antigen_rmse"] = antigen_rmse
    out["blocked_rmse"] = blocked_rmse
    out["xerr_low"] = xerr_low
    out["xerr_high"] = xerr_high
    out["yerr_low"] = yerr_low
    out["yerr_high"] = yerr_high
    out["grey_flag"] = (out["group"] == "Patient") & (out["blocked_ec50"] > out["antigen_ec50"] / 2.0)
    out["blocked_gt_half_antigen"] = np.where(out["grey_flag"], "Yes", "No")
    out["log10_ratio_antigen_blocked_ec50"] = np.log10(out["antigen_ec50"] / out["blocked_ec50"])

    for idx in range(n_dilutions):
        out[f"antigen_over_blocked_ratio_d{idx + 1}"] = [
            float(r[idx]) if np.isfinite(r[idx]) else np.nan for r in ratio_by_dilution
        ]

    return out


def save_basic_data_tables(df: pd.DataFrame, outdir: Path, antigen_name: str) -> None:
    df.to_csv(outdir / "parsed_ec50_data.csv", index=False)
    patient_table = (
        df[df["group"] == "Patient"][["patient_id", "antigen_ec50", "blocked_ec50", "blocked_gt_half_antigen"]]
        .rename(
            columns={
                "patient_id": "Patient ID",
                "antigen_ec50": f"EC50 {antigen_name}",
                "blocked_ec50": "EC50 blocked",
                "blocked_gt_half_antigen": f"Blocked > {antigen_name}/2",
            }
        )
        .sort_values("Patient ID")
        .copy()
    )
    patient_table[f"EC50 {antigen_name}"] = patient_table[f"EC50 {antigen_name}"].round(4)
    patient_table["EC50 blocked"] = patient_table["EC50 blocked"].round(4)
    patient_table.to_csv(outdir / "patient_ec50_table.csv", index=False)
    with pd.ExcelWriter(outdir / "patient_ec50_table.xlsx", engine="openpyxl") as writer:
        patient_table.to_excel(writer, sheet_name="Patient EC50 Table", index=False)


def log_axis_limits(df: pd.DataFrame, dilution_factors: np.ndarray) -> Tuple[float, float]:
    values = np.concatenate(
        [
            dilution_factors,
            df["antigen_ec50"].astype(float).to_numpy(),
            df["blocked_ec50"].astype(float).to_numpy(),
        ]
    )
    values = values[np.isfinite(values) & (values > 0)]
    return 10 ** math.floor(math.log10(float(values.min())) - 0.2), 10 ** math.ceil(math.log10(float(values.max())) + 0.2)


def configure_log_ec50_axes(ax: plt.Axes, df: pd.DataFrame, dilution_factors: np.ndarray) -> None:
    lower, upper = log_axis_limits(df, dilution_factors)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lower, upper)
    ax.set_ylim(lower, upper)
    ax.grid(True, alpha=0.25)


def make_log_histogram_bins(x_values: np.ndarray, dilution_factors: np.ndarray) -> np.ndarray:
    base = np.sort(np.unique(dilution_factors))
    if len(base) < 2:
        return np.logspace(np.log10(x_values.min() / 2.0), np.log10(x_values.max() * 2.0), 8)

    ratios = base[1:] / base[:-1]
    typical_ratio = float(np.median(ratios))
    bins = [base[0] / math.sqrt(ratios[0])]
    bins.extend(np.sqrt(base[:-1] * base[1:]))
    bins.append(base[-1] * math.sqrt(ratios[-1]))
    while bins[0] > x_values.min() / 1.2:
        bins.insert(0, bins[0] / typical_ratio)
    while bins[-1] < x_values.max() * 1.2:
        bins.append(bins[-1] * typical_ratio)
    return np.asarray(bins, dtype=float)

def plot_figure_1_histogram(
    df: pd.DataFrame,
    outdir: Path,
    antigen_name: str,
    antigen_slug: str,
    dilution_factors: np.ndarray,
) -> None:
    patients = df[df["group"] == "Patient"].copy()
    controls = df[df["group"] == "Control"].copy()
    if patients.empty:
        return

    x = patients["antigen_ec50"].astype(float).to_numpy()
    x = x[np.isfinite(x) & (x > 0)]
    if len(x) == 0:
        return

    bins = make_log_histogram_bins(x, dilution_factors)
    colors = ["#1f77b4", "#d62728", "#2ca02c", "#ff7f0e"]

    fig, ax = plt.subplots(figsize=(8, 5.625), dpi=220)
    ax.hist(x, bins=bins, color="#c8c8c8", edgecolor="#5a5a5a", linewidth=1.0)
    ymax = ax.get_ylim()[1]
    for idx, (_, row) in enumerate(controls.iterrows(), start=1):
        val = float(row["antigen_ec50"])
        color = colors[(idx - 1) % len(colors)]
        ax.axvline(val, linewidth=2, color=color)
        ax.text(
            val,
            ymax * max(0.55, 0.95 - 0.08 * (idx - 1)),
            f"{row['sample_id']}\n{val:.2f}",
            rotation=90,
            va="top",
            ha="center",
            color=color,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=1.5),
        )

    patient_median = float(np.median(x))
    ax.axvline(patient_median, color="black", linestyle="--", linewidth=1.5)
    ax.text(
        patient_median,
        ymax * 0.82,
        f"patient median\n{patient_median:.2f}",
        rotation=90,
        va="top",
        ha="center",
        color="black",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=1.5),
    )
    ax.set_xscale("log")
    ax.set_xlabel("EC50 (dilution factor)")
    ax.set_ylabel("Patients (count)")
    ax.set_title(f"Figure 1. {antigen_name} EC50 distribution")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(outdir / f"figure_1_histogram_{antigen_slug}_ec50.svg", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5.625), dpi=220)
    ax.hist(x, bins=bins, color="#c8c8c8", edgecolor="#5a5a5a", linewidth=1.0)
    if len(x) >= 2:
        logx = np.log10(x)
        mu = float(np.mean(logx))
        sigma = float(np.std(logx, ddof=1))
        log_bins = np.log10(bins)
        binw = np.diff(log_bins).mean()
        ygrid = np.linspace(log_bins.min(), log_bins.max(), 500)
        pdf = (1.0 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((ygrid - mu) / sigma) ** 2)
        ax.plot(10 ** ygrid, len(logx) * binw * pdf, linewidth=2.0, label="Gaussian fit (log10 EC50)")

    ymax = ax.get_ylim()[1]
    for idx, (_, row) in enumerate(controls.iterrows(), start=1):
        val = float(row["antigen_ec50"])
        color = colors[(idx - 1) % len(colors)]
        ax.axvline(val, linewidth=2, color=color)
        ax.text(
            val,
            ymax * max(0.55, 0.95 - 0.08 * (idx - 1)),
            f"{row['sample_id']}\n{val:.2f}",
            rotation=90,
            va="top",
            ha="center",
            color=color,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=1.5),
        )

    ax.axvline(patient_median, color="black", linestyle="--", linewidth=1.5)
    ax.set_xscale("log")
    ax.set_xlabel("EC50 (dilution factor)")
    ax.set_ylabel("Patients (count)")
    ax.set_title(f"{antigen_name} EC50 distribution with Gaussian fit")
    ax.grid(axis="y", alpha=0.3)
    if len(x) >= 2:
        ax.legend(frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(outdir / f"figure_1_histogram_{antigen_slug}_ec50_with_gaussian.svg", bbox_inches="tight")
    plt.close(fig)


def plot_figure_2_scatter(df: pd.DataFrame, outdir: Path, antigen_name: str, dilution_factors: np.ndarray) -> None:
    patients_grey = df[(df["group"] == "Patient") & (df["grey_flag"])].copy()
    patients_black = df[(df["group"] == "Patient") & (~df["grey_flag"])].copy()
    controls = df[df["group"] == "Control"].copy()

    fig, ax = plt.subplots(figsize=(6.5, 6.5), dpi=220)
    ax.scatter(
        patients_grey["blocked_ec50"],
        patients_grey["antigen_ec50"],
        s=28,
        color="lightgrey",
        label=f"Patients: blocked > {antigen_name}/2",
    )
    ax.scatter(
        patients_black["blocked_ec50"],
        patients_black["antigen_ec50"],
        s=28,
        color="black",
        label="Other patients",
    )
    if not controls.empty:
        ax.scatter(controls["blocked_ec50"], controls["antigen_ec50"], s=42, color="red", label="Controls")

    for _, row in controls.iterrows():
        ax.annotate(
            row["sample_id"],
            (row["blocked_ec50"], row["antigen_ec50"]),
            xytext=(6, 6),
            textcoords="offset points",
            color="red",
        )

    configure_log_ec50_axes(ax, df, dilution_factors)
    ax.set_xlabel("EC50 blocked")
    ax.set_ylabel(f"EC50 {antigen_name}")
    ax.set_title(f"Figure 2. {antigen_name} EC50 versus blocked EC50")
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(outdir / "figure_2_scatter_grey_subset.svg", bbox_inches="tight")
    plt.close(fig)


def plot_figure_2_with_uncertainty(
    df: pd.DataFrame,
    outdir: Path,
    antigen_name: str,
    dilution_factors: np.ndarray,
) -> None:
    patients_grey = df[(df["group"] == "Patient") & (df["grey_flag"])].copy()
    patients_black = df[(df["group"] == "Patient") & (~df["grey_flag"])].copy()
    controls = df[df["group"] == "Control"].copy()
    fig, ax = plt.subplots(figsize=(6.5, 6.5), dpi=220)

    def errplot(subset: pd.DataFrame, color: str, ms: float, zorder: int) -> None:
        ax.errorbar(
            subset["blocked_ec50"],
            subset["antigen_ec50"],
            xerr=np.vstack([subset["xerr_low"], subset["xerr_high"]]),
            yerr=np.vstack([subset["yerr_low"], subset["yerr_high"]]),
            fmt="o",
            markersize=ms,
            color=color,
            ecolor=color,
            elinewidth=0.8 if color != "red" else 1.0,
            capsize=0,
            linestyle="none",
            zorder=zorder,
        )

    errplot(patients_grey, "lightgrey", 4.5, 1)
    errplot(patients_black, "black", 4.5, 2)
    if not controls.empty:
        errplot(controls, "red", 5.5, 3)

    for _, row in controls.iterrows():
        ax.annotate(
            row["sample_id"],
            (row["blocked_ec50"], row["antigen_ec50"]),
            xytext=(6, 6),
            textcoords="offset points",
            color="red",
            fontsize=9,
        )

    configure_log_ec50_axes(ax, df, dilution_factors)
    ax.set_xlabel("EC50 blocked")
    ax.set_ylabel(f"EC50 {antigen_name}")
    ax.set_title(f"Figure 2. {antigen_name} EC50 versus blocked EC50 with uncertainty")
    fig.tight_layout()
    fig.savefig(outdir / "figure_2_scatter_with_uncertainty.svg", bbox_inches="tight")
    plt.close(fig)

    df[
        [
            "sample_id",
            "group",
            "antigen_ec50",
            "blocked_ec50",
            "antigen_se_log10_ec50",
            "blocked_se_log10_ec50",
            "antigen_rmse",
            "blocked_rmse",
            "xerr_low",
            "xerr_high",
            "yerr_low",
            "yerr_high",
            "grey_flag",
        ]
    ].to_csv(outdir / "figure_2_scatter_with_uncertainty_source_data.csv", index=False)

def plot_exploratory_figures(
    df: pd.DataFrame,
    outdir: Path,
    antigen_name: str,
    dilution_factors: np.ndarray,
) -> None:
    patients = df[df["group"] == "Patient"].copy()
    controls = df[df["group"] == "Control"].copy()
    if patients.empty:
        return

    n_dilutions = len(dilution_factors)
    antigen_mat = np.vstack([dilution_values(row, "antigen", n_dilutions) for _, row in patients.iterrows()])
    blocked_mat = np.vstack([dilution_values(row, "blocked", n_dilutions) for _, row in patients.iterrows()])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        antigen_med = np.nanmedian(antigen_mat, axis=0)
        blocked_med = np.nanmedian(blocked_mat, axis=0)

    fig, ax = plt.subplots(figsize=(7, 5), dpi=220)
    ax.plot(dilution_factors, antigen_med, marker="o", label=f"Patient median {antigen_name}")
    ax.plot(dilution_factors, blocked_med, marker="o", label="Patient median blocked")
    ax.set_xscale("log")
    ax.set_xlabel("Dilution factor")
    ax.set_ylabel("Optical density")
    ax.set_title("Exploratory: cohort median curves")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outdir / "exploratory_cohort_median_curves.svg", bbox_inches="tight")
    plt.close(fig)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        antigen_norm = antigen_mat / np.nanmax(antigen_mat, axis=1, keepdims=True)
        blocked_norm = blocked_mat / np.nanmax(blocked_mat, axis=1, keepdims=True)
        antigen_norm_med = np.nanmedian(antigen_norm, axis=0)
        blocked_norm_med = np.nanmedian(blocked_norm, axis=0)

    fig, ax = plt.subplots(figsize=(7, 5), dpi=220)
    ax.plot(dilution_factors, antigen_norm_med, marker="o", label=f"Patient median {antigen_name} (normalized)")
    ax.plot(dilution_factors, blocked_norm_med, marker="o", label="Patient median blocked (normalized)")
    ax.set_xscale("log")
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("Dilution factor")
    ax.set_ylabel("Normalized optical density")
    ax.set_title("Exploratory: normalized shape comparison")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outdir / "exploratory_normalized_shape_comparison.svg", bbox_inches="tight")
    plt.close(fig)

    with np.errstate(divide="ignore", invalid="ignore"):
        ratio_mat = antigen_mat / blocked_mat
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        ratio_med = np.nanmedian(ratio_mat, axis=0)

    fig, ax = plt.subplots(figsize=(7, 5), dpi=220)
    ax.plot(dilution_factors, ratio_med, marker="o")
    ax.axhline(1.0, color="gray", linestyle=":")
    ax.set_xscale("log")
    ax.set_xlabel("Dilution factor")
    ax.set_ylabel(f"{antigen_name} / blocked ratio")
    ax.set_title("Exploratory: patient median antigen/blocked ratio by dilution")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / "exploratory_ratio_vs_dilution.svg", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.5, 5.5), dpi=220)
    ax.scatter(patients["log10_ratio_antigen_blocked_ec50"], patients["net_auc"], s=24, color="black")
    if not controls.empty:
        ax.scatter(controls["log10_ratio_antigen_blocked_ec50"], controls["net_auc"], s=36, color="red")
    for _, row in controls.iterrows():
        ax.annotate(
            row["sample_id"],
            (row["log10_ratio_antigen_blocked_ec50"], row["net_auc"]),
            xytext=(6, 6),
            textcoords="offset points",
            color="red",
        )
    ax.axvline(0, color="gray", linestyle=":")
    ax.set_xlabel(f"log10(EC50 {antigen_name} / EC50 blocked)")
    ax.set_ylabel(f"AUC({antigen_name} - blocked)")
    ax.set_title("Exploratory: specificity summary scatter")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / "exploratory_specificity_scatter.svg", bbox_inches="tight")
    plt.close(fig)

    pd.DataFrame(
        {
            "dilution_factor": dilution_factors,
            "patient_median_antigen": antigen_med,
            "patient_median_blocked": blocked_med,
            "patient_median_antigen_norm": antigen_norm_med,
            "patient_median_blocked_norm": blocked_norm_med,
            "patient_median_antigen_over_blocked": ratio_med,
        }
    ).to_csv(outdir / "exploratory_cohort_medians.csv", index=False)


def write_legends(outdir: Path, antigen_name: str, n_dilutions: int) -> None:
    legend_text = f"""Figure 1. Distribution of {antigen_name} EC50 values in the patient cohort.

This figure shows the distribution of EC50 values obtained in the {antigen_name} ELISA for the patient cohort made available for analysis. The histogram includes only patient samples. If control samples are present in both paired workbooks, they are shown separately and labeled with their individual EC50 values so that their positions can be compared directly with the patient distribution.

In this analysis, EC50 refers to the dilution factor at which the fitted ELISA signal reaches one-half of its dynamic range for a given sample. The paired workbooks provide a fitted EC50 term on a log10 scale, and the script converts that quantity back to dilution-factor units for plotting and downstream comparisons. The optical-density values in the selected measurement columns are used for summary curves, ratio plots, and approximate uncertainty estimates.

Figure 2. Comparison of EC50 values obtained on {antigen_name}-coated plates and blocked control plates, including per-sample uncertainty estimates.

This figure compares, for each individual sample, the EC50 obtained in the {antigen_name}-coated ELISA with the EC50 obtained on the corresponding blocked control plate. The y-axis shows the EC50 from the antigen-coated assay, and the x-axis shows the EC50 from the matched blocked plate that contains no antigen. The blocked plate is used as an estimate of background or nonspecific binding.

Each point represents one sample analyzed under both conditions. Patient samples are plotted as black symbols and, when present, control samples are plotted as red symbols. A subset of patient samples is highlighted in light gray. These are the samples for which EC50(blocked) > 0.5 x EC50({antigen_name}). These light-gray samples are flagged because they may not represent convincing antigen-specific reactivity.

Per-sample uncertainty was approximated by refitting a three-parameter logistic model to the selected {n_dilutions}-point dilution series while keeping the reported EC50 fixed. Residual deviation between the observed points and the constrained fitted curve was summarized as a root mean square error (RMSE). This residual scatter in signal space was converted into an approximate uncertainty in log10(EC50) using the local slope of the logistic curve at the inflection point. When too few finite measurement points were available for this constrained refit, the uncertainty terms were left undefined.
"""
    (outdir / "detailed_legends.txt").write_text(legend_text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare paired ELISA curve exports for an antigen-coated plate and a blocked control plate."
    )
    parser.add_argument("--antigen-input", help="Path to the antigen workbook. If omitted, a file picker opens.")
    parser.add_argument("--blocked-input", help="Path to the blocked workbook. If omitted, a file picker opens.")
    parser.add_argument("--antigen-name", help="Label to use in figure titles and tables.")
    parser.add_argument("--measurement-columns", default="K:N", help="Measurement columns to read from each workbook.")
    parser.add_argument("--outdir", help="Output directory. Defaults to a sibling folder next to the antigen workbook.")
    args = parser.parse_args()

    antigen_path = Path(args.antigen_input) if args.antigen_input else choose_input_workbook("Select antigen workbook")
    blocked_path = Path(args.blocked_input) if args.blocked_input else choose_input_workbook("Select blocked-control workbook")
    antigen_name = args.antigen_name or choose_antigen_name(antigen_path)
    antigen_slug = sanitize_slug(antigen_name)
    measurement_cols = parse_column_range(args.measurement_columns)
    outdir = Path(args.outdir) if args.outdir else antigen_path.with_name(f"{antigen_slug}_vs_blocked_output")
    outdir.mkdir(parents=True, exist_ok=True)

    antigen_df, antigen_dilutions = load_condition_workbook(antigen_path, measurement_cols, "antigen")
    blocked_df, blocked_dilutions = load_condition_workbook(blocked_path, measurement_cols, "blocked")
    if len(antigen_dilutions) != len(blocked_dilutions) or not np.allclose(
        antigen_dilutions, blocked_dilutions, rtol=1e-3, atol=1e-6
    ):
        raise SystemExit("The antigen and blocked workbooks do not share the same dilution headers.")

    merged = compute_derived_metrics(merge_condition_data(antigen_df, blocked_df), antigen_dilutions)
    save_basic_data_tables(merged, outdir, antigen_name)
    plot_figure_1_histogram(merged, outdir, antigen_name, antigen_slug, antigen_dilutions)
    plot_figure_2_scatter(merged, outdir, antigen_name, antigen_dilutions)
    plot_figure_2_with_uncertainty(merged, outdir, antigen_name, antigen_dilutions)
    plot_exploratory_figures(merged, outdir, antigen_name, antigen_dilutions)
    write_legends(outdir, antigen_name, len(antigen_dilutions))

    print(f"Antigen workbook: {antigen_path.resolve()}")
    print(f"Blocked workbook: {blocked_path.resolve()}")
    print(f"Samples analyzed: {len(merged)}")
    print(f"Done. Outputs written to: {outdir.resolve()}")


if __name__ == "__main__":
    main()
