#!/usr/bin/env python3
"""
Extract the analysis panel from the full PRAMStat master dataset.

This script produces data/analysis_panel.csv — the wide-format state-year
panel used by run_pcha_analysis.py. If you have the full master dataset,
run this script first. If you received only the repository, the pre-built
panel is included in data/.

Usage:
    python scripts/extract_panel_from_master.py [--master PATH]
"""
import argparse
import hashlib
from pathlib import Path
import pandas as pd
import numpy as np


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--master",
        default=str(Path(__file__).resolve().parent.parent / "data" / "analysis_master.csv"),
        help="Path to the full analysis_master.csv (6.4M rows)",
    )
    args = parser.parse_args()
    master_path = Path(args.master)
    out_dir = Path(__file__).resolve().parent.parent / "data"

    print(f"Loading master dataset: {master_path}")
    print(f"SHA-256: {sha256_file(master_path)[:16]}...")

    df = pd.read_csv(master_path, low_memory=False)
    print(f"Loaded {len(df):,} rows")

    # Filter: unstratified, YES response
    mask = (
        df["break_out"].isna() &
        df["break_out_category"].isna() &
        (df["binary_response_class"] == "yes")
    )
    sub = df[mask].copy()
    sub["estimate"] = pd.to_numeric(sub["estimate"], errors="coerce")
    sub["std_err"] = pd.to_numeric(sub["std_err"], errors="coerce")

    # Wide pivot: estimate
    panel = sub.pivot_table(
        index=["location_abbr", "year"],
        columns="question_id",
        values="estimate",
        aggfunc="first",
    ).reset_index()
    panel.columns.name = None

    # SE pivot
    se_panel = sub.pivot_table(
        index=["location_abbr", "year"],
        columns="question_id",
        values="std_err",
        aggfunc="first",
    ).reset_index()
    se_panel.columns.name = None

    # Harmonized PPD outcome
    if "QUO74" in panel.columns and "QUO219" in panel.columns:
        panel["outcome_ppd"] = panel["QUO219"].combine_first(panel["QUO74"])
    elif "QUO219" in panel.columns:
        panel["outcome_ppd"] = panel["QUO219"]
    elif "QUO74" in panel.columns:
        panel["outcome_ppd"] = panel["QUO74"]

    panel["ppd_era"] = np.where(
        panel.get("QUO219", pd.Series(dtype=float)).notna(), "2009-2011",
        np.where(panel.get("QUO74", pd.Series(dtype=float)).notna(), "2004-2008", "none")
    )

    # Merge SE for PPD
    if "QUO74" in se_panel.columns and "QUO219" in se_panel.columns:
        se_panel["se_ppd"] = se_panel["QUO219"].combine_first(se_panel["QUO74"])
    elif "QUO219" in se_panel.columns:
        se_panel["se_ppd"] = se_panel["QUO219"]
    elif "QUO74" in se_panel.columns:
        se_panel["se_ppd"] = se_panel["QUO74"]

    if "se_ppd" in se_panel.columns:
        panel = panel.merge(
            se_panel[["location_abbr", "year", "se_ppd"]],
            on=["location_abbr", "year"], how="left",
        )

    # Question dictionary
    qdict = (
        sub.groupby("question_id")
        .agg(
            question_text=("question", "first"),
            min_year=("year", "min"),
            max_year=("year", "max"),
            n_rows=("estimate", "size"),
        )
        .reset_index()
    )

    panel.to_csv(out_dir / "analysis_panel.csv", index=False)
    se_panel.to_csv(out_dir / "se_panel.csv", index=False)
    qdict.to_csv(out_dir / "question_dictionary.csv", index=False)

    print(f"\nPanel saved: {len(panel)} rows × {panel.shape[1]} cols -> {out_dir / 'analysis_panel.csv'}")
    print(f"SE panel saved: {se_panel.shape} -> {out_dir / 'se_panel.csv'}")
    print(f"Question dictionary: {len(qdict)} entries -> {out_dir / 'question_dictionary.csv'}")


if __name__ == "__main__":
    main()
