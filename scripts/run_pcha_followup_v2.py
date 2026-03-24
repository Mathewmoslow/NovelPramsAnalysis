#!/usr/bin/env python3
"""
PCHA Follow-Up Analysis v2
===========================
Addresses all outstanding items from the v1 sanity check:

CRITICAL:
  1. Compute Cronbach's alpha + parallel analysis on complete-case subsets
  2. Fix predicted-benefit sign/label
  3. Strict pre-conception-only PCHA sensitivity (excluding postpartum components)

IMPORTANT:
  4. Report VIF for interaction models
  5. Investigate QUO16 absence
  6. Strict pre-conception PCHA (the highest-value-add analysis)
  7. Within-state partial permutation test

All outputs go to output/<run_id>/ alongside v1 outputs (never overwrites).

Usage:
    python scripts/run_pcha_followup_v2.py [--panel PATH] [--run-id ID]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Constants (must match v1 for consistency)
# ---------------------------------------------------------------------------
RANDOM_SEED = 42
N_PERMUTATIONS = 5000
N_BOOTSTRAP = 5000
MIN_N = 30

RISK_ANCHORS = {
    "QUO197": "Argued with husband/partner more than usual (12 mo before pregnancy)",
    "QUO210": "Any partner-related stressors reported",
    "QUO313": "Physical IPV by ex-partner (before pregnancy)",
    "QUO315": "Physical IPV by ex-partner (during pregnancy)",
}

# FULL PCHA (same as v1)
PCHA_ALL = {
    "QUO179": "Exercised 3+ days/week in 3 months before pregnancy",
    "QUO41":  "Took multivitamins >4x/week month before pregnancy",
    "QUO65":  "Took daily multivitamin in month before pregnancy",
    "QUO249": "Had teeth cleaned by dentist before pregnancy",
    "QUO75":  "Had teeth cleaned during pregnancy",
    "QUO257": "Was trying to become pregnant",
    "QUO16":  "Pregnancy was intended",
    "QUO296": "Got prenatal care as early as wanted (2000-2008)",
    "QUO297": "Prenatal care began in first trimester (2009-2011)",
    "QUO4":   "Still breastfeeding at 4 weeks postpartum",
    "QUO5":   "Ever breastfed or pumped breast milk",
    "QUO44":  "Still breastfeeding at 8 weeks postpartum",
    "QUO101": "Baby had checkup/exam within first week",
}

# STRICT PRE-CONCEPTION PCHA: only items measured BEFORE or DURING pregnancy
# Excludes all postpartum measurements (breastfeeding, infant checkup)
PCHA_STRICT_PRECONCEPTION = {
    "QUO179": "Exercised 3+ days/week in 3 months before pregnancy",
    "QUO41":  "Took multivitamins >4x/week month before pregnancy",
    "QUO65":  "Took daily multivitamin in month before pregnancy",
    "QUO249": "Had teeth cleaned by dentist before pregnancy",
    "QUO75":  "Had teeth cleaned during pregnancy",
    "QUO257": "Was trying to become pregnant",
    "QUO296": "Got prenatal care as early as wanted (2000-2008)",
    "QUO297": "Prenatal care began in first trimester (2009-2011)",
}

# NARROW PRE-CONCEPTION: only items measured BEFORE pregnancy
# (the strictest temporal test — nothing prenatal or postpartum)
PCHA_NARROW_PREPREG = {
    "QUO179": "Exercised 3+ days/week in 3 months before pregnancy",
    "QUO41":  "Took multivitamins >4x/week month before pregnancy",
    "QUO65":  "Took daily multivitamin in month before pregnancy",
    "QUO249": "Had teeth cleaned by dentist before pregnancy",
    "QUO257": "Was trying to become pregnant",
}

# 7-component subset available for all 188 rows (for psychometrics)
PCHA_FULL_COVERAGE = {
    "QUO41":  "Took multivitamins >4x/week month before pregnancy",
    "QUO65":  "Took daily multivitamin in month before pregnancy",
    "QUO257": "Was trying to become pregnant",
    "QUO4":   "Still breastfeeding at 4 weeks postpartum",
    "QUO5":   "Ever breastfed or pumped breast milk",
    "QUO44":  "Still breastfeeding at 8 weeks postpartum",
    "QUO101": "Baby had checkup/exam within first week",
}


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------
def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def zscore(x: np.ndarray) -> np.ndarray:
    s = float(np.nanstd(x, ddof=1))
    if not math.isfinite(s) or s < 1e-12:
        return np.full_like(x, np.nan, dtype=float)
    return (x - float(np.nanmean(x))) / s


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = np.empty(n, dtype=float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = ranked[i] * n / rank
        prev = min(prev, val)
        q[i] = prev
    out = np.empty(n, dtype=float)
    out[order] = np.clip(q, 0, 1)
    return out


def fit_ols(y, X_cols, col_names=None, weights=None):
    X = np.column_stack([np.ones(len(y)), *X_cols])
    n, k = X.shape
    if n <= k or n < 8:
        return None
    if col_names is None:
        col_names = ["intercept"] + [f"x{i}" for i in range(len(X_cols))]
    else:
        col_names = ["intercept"] + list(col_names)
    if weights is not None:
        W = np.diag(np.sqrt(weights))
        Xw = W @ X
        yw = W @ y
    else:
        Xw = X
        yw = y
    try:
        beta, residuals, rank, sv = np.linalg.lstsq(Xw, yw, rcond=None)
    except np.linalg.LinAlgError:
        return None
    if rank < k:
        return None
    resid = y - X @ beta
    rss = float(np.sum(resid ** 2))
    tss = float(np.sum((y - np.mean(y)) ** 2))
    if tss <= 0:
        return None
    r2 = 1.0 - rss / tss
    adj_r2 = 1.0 - (1.0 - r2) * ((n - 1) / max(n - k, 1))
    rmse = math.sqrt(rss / n)
    aic = n * math.log(max(rss / n, 1e-12)) + 2 * k
    bic = n * math.log(max(rss / n, 1e-12)) + k * math.log(n)
    s2 = rss / max(n - k, 1)
    try:
        XtX_inv = np.linalg.inv(Xw.T @ Xw)
    except np.linalg.LinAlgError:
        return None
    cov = s2 * XtX_inv
    se = np.sqrt(np.maximum(np.diag(cov), 0))
    tvals = np.where(se > 0, beta / se, 0)
    pvals = 2 * (1 - stats.t.cdf(np.abs(tvals), df=max(n - k, 1)))
    return {
        "beta": beta, "se": se, "pvals": pvals, "n": n, "k": k,
        "r2": r2, "adj_r2": adj_r2, "rmse": rmse, "aic": aic, "bic": bic,
        "rss": rss, "tss": tss, "col_names": col_names, "XtX_inv": XtX_inv,
    }


def fit_ols_cluster_robust(y, X_cols, cluster_ids, col_names=None):
    X = np.column_stack([np.ones(len(y)), *X_cols])
    n, k = X.shape
    if n <= k or n < 8:
        return None
    if col_names is None:
        col_names = ["intercept"] + [f"x{i}" for i in range(len(X_cols))]
    else:
        col_names = ["intercept"] + list(col_names)
    try:
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    except np.linalg.LinAlgError:
        return None
    resid = y - X @ beta
    rss = float(np.sum(resid ** 2))
    tss = float(np.sum((y - np.mean(y)) ** 2))
    if tss <= 0:
        return None
    r2 = 1.0 - rss / tss
    adj_r2 = 1.0 - (1.0 - r2) * ((n - 1) / max(n - k, 1))
    try:
        XtX_inv = np.linalg.inv(X.T @ X)
    except np.linalg.LinAlgError:
        return None
    unique_clusters = np.unique(cluster_ids)
    G = len(unique_clusters)
    meat = np.zeros((k, k))
    for c in unique_clusters:
        mask = cluster_ids == c
        Xc = X[mask]
        ec = resid[mask]
        score_c = Xc.T @ ec
        meat += np.outer(score_c, score_c)
    correction = (G / (G - 1)) * ((n - 1) / (n - k))
    V_cr = correction * XtX_inv @ meat @ XtX_inv
    se = np.sqrt(np.maximum(np.diag(V_cr), 0))
    tvals = np.where(se > 0, beta / se, 0)
    df_cr = max(G - 1, 1)
    pvals = 2 * (1 - stats.t.cdf(np.abs(tvals), df=df_cr))
    return {
        "beta": beta, "se": se, "pvals": pvals, "n": n, "k": k,
        "r2": r2, "adj_r2": adj_r2, "col_names": col_names,
    }


def make_fe_dummies(panel):
    states = panel["location_abbr"].values
    years = panel["year"].values.astype(int)
    u_states = sorted(set(states))
    u_years = sorted(set(years))
    dummies = []
    names = []
    for s in u_states[1:]:
        dummies.append((states == s).astype(float))
        names.append(f"FE_state_{s}")
    for yr in u_years[1:]:
        dummies.append((years == yr).astype(float))
        names.append(f"FE_year_{yr}")
    return dummies, names


def parallel_analysis(data_matrix, n_iter=1000, seed=RANDOM_SEED):
    rng = np.random.RandomState(seed)
    n, p = data_matrix.shape
    all_eigs = np.zeros((n_iter, p))
    for i in range(n_iter):
        rand_data = rng.standard_normal((n, p))
        corr = np.corrcoef(rand_data, rowvar=False)
        eigs = np.linalg.eigvalsh(corr)[::-1]
        all_eigs[i] = eigs
    return np.percentile(all_eigs, 95, axis=0)


def compute_cronbach_alpha(df_items):
    """Cronbach's alpha on a DataFrame of items (rows=observations, cols=items)."""
    df_clean = df_items.dropna()
    if len(df_clean) < 10 or df_clean.shape[1] < 2:
        return float("nan"), len(df_clean), df_clean.shape[1]
    k = df_clean.shape[1]
    item_vars = df_clean.var(ddof=1).values
    total_var = df_clean.sum(axis=1).var(ddof=1)
    if total_var <= 0:
        return float("nan"), len(df_clean), k
    alpha = (k / (k - 1)) * (1 - item_vars.sum() / total_var)
    return float(alpha), len(df_clean), k


def compute_vif(X_cols, col_names):
    """Compute VIF for each predictor column."""
    X = np.column_stack(X_cols)
    n, p = X.shape
    vifs = {}
    for j in range(p):
        y_j = X[:, j]
        others = [X[:, i] for i in range(p) if i != j]
        if not others:
            vifs[col_names[j]] = 1.0
            continue
        Xo = np.column_stack([np.ones(n)] + others)
        try:
            beta, *_ = np.linalg.lstsq(Xo, y_j, rcond=None)
        except np.linalg.LinAlgError:
            vifs[col_names[j]] = float("inf")
            continue
        resid = y_j - Xo @ beta
        rss = float(np.sum(resid ** 2))
        tss = float(np.sum((y_j - np.mean(y_j)) ** 2))
        if tss <= 0:
            vifs[col_names[j]] = float("inf")
            continue
        r2_j = 1.0 - rss / tss
        vifs[col_names[j]] = 1.0 / max(1.0 - r2_j, 1e-12)
    return vifs


def build_pcha_index(panel, component_dict, min_n=MIN_N, min_components=3):
    """Build a z-score composite index from a dict of QUO->description."""
    avail = [q for q in component_dict if q in panel.columns
             and panel[q].notna().sum() >= min_n]
    z_cols = {}
    for qid in avail:
        z_cols[qid] = zscore(panel[qid].values.astype(float))
    if not z_cols:
        return np.full(len(panel), np.nan), avail, z_cols
    z_matrix = np.column_stack([z_cols[q] for q in avail])
    n_valid = np.sum(~np.isnan(z_matrix), axis=1)
    index = np.nanmean(z_matrix, axis=1)
    index[n_valid < min_components] = np.nan
    return index, avail, z_cols


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="PCHA Follow-Up Analysis v2")
    parser.add_argument(
        "--panel",
        default=str(Path(__file__).resolve().parent.parent / "data" / "analysis_panel.csv"),
    )
    parser.add_argument("--run-id", default="")
    parser.add_argument("--seed", type=int, default=RANDOM_SEED)
    args = parser.parse_args()

    run_id = args.run_id or datetime.now(timezone.utc).strftime("run_%Y%m%dT%H%M%SZ_pcha_v2")
    seed = args.seed
    rng = np.random.RandomState(seed)

    base_dir = Path(__file__).resolve().parent.parent
    out_dir = base_dir / "output" / run_id
    fig_dir = base_dir / "figures" / run_id
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    panel_path = Path(args.panel)
    panel_hash = sha256_file(panel_path)
    panel = pd.read_csv(panel_path, low_memory=False)

    ppd_panel = panel[panel["outcome_ppd"].notna()].copy()
    ppd_panel["outcome_ppd"] = ppd_panel["outcome_ppd"].astype(float)

    print(f"[PCHA v2] Run ID: {run_id}")
    print(f"[PCHA v2] Seed: {seed}")
    print(f"[PCHA v2] Panel: {len(panel)} rows, PPD panel: {len(ppd_panel)} rows")
    print(f"[PCHA v2] SHA-256: {panel_hash[:16]}...")

    metadata = {
        "run_id": run_id,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": sys.version,
        "platform": platform.platform(),
        "panel_sha256": panel_hash,
        "panel_rows": len(panel),
        "ppd_panel_rows": len(ppd_panel),
        "random_seed": seed,
        "purpose": "Follow-up addressing v1 sanity check items",
        "parent_run": "run_20260324T_pcha_v1",
    }

    # ===================================================================
    # ITEM 1: CRONBACH'S ALPHA + PARALLEL ANALYSIS
    # ===================================================================
    print("\n" + "=" * 60)
    print("ITEM 1: PSYCHOMETRIC DIAGNOSTICS")
    print("=" * 60)

    psychometric_results = []

    # --- 1A: 7-component full-coverage subset (n=188) ---
    print("\n--- 1A: 7-component full-coverage subset ---")
    fc_avail = [q for q in PCHA_FULL_COVERAGE if q in ppd_panel.columns]
    fc_z_df = pd.DataFrame({q: zscore(ppd_panel[q].values.astype(float)) for q in fc_avail},
                           index=ppd_panel.index)
    alpha_7, n_alpha_7, k_alpha_7 = compute_cronbach_alpha(fc_z_df)
    print(f"  Cronbach's alpha (7 components, n={n_alpha_7}): {alpha_7:.4f}")

    # Parallel analysis on 7-component subset
    fc_clean = fc_z_df.dropna()
    if len(fc_clean) >= 20:
        data_mat = fc_clean.values
        corr_mat = np.corrcoef(data_mat, rowvar=False)
        actual_eigs = np.sort(np.linalg.eigvalsh(corr_mat))[::-1]
        pa_thresholds = parallel_analysis(data_mat, n_iter=1000, seed=seed)

        pa_7_rows = []
        for i in range(len(actual_eigs)):
            pa_7_rows.append({
                "component": i + 1,
                "actual_eigenvalue": float(actual_eigs[i]),
                "pa_95th_threshold": float(pa_thresholds[i]),
                "retain": actual_eigs[i] > pa_thresholds[i],
            })
        pa_7_df = pd.DataFrame(pa_7_rows)
        pa_7_df.to_csv(out_dir / "parallel_analysis_7comp.csv", index=False)
        n_retain_7 = int(pa_7_df["retain"].sum())
        print(f"  Parallel analysis (7 comp): retain {n_retain_7} factors")
        for _, r in pa_7_df.iterrows():
            marker = " ***" if r["retain"] else ""
            print(f"    PC{int(r['component'])}: eig={r['actual_eigenvalue']:.3f} "
                  f"vs thresh={r['pa_95th_threshold']:.3f}{marker}")
    else:
        pa_7_df = pd.DataFrame()
        n_retain_7 = 0

    psychometric_results.append({
        "subset": "7-component full-coverage (n=188)",
        "components": ", ".join(fc_avail),
        "n_obs": n_alpha_7,
        "k_items": k_alpha_7,
        "cronbach_alpha": alpha_7,
        "pa_factors_retained": n_retain_7,
    })

    # --- 1B: All 12 components on complete-case subset ---
    print("\n--- 1B: All 12 components, complete-case subset ---")
    all_avail = [q for q in PCHA_ALL if q in ppd_panel.columns]
    all_z_df = pd.DataFrame({q: zscore(ppd_panel[q].values.astype(float)) for q in all_avail},
                            index=ppd_panel.index)
    alpha_12, n_alpha_12, k_alpha_12 = compute_cronbach_alpha(all_z_df)
    print(f"  Cronbach's alpha ({k_alpha_12} components, n={n_alpha_12}): {alpha_12:.4f}")

    # Parallel analysis on 12-component complete case
    all_clean = all_z_df.dropna()
    if len(all_clean) >= 20:
        data_mat_12 = all_clean.values
        corr_mat_12 = np.corrcoef(data_mat_12, rowvar=False)
        actual_eigs_12 = np.sort(np.linalg.eigvalsh(corr_mat_12))[::-1]
        pa_thresh_12 = parallel_analysis(data_mat_12, n_iter=1000, seed=seed)

        pa_12_rows = []
        for i in range(len(actual_eigs_12)):
            pa_12_rows.append({
                "component": i + 1,
                "actual_eigenvalue": float(actual_eigs_12[i]),
                "pa_95th_threshold": float(pa_thresh_12[i]),
                "retain": actual_eigs_12[i] > pa_thresh_12[i],
            })
        pa_12_df = pd.DataFrame(pa_12_rows)
        pa_12_df.to_csv(out_dir / "parallel_analysis_12comp.csv", index=False)
        n_retain_12 = int(pa_12_df["retain"].sum())
        print(f"  Parallel analysis ({k_alpha_12} comp): retain {n_retain_12} factors")
        for _, r in pa_12_df.iterrows():
            marker = " ***" if r["retain"] else ""
            print(f"    PC{int(r['component'])}: eig={r['actual_eigenvalue']:.3f} "
                  f"vs thresh={r['pa_95th_threshold']:.3f}{marker}")

        # PCA loadings
        evals_12, evecs_12 = np.linalg.eigh(corr_mat_12)
        idx_sort = np.argsort(evals_12)[::-1]
        evals_12 = evals_12[idx_sort]
        evecs_12 = evecs_12[:, idx_sort]
        if n_retain_12 >= 1:
            loadings_12 = evecs_12[:, :n_retain_12] * np.sqrt(evals_12[:n_retain_12])
            loading_df = pd.DataFrame(
                loadings_12,
                index=[q for q in all_avail if q in all_clean.columns],
                columns=[f"PC{i+1}" for i in range(n_retain_12)],
            )
            loading_df.to_csv(out_dir / "pca_loadings_12comp.csv")
            print(f"  PCA loadings saved ({n_retain_12} factors)")
    else:
        pa_12_df = pd.DataFrame()
        n_retain_12 = 0
        print(f"  Insufficient complete cases (n={len(all_clean)}) for parallel analysis")

    psychometric_results.append({
        "subset": f"{k_alpha_12}-component complete-case",
        "components": ", ".join(all_avail),
        "n_obs": n_alpha_12,
        "k_items": k_alpha_12,
        "cronbach_alpha": alpha_12,
        "pa_factors_retained": n_retain_12,
    })

    # --- 1C: Strict pre-conception components ---
    print("\n--- 1C: Strict pre-conception components ---")
    strict_avail = [q for q in PCHA_STRICT_PRECONCEPTION if q in ppd_panel.columns]
    strict_z_df = pd.DataFrame({q: zscore(ppd_panel[q].values.astype(float)) for q in strict_avail},
                               index=ppd_panel.index)
    alpha_strict, n_alpha_strict, k_alpha_strict = compute_cronbach_alpha(strict_z_df)
    print(f"  Cronbach's alpha ({k_alpha_strict} strict-precon, n={n_alpha_strict}): {alpha_strict:.4f}")

    psychometric_results.append({
        "subset": f"{k_alpha_strict}-component strict pre-conception",
        "components": ", ".join(strict_avail),
        "n_obs": n_alpha_strict,
        "k_items": k_alpha_strict,
        "cronbach_alpha": alpha_strict,
        "pa_factors_retained": None,
    })

    pd.DataFrame(psychometric_results).to_csv(out_dir / "psychometric_diagnostics.csv", index=False)

    # Generate parallel analysis scree plots
    for label, pa_data in [("7comp", pa_7_df), ("12comp", pa_12_df if len(pa_12_df) > 0 else pd.DataFrame())]:
        if len(pa_data) == 0:
            continue
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(pa_data["component"], pa_data["actual_eigenvalue"], "o-", color="#4C72B0",
                label="Actual eigenvalues", markersize=8)
        ax.plot(pa_data["component"], pa_data["pa_95th_threshold"], "s--", color="#C44E52",
                label="95th percentile (random)", markersize=6)
        ax.axhline(1.0, color="gray", linestyle=":", alpha=0.5, label="Kaiser criterion (=1)")
        ax.set_xlabel("Component Number")
        ax.set_ylabel("Eigenvalue")
        ax.set_title(f"Parallel Analysis: PCHA ({label})")
        ax.legend()
        ax.set_xticks(pa_data["component"].values)
        plt.tight_layout()
        fig.savefig(fig_dir / f"fig5_parallel_analysis_{label}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  Figure: parallel analysis {label} saved")

    # ===================================================================
    # ITEM 2: FIX BENEFIT SIGN AND COMPUTE CORRECTED BENEFITS
    # ===================================================================
    print("\n" + "=" * 60)
    print("ITEM 2: CORRECTED PREDICTED-BENEFIT CALCULATIONS")
    print("=" * 60)

    benefit_results = []
    for risk_qid in ["QUO197", "QUO210"]:
        if risk_qid not in ppd_panel.columns:
            continue
        # Rebuild the full PCHA index (same as v1)
        pcha_idx, pcha_avail, pcha_z_cols = build_pcha_index(ppd_panel, PCHA_ALL)
        ppd_panel["pcha_index_full"] = pcha_idx

        mask = (ppd_panel["outcome_ppd"].notna() &
                ppd_panel[risk_qid].notna() &
                ppd_panel["pcha_index_full"].notna())
        work = ppd_panel[mask].copy()
        if len(work) < MIN_N:
            continue

        y = work["outcome_ppd"].values
        risk_z = zscore(work[risk_qid].values)
        pcha_z = zscore(work["pcha_index_full"].values)

        m = fit_ols(y, [risk_z, pcha_z, risk_z * pcha_z],
                    col_names=["risk_z", "pcha_z", "risk_x_pcha"])
        if m:
            b0 = float(m["beta"][0])
            b1 = float(m["beta"][m["col_names"].index("risk_z")])
            b2 = float(m["beta"][m["col_names"].index("pcha_z")])
            b3 = float(m["beta"][m["col_names"].index("risk_x_pcha")])

            # At risk_z=+1, pcha_z=-1: predicted = b0 + b1(1) + b2(-1) + b3(1)(-1)
            pred_high_risk_low_pcha = b0 + b1 - b2 - b3
            # At risk_z=+1, pcha_z=+1: predicted = b0 + b1(1) + b2(+1) + b3(1)(+1)
            pred_high_risk_high_pcha = b0 + b1 + b2 + b3
            # Reduction = low_pcha - high_pcha (positive = PCHA reduces PPD)
            ppd_reduction = pred_high_risk_low_pcha - pred_high_risk_high_pcha

            print(f"  {risk_qid}:")
            print(f"    Predicted PPD at +1SD risk, -1SD PCHA: {pred_high_risk_low_pcha:.2f}%")
            print(f"    Predicted PPD at +1SD risk, +1SD PCHA: {pred_high_risk_high_pcha:.2f}%")
            print(f"    PPD reduction (high vs low PCHA at high risk): {ppd_reduction:.2f} pp")

            benefit_results.append({
                "risk_qid": risk_qid,
                "b0_intercept": b0, "b1_risk": b1, "b2_pcha": b2, "b3_interaction": b3,
                "predicted_ppd_highrisk_lowpcha": pred_high_risk_low_pcha,
                "predicted_ppd_highrisk_highpcha": pred_high_risk_high_pcha,
                "ppd_reduction_pp": ppd_reduction,
                "n": m["n"],
            })

    pd.DataFrame(benefit_results).to_csv(out_dir / "corrected_benefit_estimates.csv", index=False)

    # ===================================================================
    # ITEM 3: STRICT PRE-CONCEPTION PCHA (THE KEY TEST)
    # ===================================================================
    print("\n" + "=" * 60)
    print("ITEM 3: STRICT PRE-CONCEPTION PCHA SENSITIVITY")
    print("=" * 60)

    strict_results = []
    for label, component_dict in [
        ("strict_preconception", PCHA_STRICT_PRECONCEPTION),
        ("narrow_prepreg_only", PCHA_NARROW_PREPREG),
    ]:
        print(f"\n--- {label} ---")
        idx, avail, z_cols = build_pcha_index(ppd_panel, component_dict)
        ppd_panel[f"pcha_{label}"] = idx
        n_valid = int(np.sum(~np.isnan(idx)))
        print(f"  Components: {avail}")
        print(f"  Valid rows: {n_valid}")

        if n_valid < MIN_N:
            print(f"  SKIP: insufficient valid rows ({n_valid} < {MIN_N})")
            continue

        # Cronbach's alpha for this subset
        z_df_sub = pd.DataFrame(z_cols, index=ppd_panel.index)
        alpha_sub, n_sub, k_sub = compute_cronbach_alpha(z_df_sub)
        print(f"  Cronbach's alpha: {alpha_sub:.4f} (k={k_sub}, n={n_sub})")

        for risk_qid in ["QUO197", "QUO210", "QUO313", "QUO315"]:
            if risk_qid not in ppd_panel.columns:
                continue
            mask = (ppd_panel["outcome_ppd"].notna() &
                    ppd_panel[risk_qid].notna() &
                    ppd_panel[f"pcha_{label}"].notna())
            work = ppd_panel[mask].copy()
            if len(work) < MIN_N:
                print(f"  {risk_qid}: n={len(work)} < {MIN_N}, skip")
                continue

            y = work["outcome_ppd"].values
            risk_z = zscore(work[risk_qid].values)
            pcha_z = zscore(work[f"pcha_{label}"].values)
            interaction = risk_z * pcha_z

            # OLS interaction
            m = fit_ols(y, [risk_z, pcha_z, interaction],
                        col_names=["risk_z", "pcha_z", "risk_x_pcha"])

            # FE interaction
            fe_d, fe_n = make_fe_dummies(work)
            m_fe = fit_ols(y, fe_d + [risk_z, pcha_z, interaction],
                           col_names=fe_n + ["risk_z", "pcha_z", "risk_x_pcha"])

            # Cluster-robust
            m_cr = fit_ols_cluster_robust(
                y, fe_d + [risk_z, pcha_z, interaction],
                cluster_ids=work["location_abbr"].values,
                col_names=fe_n + ["risk_z", "pcha_z", "risk_x_pcha"],
            )

            # Permutation test
            print(f"  {risk_qid} × {label}: running {N_PERMUTATIONS} permutations...")
            if m:
                int_idx = m["col_names"].index("risk_x_pcha")
                obs_beta = float(m["beta"][int_idx])
                perm_betas = np.zeros(N_PERMUTATIONS)
                for i in range(N_PERMUTATIONS):
                    pcha_perm = rng.permutation(pcha_z)
                    int_perm = risk_z * pcha_perm
                    m_perm = fit_ols(y, [risk_z, pcha_perm, int_perm],
                                     col_names=["risk_z", "pcha_z", "interaction"])
                    if m_perm:
                        perm_betas[i] = float(m_perm["beta"][3])
                perm_p = float(np.mean(np.abs(perm_betas) >= np.abs(obs_beta)))

                # Bootstrap CI
                boot_betas = np.zeros(N_BOOTSTRAP)
                for i in range(N_BOOTSTRAP):
                    bidx = rng.choice(len(y), size=len(y), replace=True)
                    m_boot = fit_ols(y[bidx], [risk_z[bidx], pcha_z[bidx], interaction[bidx]],
                                     col_names=["risk_z", "pcha_z", "interaction"])
                    if m_boot:
                        boot_betas[i] = float(m_boot["beta"][3])
                    else:
                        boot_betas[i] = obs_beta
                boot_lo = float(np.percentile(boot_betas, 2.5))
                boot_hi = float(np.percentile(boot_betas, 97.5))
            else:
                obs_beta = perm_p = boot_lo = boot_hi = None

            # Within-state permutation (partial permutation)
            if m and risk_qid in ["QUO197", "QUO210"]:
                print(f"    Running within-state permutation ({N_PERMUTATIONS} iter)...")
                within_perm_betas = np.zeros(N_PERMUTATIONS)
                states = work["location_abbr"].values
                for i in range(N_PERMUTATIONS):
                    pcha_within = pcha_z.copy()
                    for st in np.unique(states):
                        st_mask = states == st
                        pcha_within[st_mask] = rng.permutation(pcha_within[st_mask])
                    int_wp = risk_z * pcha_within
                    m_wp = fit_ols(y, [risk_z, pcha_within, int_wp],
                                   col_names=["risk_z", "pcha_z", "interaction"])
                    if m_wp:
                        within_perm_betas[i] = float(m_wp["beta"][3])
                within_perm_p = float(np.mean(np.abs(within_perm_betas) >= np.abs(obs_beta)))
                print(f"    Within-state permutation p: {within_perm_p:.4f}")
            else:
                within_perm_p = None

            # Extract results
            def extract(model, name):
                if model is None:
                    return None, None, None
                if name not in model["col_names"]:
                    return None, None, None
                i = model["col_names"].index(name)
                return float(model["beta"][i]), float(model["se"][i]), float(model["pvals"][i])

            ols_b, ols_se, ols_p = extract(m, "risk_x_pcha")
            fe_b, fe_se, fe_p = extract(m_fe, "risk_x_pcha")
            cr_b, cr_se, cr_p = extract(m_cr, "risk_x_pcha")

            print(f"    OLS: beta={ols_b:.4f}, p={ols_p:.4f}" if ols_b else "    OLS: failed")
            if perm_p is not None:
                print(f"    Perm p={perm_p:.4f}, Boot CI=[{boot_lo:.4f}, {boot_hi:.4f}]")
            if fe_b is not None:
                print(f"    FE:  beta={fe_b:.4f}, p={fe_p:.4f}")
            if cr_b is not None:
                print(f"    CR:  beta={cr_b:.4f}, p={cr_p:.4f}")

            strict_results.append({
                "pcha_variant": label,
                "n_components": len(avail),
                "cronbach_alpha": alpha_sub,
                "risk_qid": risk_qid,
                "n": len(work),
                "ols_interaction_beta": ols_b,
                "ols_interaction_se": ols_se,
                "ols_interaction_p": ols_p,
                "permutation_p": perm_p,
                "within_state_permutation_p": within_perm_p,
                "bootstrap_ci_lo": boot_lo,
                "bootstrap_ci_hi": boot_hi,
                "fe_interaction_beta": fe_b,
                "fe_interaction_p": fe_p,
                "cr_interaction_beta": cr_b,
                "cr_interaction_p": cr_p,
                "adj_r2": m["adj_r2"] if m else None,
                "fe_adj_r2": m_fe["adj_r2"] if m_fe else None,
            })

    strict_df = pd.DataFrame(strict_results)
    strict_df.to_csv(out_dir / "strict_preconception_sensitivity.csv", index=False)

    # ===================================================================
    # ITEM 4: VIF ANALYSIS
    # ===================================================================
    print("\n" + "=" * 60)
    print("ITEM 4: VARIANCE INFLATION FACTORS")
    print("=" * 60)

    vif_results = []
    pcha_idx_full, _, _ = build_pcha_index(ppd_panel, PCHA_ALL)
    ppd_panel["pcha_index_full"] = pcha_idx_full

    for risk_qid in ["QUO197", "QUO210"]:
        if risk_qid not in ppd_panel.columns:
            continue
        mask = (ppd_panel["outcome_ppd"].notna() &
                ppd_panel[risk_qid].notna() &
                ppd_panel["pcha_index_full"].notna())
        work = ppd_panel[mask].copy()
        if len(work) < MIN_N:
            continue

        risk_z = zscore(work[risk_qid].values)
        pcha_z = zscore(work["pcha_index_full"].values)
        interaction = risk_z * pcha_z

        vifs = compute_vif([risk_z, pcha_z, interaction],
                           ["risk_z", "pcha_z", "risk_x_pcha"])

        print(f"  {risk_qid}:")
        for name, vif_val in vifs.items():
            print(f"    VIF({name}) = {vif_val:.2f}")
            vif_results.append({
                "risk_qid": risk_qid,
                "variable": name,
                "vif": vif_val,
            })

    pd.DataFrame(vif_results).to_csv(out_dir / "vif_analysis.csv", index=False)

    # ===================================================================
    # ITEM 5: INVESTIGATE QUO16 ABSENCE
    # ===================================================================
    print("\n" + "=" * 60)
    print("ITEM 5: QUO16 INVESTIGATION")
    print("=" * 60)

    if "QUO16" in panel.columns:
        n_quo16 = int(panel["QUO16"].notna().sum())
        n_quo16_ppd = int(ppd_panel["QUO16"].notna().sum()) if "QUO16" in ppd_panel.columns else 0
        print(f"  QUO16 in full panel: {n_quo16} non-null values")
        print(f"  QUO16 in PPD panel: {n_quo16_ppd} non-null values")
        quo16_note = f"QUO16 present in panel with {n_quo16} values ({n_quo16_ppd} in PPD subset)"
    else:
        # Check if it exists under a different response class
        print("  QUO16 NOT found in panel columns.")
        print("  This likely means QUO16 did not have unstratified YES-response rows")
        print("  in the master dataset, possibly because it uses a different response")
        print("  coding (e.g., 'intended' instead of 'yes').")
        quo16_note = ("QUO16 absent from panel — likely filtered out during panel construction "
                      "because it lacks unstratified YES-class rows in the master data. "
                      "QUO257 (trying to become pregnant) serves as the intentionality proxy.")

    with open(out_dir / "quo16_investigation.txt", "w") as f:
        f.write(quo16_note)
    print(f"  Note saved: {quo16_note[:80]}...")

    # ===================================================================
    # ITEM 6: ERA-SPECIFIC STRICT PRE-CONCEPTION
    # ===================================================================
    print("\n" + "=" * 60)
    print("ITEM 6: ERA-SPECIFIC STRICT PRE-CONCEPTION SENSITIVITY")
    print("=" * 60)

    era_strict_results = []
    for label, component_dict in [("strict_preconception", PCHA_STRICT_PRECONCEPTION)]:
        idx, avail, z_cols = build_pcha_index(ppd_panel, component_dict)

        for era in ["2004-2008", "2009-2011"]:
            era_panel = ppd_panel[ppd_panel["ppd_era"] == era].copy()
            era_idx, era_avail, _ = build_pcha_index(era_panel, component_dict)
            era_panel[f"pcha_{label}"] = era_idx

            for risk_qid in ["QUO197", "QUO210"]:
                if risk_qid not in era_panel.columns:
                    continue
                mask = (era_panel["outcome_ppd"].notna() &
                        era_panel[risk_qid].notna() &
                        era_panel[f"pcha_{label}"].notna())
                work = era_panel[mask].copy()
                if len(work) < MIN_N:
                    print(f"  {era} {risk_qid}: n={len(work)} < {MIN_N}, skip")
                    continue

                y = work["outcome_ppd"].values
                risk_z = zscore(work[risk_qid].values)
                pcha_z = zscore(work[f"pcha_{label}"].values)

                m = fit_ols(y, [risk_z, pcha_z, risk_z * pcha_z],
                            col_names=["risk_z", "pcha_z", "risk_x_pcha"])
                if m:
                    i = m["col_names"].index("risk_x_pcha")
                    print(f"  {era} {risk_qid}: beta={m['beta'][i]:.4f}, p={m['pvals'][i]:.4f}, n={m['n']}")
                    era_strict_results.append({
                        "era": era,
                        "pcha_variant": label,
                        "risk_qid": risk_qid,
                        "n": m["n"],
                        "interaction_beta": float(m["beta"][i]),
                        "interaction_p": float(m["pvals"][i]),
                        "adj_r2": m["adj_r2"],
                    })

    pd.DataFrame(era_strict_results).to_csv(out_dir / "era_strict_preconception.csv", index=False)

    # ===================================================================
    # ITEM 7: COMPARISON FIGURE — FULL vs STRICT vs NARROW
    # ===================================================================
    print("\n" + "=" * 60)
    print("ITEM 7: COMPARISON FIGURE")
    print("=" * 60)

    # Build comparison data from strict_df
    if len(strict_df) > 0:
        # Add the v1 full-PCHA results for comparison
        v1_results = pd.read_csv(
            base_dir / "output" / "run_20260324T_pcha_v1" / "pcha_interaction_models.csv"
        )
        comparison_rows = []
        for _, row in v1_results.iterrows():
            comparison_rows.append({
                "pcha_variant": "full_12comp",
                "risk_qid": row["risk_qid"],
                "ols_interaction_beta": row["interaction_beta"],
                "ols_interaction_p": row["interaction_p_parametric"],
                "permutation_p": row["interaction_p_permutation"],
                "bootstrap_ci_lo": row["interaction_boot_ci_lo"],
                "bootstrap_ci_hi": row["interaction_boot_ci_hi"],
                "fe_interaction_beta": row["fe_interaction_beta"],
                "fe_interaction_p": row["fe_interaction_p"],
                "cr_interaction_beta": row["cr_interaction_beta"],
                "cr_interaction_p": row["cr_interaction_p"],
                "n": row["n"],
            })
        v1_df = pd.DataFrame(comparison_rows)
        compare_df = pd.concat([v1_df, strict_df], ignore_index=True)
        compare_df.to_csv(out_dir / "variant_comparison.csv", index=False)

        # Figure: betas by variant for QUO197 and QUO210
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        plt.rcParams.update({"font.family": "serif", "font.size": 11})

        for ax_idx, risk_qid in enumerate(["QUO197", "QUO210"]):
            ax = axes[ax_idx]
            sub = compare_df[compare_df["risk_qid"] == risk_qid].copy()
            sub = sub[sub["ols_interaction_beta"].notna()]
            if len(sub) == 0:
                continue

            variants = sub["pcha_variant"].values
            betas = sub["ols_interaction_beta"].values
            ci_lo = sub["bootstrap_ci_lo"].fillna(0).values
            ci_hi = sub["bootstrap_ci_hi"].fillna(0).values
            perm_ps = sub["permutation_p"].fillna(1).values

            y_pos = np.arange(len(variants))
            colors = []
            for p in perm_ps:
                if p < 0.05:
                    colors.append("#4C72B0")
                elif p < 0.10:
                    colors.append("#55A868")
                else:
                    colors.append("#8C8C8C")

            ax.barh(y_pos, betas, color=colors, height=0.5, edgecolor="white", alpha=0.8)
            ax.errorbar(betas, y_pos, xerr=[betas - ci_lo, ci_hi - betas],
                       fmt="none", ecolor="black", capsize=3)
            ax.set_yticks(y_pos)
            ax.set_yticklabels([v.replace("_", "\n") for v in variants], fontsize=9)
            ax.axvline(0, color="black", linewidth=0.8)
            ax.set_xlabel("Interaction Beta (Risk × PCHA)")
            ax.set_title(f"{risk_qid}: PCHA Variant Comparison")

            # Annotate with p-values
            for i, (b, p) in enumerate(zip(betas, perm_ps)):
                label = f"perm p={p:.3f}" if np.isfinite(p) and p < 1 else ""
                ax.text(b + 0.02 if b < 0 else b - 0.02, i, label,
                       va="center", ha="left" if b < 0 else "right", fontsize=8)

        plt.tight_layout()
        fig.savefig(fig_dir / "fig7_variant_comparison.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        print("  Figure 7: variant comparison saved")

    # ===================================================================
    # ITEM 8: LOSO FOR STRICT PRE-CONCEPTION
    # ===================================================================
    print("\n" + "=" * 60)
    print("ITEM 8: LOSO FOR STRICT PRE-CONCEPTION PCHA")
    print("=" * 60)

    loso_strict_results = []
    for label, component_dict in [("strict_preconception", PCHA_STRICT_PRECONCEPTION)]:
        idx, avail, _ = build_pcha_index(ppd_panel, component_dict)
        ppd_panel[f"pcha_{label}"] = idx

        for risk_qid in ["QUO197", "QUO210"]:
            if risk_qid not in ppd_panel.columns:
                continue
            mask = (ppd_panel["outcome_ppd"].notna() &
                    ppd_panel[risk_qid].notna() &
                    ppd_panel[f"pcha_{label}"].notna())
            work = ppd_panel[mask].copy()
            if len(work) < MIN_N:
                continue

            states = work["location_abbr"].unique()
            for st in states:
                sub = work[work["location_abbr"] != st]
                if len(sub) < MIN_N:
                    continue
                y = sub["outcome_ppd"].values
                rz = zscore(sub[risk_qid].values)
                pz = zscore(sub[f"pcha_{label}"].values)
                m = fit_ols(y, [rz, pz, rz * pz],
                            col_names=["risk_z", "pcha_z", "risk_x_pcha"])
                if m:
                    i = m["col_names"].index("risk_x_pcha")
                    loso_strict_results.append({
                        "pcha_variant": label,
                        "risk_qid": risk_qid,
                        "excluded_state": st,
                        "n": m["n"],
                        "interaction_beta": float(m["beta"][i]),
                        "interaction_p": float(m["pvals"][i]),
                    })

    loso_strict_df = pd.DataFrame(loso_strict_results)
    loso_strict_df.to_csv(out_dir / "loso_strict_preconception.csv", index=False)
    if len(loso_strict_df) > 0:
        for rq in loso_strict_df["risk_qid"].unique():
            sub = loso_strict_df[loso_strict_df["risk_qid"] == rq]
            betas = sub["interaction_beta"].values
            print(f"  {rq} strict LOSO: range=[{betas.min():.4f}, {betas.max():.4f}], "
                  f"mean={betas.mean():.4f}, SD={betas.std():.4f}")

    # ===================================================================
    # SUMMARY
    # ===================================================================
    print("\n" + "=" * 60)
    print("GENERATING SUMMARY")
    print("=" * 60)

    summary_lines = [
        f"# PCHA Follow-Up Analysis v2 Summary",
        f"",
        f"**Run ID**: {run_id}",
        f"**Timestamp**: {datetime.now(timezone.utc).isoformat()}",
        f"**Parent Run**: run_20260324T_pcha_v1",
        f"**Seed**: {seed}",
        f"**Panel SHA-256**: {panel_hash[:16]}...",
        f"",
        f"## 1. Psychometric Diagnostics",
        f"",
    ]
    for r in psychometric_results:
        summary_lines.append(f"- **{r['subset']}**: alpha={r['cronbach_alpha']:.4f}, "
                           f"PA retained {r['pa_factors_retained']} factors")
    summary_lines.append("")

    summary_lines.extend([
        f"## 2. Corrected Benefit Estimates",
        f"",
    ])
    for r in benefit_results:
        summary_lines.append(
            f"- {r['risk_qid']}: At +1SD risk, moving from low to high PCHA reduces "
            f"PPD by **{r['ppd_reduction_pp']:.2f} pp** "
            f"({r['predicted_ppd_highrisk_lowpcha']:.1f}% → {r['predicted_ppd_highrisk_highpcha']:.1f}%)"
        )
    summary_lines.append("")

    summary_lines.extend([
        f"## 3. Strict Pre-Conception PCHA Results",
        f"",
    ])
    if len(strict_df) > 0:
        for _, r in strict_df.iterrows():
            sig = "***" if (r.get("ols_interaction_p") or 1) < 0.05 else ""
            perm_str = f"{r['permutation_p']:.4f}" if r['permutation_p'] is not None and not pd.isna(r['permutation_p']) else "N/A"
            summary_lines.append(
                f"- {r['pcha_variant']} x {r['risk_qid']} (n={int(r['n'])}): "
                f"beta={r['ols_interaction_beta']:.4f}, p={r['ols_interaction_p']:.4f}, "
                f"perm p={perm_str} {sig}"
            )
    summary_lines.append("")

    summary_lines.extend([
        f"## 4. VIF Results",
        f"",
    ])
    for r in vif_results:
        flag = " (HIGH)" if r["vif"] > 5 else ""
        summary_lines.append(f"- {r['risk_qid']} {r['variable']}: VIF={r['vif']:.2f}{flag}")
    summary_lines.append("")

    summary_lines.extend([
        f"## 5. QUO16 Investigation",
        f"",
        f"- {quo16_note}",
        f"",
    ])

    with open(out_dir / "v2_summary.md", "w") as f:
        f.write("\n".join(summary_lines))

    metadata["outputs"] = [str(p.name) for p in out_dir.iterdir()]
    with open(out_dir / "run_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)

    print(f"\n{'=' * 60}")
    print(f"ANALYSIS COMPLETE")
    print(f"Run ID: {run_id}")
    print(f"Output: {out_dir}")
    print(f"Figures: {fig_dir}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
