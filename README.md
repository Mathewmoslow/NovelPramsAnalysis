# Pre-Conception Health Agency (PCHA) and Postpartum Depression: A Novel Ecological Analysis

## Overview

This repository contains a fully reproducible ecological panel analysis testing whether
a latent "Pre-Conception Health Agency" (PCHA) construct — capturing pre-pregnancy
exercise, multivitamin use, dental care, pregnancy intentionality, and early prenatal
care engagement — moderates the well-established partner-stress → postpartum depressive
symptom (PPD) pathway at the US state-year level.

## Data Source

CDC PRAMStat (Pregnancy Risk Assessment Monitoring System Statistics), 2000–2011.
State-year aggregate prevalence estimates for 141+ behavioral and health indicators
across 40 US locations.

The master dataset is built from publicly available PRAMStat downloads via the
parent project pipeline. See `data/DATA_PROVENANCE.md` for full chain of custody.

## Reproduction

```bash
# Requirements: Python 3.10+, pandas, numpy, scipy, matplotlib
pip install pandas numpy scipy matplotlib

# Run the full analysis (deterministic, seeded)
python scripts/run_pcha_analysis.py

# All outputs land in output/<run_id>/
# Figures land in figures/
# Report is generated in reports/
```

## Directory Structure

```
NovelPramsAnalysis/
├── README.md                  # This file
├── data/                      # Input data + provenance documentation
│   ├── analysis_master.csv    # Symlink or copy of master dataset
│   └── DATA_PROVENANCE.md     # Full data lineage documentation
├── scripts/                   # All analysis code
│   └── run_pcha_analysis.py   # Main reproducible analysis script
├── output/                    # Timestamped run outputs (CSV, JSON)
├── figures/                   # All generated figures
└── reports/                   # Final report and supplementary materials
```

## Key Research Questions

1. Does a theory-driven PCHA latent construct exist in PRAMStat ecological data?
2. Does PCHA moderate the partner-stress → PPD pathway?
3. Is this moderation distinct from insurance/structural access (LV3/LV6)?
4. Does the finding survive constant-sample, precision-weighted, and permutation tests?

## Addressing Prior Board Critique

This analysis explicitly addresses all 9 criticisms from the prior rigor review:
- Constant-sample sensitivity analyses
- Fixed-effects-only baseline R² reported
- Inverse-variance (precision) weighting
- Full variable dictionary with question wording
- Cluster-robust standard errors
- Parallel analysis for factor retention
- Complete QUO code documentation
- Transparent forking-path documentation
- Split-era sensitivity checks

## License

Academic use. Data sourced from CDC public domain PRAMStat system.
