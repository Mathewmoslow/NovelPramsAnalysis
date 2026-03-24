# PCHA Analysis Summary

**Run ID**: run_20260324T_pcha_v1
**Timestamp**: 2026-03-24T22:01:09.127624+00:00
**Random Seed**: 42
**Panel SHA-256**: db32d98dff8d81b1...

## Data
- Panel source: /Users/mathewmoslow/Documents/AHU/PPD_Analysis/NovelPramsAnalysis/data/analysis_panel.csv
- Panel: 329 state-year rows × 233 columns
- PPD panel: 188 rows with outcome
- Era breakdown: {'2004-2008': 102, '2009-2011': 86}

## PCHA Index
- Components used: 12
  - QUO179, QUO41, QUO65, QUO249, QUO75, QUO257, QUO296, QUO297, QUO4, QUO5, QUO44, QUO101
- Valid rows: 188
- Cronbach's alpha: nan
- Mean: -0.011, SD: 0.726

## FE-Only Baseline
- R²: 0.8043
- Adj-R²: 0.7476
- n: 188, k: 43

## Parallel Analysis
- Retained components: 0

## Interaction Results

### QUO197: Argued with husband/partner more than usual (12 mo before pregnancy)
- n = 188
- Risk-only adj-R²: 0.3717
- Additive adj-R²: 0.4207
- Interaction adj-R²: 0.4497
- **Interaction beta: -0.3819**
- Parametric p: 0.0012
- **Permutation p: 0.0914** (5000 iterations)
- Bootstrap 95% CI: [-0.6502, -0.0421]
- Attenuation signature: True
- Predicted high-risk benefit: -2.74 pp
- FE interaction: beta=-0.4741, p=0.0141
- CR-robust: beta=-0.4741, p=0.0399

### QUO210: Any partner-related stressors reported
- n = 188
- Risk-only adj-R²: 0.4199
- Additive adj-R²: 0.4433
- Interaction adj-R²: 0.4720
- **Interaction beta: -0.3735**
- Parametric p: 0.0011
- **Permutation p: 0.0790** (5000 iterations)
- Bootstrap 95% CI: [-0.6362, -0.0823]
- Attenuation signature: True
- Predicted high-risk benefit: -2.25 pp
- FE interaction: beta=-0.5476, p=0.0043
- CR-robust: beta=-0.5476, p=0.0306

### QUO313: Physical IPV by ex-partner (before pregnancy)
- n = 111
- Risk-only adj-R²: 0.4709
- Additive adj-R²: 0.5803
- Interaction adj-R²: 0.5790
- **Interaction beta: -0.0917**
- Parametric p: 0.4139
- **Permutation p: 0.7276** (5000 iterations)
- Bootstrap 95% CI: [-0.9655, 0.1457]
- Attenuation signature: True
- Predicted high-risk benefit: -2.79 pp
- FE interaction: beta=0.2384, p=0.4117
- CR-robust: beta=0.2384, p=0.4687

### QUO315: Physical IPV by ex-partner (during pregnancy)
- n = 111
- Risk-only adj-R²: 0.4990
- Additive adj-R²: 0.5946
- Interaction adj-R²: 0.5928
- **Interaction beta: -0.0825**
- Parametric p: 0.4717
- **Permutation p: 0.7160** (5000 iterations)
- Bootstrap 95% CI: [-0.6749, 0.3672]
- Attenuation signature: True
- Predicted high-risk benefit: -2.61 pp
- FE interaction: beta=-0.1396, p=0.5881
- CR-robust: beta=-0.1396, p=0.5418