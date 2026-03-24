# PCHA Follow-Up Analysis v2 Summary

**Run ID**: run_20260324T_pcha_v2
**Timestamp**: 2026-03-24T22:17:39.844114+00:00
**Parent Run**: run_20260324T_pcha_v1
**Seed**: 42
**Panel SHA-256**: db32d98dff8d81b1...

## 1. Psychometric Diagnostics

- **7-component full-coverage (n=188)**: alpha=0.8711, PA retained 2 factors
- **12-component complete-case**: alpha=nan, PA retained 0 factors
- **8-component strict pre-conception**: alpha=nan, PA retained None factors

## 2. Corrected Benefit Estimates

- QUO197: At +1SD risk, moving from low to high PCHA reduces PPD by **2.74 pp** (14.9% → 12.2%)
- QUO210: At +1SD risk, moving from low to high PCHA reduces PPD by **2.25 pp** (14.9% → 12.7%)

## 3. Strict Pre-Conception PCHA Results

- strict_preconception x QUO197 (n=188): beta=-0.3871, p=0.0021, perm p=0.0894 ***
- strict_preconception x QUO210 (n=188): beta=-0.3921, p=0.0009, perm p=0.0628 ***
- strict_preconception x QUO313 (n=111): beta=-0.0967, p=0.4376, perm p=0.7164 
- strict_preconception x QUO315 (n=111): beta=-0.0561, p=0.6603, perm p=0.8084 
- narrow_prepreg_only x QUO197 (n=188): beta=-0.3796, p=0.0025, perm p=0.1012 ***
- narrow_prepreg_only x QUO210 (n=188): beta=-0.3833, p=0.0013, perm p=0.0722 ***
- narrow_prepreg_only x QUO313 (n=111): beta=-0.0891, p=0.4735, perm p=0.7308 
- narrow_prepreg_only x QUO315 (n=111): beta=-0.0500, p=0.6916, perm p=0.8288 

## 4. VIF Results

- QUO197 risk_z: VIF=2.40
- QUO197 pcha_z: VIF=2.00
- QUO197 risk_x_pcha: VIF=1.35
- QUO210 risk_z: VIF=2.63
- QUO210 pcha_z: VIF=2.27
- QUO210 risk_x_pcha: VIF=1.32

## 5. QUO16 Investigation

- QUO16 absent from panel — likely filtered out during panel construction because it lacks unstratified YES-class rows in the master data. QUO257 (trying to become pregnant) serves as the intentionality proxy.
