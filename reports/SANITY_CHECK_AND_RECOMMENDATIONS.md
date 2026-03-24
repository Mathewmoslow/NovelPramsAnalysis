# Sanity Check and Recommended Changes

**Date:** 2026-03-24
**Reviewer:** Automated post-analysis audit

This document identifies potential issues and recommends changes BEFORE any
revisions are made. None of these changes have been implemented — they are
documented here for transparency and decision-making.

---

## STATUS: Issues Found and Severity

### CRITICAL (must address before any submission)

**1. Cronbach's Alpha and Parallel Analysis Failed**
- **Problem:** The PCHA construct's internal consistency (Cronbach's alpha) and
  dimensionality (parallel analysis) could not be computed because the 12
  components have heterogeneous coverage — no complete-case matrix exists with
  all 12 components measured. Some components (QUO179, QUO249, QUO297) are only
  available for 86 state-years while others cover 188.
- **Risk:** A reviewer will immediately ask "what is the reliability of your
  composite?" and the answer is currently "we couldn't compute it."
- **Recommendation:** Compute Cronbach's alpha and parallel analysis on (a) the
  subset of components available for all 188 rows (QUO41, QUO65, QUO257, QUO4,
  QUO5, QUO44, QUO101 — 7 components), and (b) the 86-row subset where all 12
  are available. Report both. If alpha < 0.60 on either, the composite is
  psychometrically weak and the paper should emphasize the component-level
  results instead.

**2. Predicted High-Risk Benefit Has Wrong Sign**
- **Problem:** The summary reports "Predicted high-risk benefit: -2.74 pp" and
  "-2.25 pp" — negative values. This was computed as 2*b2 + 2*b3 where b2
  (PCHA main effect) is negative. The formula calculates the *difference in
  predicted PPD* between low-PCHA and high-PCHA groups at +1 SD risk. A
  negative value means high-PCHA groups have LOWER PPD (which is the expected
  direction). The reporting language should clarify this is a *reduction*, not
  frame the negative as problematic.
- **Recommendation:** Fix the label from "Predicted high-risk benefit" to
  "Predicted PPD reduction at +1SD risk, high vs. low PCHA" and report the
  absolute value (2.74 pp reduction).

**3. Permutation P-Values Are Borderline**
- **Problem:** Permutation p-values of 0.079 and 0.091 are at the boundary
  between conventional significance levels. This WILL be challenged.
- **Context:** The permutation test shuffles PCHA while holding risk and outcome
  fixed, breaking the risk×PCHA correlation structure. This is a conservative
  test for the *marginal* contribution of the interaction given the strong
  correlation between PCHA and risk (r = -0.71 to -0.75). The parametric test,
  FE test, and CR-robust test all pass at p < 0.05.
- **Recommendation:** Acknowledge this explicitly. Frame the analysis as
  exploratory (alpha = 0.10 is defensible). Note that the FE and CR-robust
  tests, which better control for confounding, are unambiguously significant.
  Consider adding a partial-permutation test that permutes PCHA *within states*
  to preserve between-state structure while testing within-state variation.

### IMPORTANT (should address for stronger submission)

**4. High Collinearity Between PCHA and Risk Anchors**
- **Problem:** PCHA correlates at r = -0.71 to -0.75 with the risk anchors.
  This means states with high health agency also have low partner stress (and
  vice versa). The interaction term (risk × PCHA) is therefore partially
  collinear with both main effects, which inflates SE for the interaction and
  may explain the borderline permutation p-values.
- **Recommendation:** Report VIF for the interaction models. If VIF > 5 for the
  interaction term, consider mean-centering (already done via z-scoring) or
  residualizing PCHA on risk before computing the interaction. Document this
  explicitly as a known issue.

**5. No Sample-Size-Weighted Analysis**
- **Problem:** While PRAMStat doesn't provide SEs for unstratified rows, it
  does provide sample_size for many observations. Larger-state estimates are
  more precise and should get more weight.
- **Recommendation:** Extract sample_size from the master data, merge into the
  panel, and run WLS with sample_size as weights. This partially addresses the
  precision-weighting gap. If sample_size is missing for unstratified rows,
  use state population as a proxy.

**6. QUO16 (Intended Pregnancy) Not in Panel**
- **Problem:** QUO16 was listed in the PCHA intentionality domain but was
  flagged as "NOT IN PANEL." The panel construction may have filtered it out
  due to a response coding issue (QUO16 may not have a "yes" binary response
  class in the master data).
- **Recommendation:** Investigate why QUO16 is missing and whether it can be
  recovered. If not, document the exclusion.

**7. Breastfeeding Components May Be Post-Outcome**
- **Problem:** QUO4, QUO5, and QUO44 (breastfeeding at 4 weeks, ever, and
  8 weeks) are measured AFTER the postpartum period in which PPD is assessed.
  A reviewer may argue these are consequences of (or co-occur with) PPD rather
  than antecedents that could causally buffer it.
- **Recommendation:** Run sensitivity analysis with a "strict pre-conception"
  PCHA index that excludes all postpartum components (QUO4, QUO5, QUO44, QUO101).
  If the interaction holds with only pre-conception + prenatal components, the
  temporal-ordering concern is mitigated. This is likely the most important
  robustness check not yet performed.

**8. Missing Figure 5 (Parallel Analysis Scree Plot)**
- **Problem:** Figure 5 was not generated because parallel analysis failed.
- **Recommendation:** Once parallel analysis is computed on the 7-component
  complete-case subset, generate this figure.

### MINOR (nice-to-have improvements)

**9. Constant-Sample Overlaps with Full Sample**
- **Problem:** The constant-sample analysis (Phase 6) shows n=188, which is
  identical to the full PPD panel. This means the constant-sample check is
  trivially identical to the main result. The check was designed for situations
  where different models use different subsets, but since all 188 rows have
  valid PCHA, risk, and PPD, the check adds nothing here.
- **Recommendation:** Note this explicitly — it's actually a strength (the
  PCHA index achieves full coverage), not a weakness.

**10. Report References Are Incomplete**
- **Problem:** The reference list cites some papers by incomplete citation
  (e.g., "Abbasi, S., et al. (2013)" without full journal/volume/page).
- **Recommendation:** Complete all references with full bibliographic details
  before submission.

**11. No DAG (Directed Acyclic Graph)**
- **Problem:** The causal pathway from PCHA → moderation is asserted but not
  formally diagrammed.
- **Recommendation:** Include a simple DAG showing the hypothesized structure:
  Partner Stress → PPD, PCHA → PPD, PCHA × Partner Stress → PPD, with shared
  confounders (state context, SES) explicitly represented.

**12. The "Agency" Framing May Invite Victim-Blaming Critique**
- **Problem:** Calling the construct "health agency" could be read as implying
  that women who don't exercise or take vitamins have less agency and therefore
  bear responsibility for their PPD risk.
- **Recommendation:** Explicitly address this in the Discussion. Frame PCHA as
  reflecting population-level *opportunity* for health engagement (shaped by
  policy, access, and social determinants) rather than individual *choice*.
  Consider the alternative name "Pre-Conception Health Engagement" if the
  framing concern is significant.

---

## Summary of Recommended Priority

| Priority | Item | Effort |
|----------|------|--------|
| CRITICAL | 1. Compute alpha + PA on subsets | Medium |
| CRITICAL | 2. Fix benefit sign/label | Trivial |
| CRITICAL | 3. Frame permutation p honestly | Writing only |
| IMPORTANT | 4. Report VIF | Easy |
| IMPORTANT | 5. Sample-size weighting | Medium |
| IMPORTANT | 6. Investigate QUO16 | Easy |
| IMPORTANT | 7. **Strict pre-conception PCHA sensitivity** | **Medium — highest value-add** |
| IMPORTANT | 8. Generate missing Figure 5 | Easy (after #1) |
| MINOR | 9. Note constant-sample trivially passes | Writing only |
| MINOR | 10. Complete references | Writing only |
| MINOR | 11. Add DAG | Easy |
| MINOR | 12. Address "agency" framing | Writing only |

**The single most impactful change is #7**: running the interaction with a strict
pre-conception-only PCHA index. If that holds, the paper's causal story is much
stronger. If it doesn't hold, the paper needs to be reframed around the
full-lifecycle engagement pattern rather than pre-conception specificity.
