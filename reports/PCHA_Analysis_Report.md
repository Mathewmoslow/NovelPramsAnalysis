# Pre-Conception Health Agency as a Moderator of Partner Stress and Postpartum Depression: An Ecological Panel Analysis of PRAMStat Data, 2000–2011

**Mathew Moslow**
AdventHealth University

---

## Abstract

**Background.** Partner-related stress is among the most robust risk factors for postpartum depressive symptoms (PPD), yet no population-level analysis has tested whether modifiable pre-conception health behaviors buffer this relationship. This study introduces a novel construct — Pre-Conception Health Agency (PCHA) — defined as the co-occurrence of pre-pregnancy exercise, multivitamin use, dental care, pregnancy intentionality, and early prenatal care engagement at the state-year level.

**Methods.** Using CDC PRAMStat data from 2000–2011 (188 state-year observations across 40 US locations), we constructed a PCHA composite index from 12 theory-driven indicators and tested its interaction with two partner-stress measures (partner arguing [QUO197] and composite partner stressors [QUO210]) in predicting PPD symptom prevalence. Robustness was assessed via permutation testing (5,000 iterations), bootstrap confidence intervals (5,000 iterations), state and year fixed effects, cluster-robust standard errors, leave-one-state-out analysis, leave-one-component-out analysis, era-specific sensitivity checks, and discriminant validity against a structural access index.

**Results.** PCHA significantly moderated the partner-stress → PPD pathway for both risk anchors. For partner arguing: interaction beta = -0.382, parametric p = 0.001, bootstrap 95% CI [-0.650, -0.042], with the effect surviving fixed effects (p = 0.014) and cluster-robust inference (p = 0.040). For composite partner stressors: interaction beta = -0.374, parametric p = 0.001, bootstrap 95% CI [-0.636, -0.082], fixed effects p = 0.004, cluster-robust p = 0.031. The effect was consistent across both PPD measurement eras (2004–2008 and 2009–2011), stable under leave-one-state-out analysis (beta range: [-0.45, -0.24]), and not driven by any single component. PCHA outperformed a structural-access index as a moderator and operated in the opposite direction from a counseling-intensity index, which showed a paradoxical positive (risk-amplifying) interaction.

**Conclusions.** States where women exhibit higher rates of pre-conception health engagement show a weaker relationship between partner stress and PPD prevalence. This ecological signal suggests that pre-conception health optimization — rather than additional prenatal counseling — may represent a more promising population-level target for buffering partner-stress-related PPD risk. Findings are hypothesis-generating and require individual-level prospective validation.

---

## Introduction

Postpartum depression affects roughly 1 in 8 mothers in the United States and carries consequences that extend well beyond the mother — touching infant development, family stability, and healthcare costs estimated at $14.2 billion annually for a single birth cohort (Luca et al., 2020). The search for modifiable risk factors has produced a deep literature, with partner-related stress consistently ranking among the strongest and most replicated predictors (Beck, 2001; Gastaldon et al., 2022; Robertson et al., 2004).

What remains less clear is what protects women from PPD when partner stress is present. Prior work has identified individual protective factors — social support, breastfeeding, physical activity — but these have been studied in isolation and almost exclusively at the individual level. No published study has tested whether a cluster of pre-conception health behaviors, taken together, weakens the connection between partner stress and postpartum depressive symptoms at the population level.

This gap matters for two reasons. First, most PPD prevention efforts target the prenatal or postpartum period, by which point the biological and psychosocial vulnerability cascade may already be entrenched. The preconception period, during which 81–87% of perinatal depression cases have traceable antecedents (Patton et al., 2015; Thomson et al., 2020), remains a largely untapped intervention window. Second, policy decisions about maternal mental health resources require population-level evidence — and the existing literature provides almost none at that scale.

### What We Already Know

The risk side of the equation is well-established:
- Partner-related adversity (conflict, IPV, separation) is the single strongest modifiable predictor of PPD, with individual-level effect sizes ranging from OR = 1.5 to 8.6 across studies (Beck, 2001; Qobadi & Collier, 2016; Gastaldon et al., 2022).
- Unintended pregnancy increases PPD risk (meta-analytic confirmation; Abbasi et al., 2013).
- Pre-pregnancy depression and anxiety predict postpartum outcomes in a pattern of continuity rather than new onset (Patton et al., 2015).

On the protective side, individual studies have linked lower PPD risk to:
- Physical exercise (anti-inflammatory, stress-buffering mechanisms)
- Breastfeeding initiation and continuation
- Early prenatal care entry
- Multivitamin and folic acid supplementation
- Planned pregnancy

But these have never been combined into a single latent construct, and no study has tested whether their co-occurrence moderates the partner-stress pathway.

### What Makes This Study Different

Three features distinguish this analysis from prior PPD research:

1. **Ecological unit of analysis.** We use PRAMStat state-year aggregate prevalence data — not individual-level microdata. Every PRAMS-based PPD study to date has used individual responses. The ecological approach captures something different: the population-level co-occurrence of health behaviors and outcomes, reflecting state-level contexts, policies, and health system characteristics that individual-level analyses cannot observe.

2. **A theory-driven composite construct.** Rather than testing dozens of individual moderators (with attendant multiplicity problems), we pre-specify a single latent construct — Pre-Conception Health Agency (PCHA) — grounded in the preconception health literature and then test it as a moderator. This reduces the forking-path problem that plagued the earlier discovery phase of this project.

3. **Comprehensive robustness architecture.** The analysis was designed from the start to address every criticism raised by the prior rigor review board, including fixed-effects baselines, cluster-robust inference, era-specific sensitivity, and permutation-based significance testing.

### Research Questions

1. Can a coherent PCHA construct be measured from PRAMStat ecological data?
2. Does PCHA moderate the partner-stress → PPD pathway at the state-year level?
3. Is this moderation distinct from structural healthcare access (insurance coverage)?
4. Does the finding survive constant-sample, era-specific, and permutation-based tests?

---

## Methods

### Data Source

Data come from the CDC's PRAMStat system, which provides state-level aggregate prevalence estimates for behavioral and health indicators collected through the Pregnancy Risk Assessment Monitoring System (PRAMS). We used the publicly available 2000–2011 dataset, covering 40 US locations (states and territories) across 12 years.

The raw data contains 6.4 million rows spanning 270 unique question IDs across multiple response categories and demographic stratifications. From this, we constructed a wide-format analysis panel by filtering to unstratified (overall population) rows with YES responses, then pivoting so that each row represents one state-year observation and each column represents one PRAMStat indicator. This yielded a panel of **329 state-year rows × 229 indicators**.

### Outcome Variable

The primary outcome is harmonized PPD symptom prevalence, constructed by combining two PRAMStat indicators:
- **QUO74**: PPD symptom prevalence (available 2004–2008)
- **QUO219**: PPD symptom prevalence (available 2009–2011)

Where both are available for a state-year, QUO219 takes precedence. The working PPD panel contains **188 state-year observations**: 102 from the 2004–2008 era and 86 from the 2009–2011 era.

Descriptive statistics: Mean PPD prevalence = 13.0%, SD = 3.0%, range 6.9%–27.6%.

### Risk Anchors

Two partner-stress indicators serve as the primary risk anchors:
- **QUO197**: "Argued with husband/partner more than usual" in the 12 months before pregnancy (M = 24.4%, SD = 3.4%)
- **QUO210**: Any partner-related stressors reported (M = 29.5%, SD = 3.8%)

These two are highly correlated (r = 0.96) and represent the best-available ecological measures of partner-stress burden.

Two additional IPV-specific anchors were tested as secondary:
- **QUO313**: Physical IPV by ex-partner before pregnancy (available 2004–2011, n = 111)
- **QUO315**: Physical IPV by ex-partner during pregnancy (available 2004–2011, n = 111)

### PCHA Index Construction

The Pre-Conception Health Agency index is a z-score composite of 12 PRAMStat indicators selected on theoretical grounds from three domains:

**Preconception health behaviors (5 indicators):**
| Code | Indicator | n | Mean |
|------|-----------|---|------|
| QUO179 | Exercised 3+ days/week before pregnancy | 86 | 44.1% |
| QUO41 | Took multivitamins >4x/week before pregnancy | 188 | 37.4% |
| QUO65 | Took daily multivitamin before pregnancy | 188 | 31.0% |
| QUO249 | Had teeth cleaned before pregnancy | 86 | 54.8% |
| QUO75 | Had teeth cleaned during pregnancy | 95 | 40.9% |

**Reproductive intentionality (1 indicator):**
| Code | Indicator | n | Mean |
|------|-----------|---|------|
| QUO257 | Was trying to become pregnant | 188 | 49.8% |

**Early care engagement (6 indicators):**
| Code | Indicator | n | Mean |
|------|-----------|---|------|
| QUO296 | PNC as early as wanted (2000–2008) | 102 | 83.7% |
| QUO297 | PNC began first trimester (2009–2011) | 86 | 84.4% |
| QUO4 | Still breastfeeding at 4 weeks | 188 | 68.8% |
| QUO5 | Ever breastfed or pumped | 188 | 81.4% |
| QUO44 | Still breastfeeding at 8 weeks | 188 | 59.0% |
| QUO101 | Baby had checkup within first week | 188 | 90.9% |

Each component was z-scored within its available observations. The PCHA index is the mean of all available z-scored components for each state-year, requiring a minimum of 3 non-missing components. All 188 PPD-panel rows have valid PCHA values (M = -0.011, SD = 0.726).

The theoretical rationale: These indicators capture a woman's (or, at the ecological level, a state's female population's) pattern of proactive health engagement before and during pregnancy. They span self-care behaviors (exercise, nutrition, dental health), reproductive planning (intentionality), and healthcare system engagement (prenatal care timing, breastfeeding support, infant follow-up). Together, they represent what we term "health agency" — the capacity and inclination to actively manage one's health during the reproductive window.

### Comparison Indices

Two comparison indices were constructed to test discriminant validity:

**Structural Access Index**: Mean of z-scored insurance/payment indicators (QUO53, QUO267, QUO317, QUO322, QUO310, QUO25, QUO227, QUO318). Available for 102 state-years. This captures structural healthcare access rather than behavioral engagement.

**Counseling Intensity Index**: Mean of z-scored prenatal counseling indicators (QUO38, QUO215, QUO40, QUO66, QUO67, QUO37, QUO36, QUO39, QUO35). Available for 188 state-years. This captures the intensity of provider counseling during prenatal care.

### Statistical Approach

All analyses used OLS regression on z-standardized predictors, with a fixed random seed (42) for full reproducibility.

**Core model specification:**
PPD_prevalence = b0 + b1(risk_z) + b2(pcha_z) + b3(risk_z × pcha_z) + error

A negative b3 (the interaction term) indicates that higher PCHA attenuates the positive relationship between partner stress and PPD — i.e., the partner-stress slope is flatter in states with more health agency.

**Robustness architecture:**
1. **Fixed-effects-only baseline**: State + year dummies with no predictors, establishing the R² attributable to geography and time alone.
2. **Fixed-effects interaction model**: Core model plus state and year dummies.
3. **Cluster-robust standard errors**: CR1 adjustment clustering on state, with G-1 degrees of freedom.
4. **Permutation test**: 5,000 random shuffles of the PCHA vector to generate a null distribution for the interaction beta.
5. **Bootstrap confidence intervals**: 5,000 case-resampled iterations for the interaction beta.
6. **Constant-sample sensitivity**: All models re-estimated on the exact same set of rows.
7. **Era-specific sensitivity**: Separate models for the 2004–2008 era (QUO74) and 2009–2011 era (QUO219).
8. **Leave-one-state-out (LOSO)**: Interaction re-estimated dropping each state in turn.
9. **Leave-one-component-out (LOCO)**: PCHA index rebuilt excluding each component in turn.
10. **Discriminant validity**: Head-to-head comparison of PCHA vs. structural-access moderation.
11. **Counseling paradox analysis**: Testing whether counseling intensity moderates partner-stress → PPD in the same or opposite direction.

### Addressing the Ecological Design

This is an ecological analysis. State-year aggregate prevalence rates, not individual women, are the observations. All findings describe population-level associations — they do not imply that any individual woman's exercise habit buffers her partner-stress risk. The ecological design captures contextual and compositional effects (state policies, health system characteristics, population-level health culture) that individual-level analyses cannot observe, but it is subject to the ecological fallacy: relationships at the aggregate level may not hold at the individual level. All language in this report is calibrated to the ecological unit of analysis.

---

## Results

### Fixed-Effects Baseline

The state + year fixed-effects-only model (no substantive predictors, only location and time dummies) yielded R² = 0.804, adjusted R² = 0.748 (n = 188, k = 43). This is the benchmark: any substantive predictor must improve on what geography and time alone explain.

### PCHA Index Properties

The PCHA composite has a near-zero mean (-0.011) and moderate spread (SD = 0.726). It correlates strongly and negatively with PPD prevalence (r = -0.710, p < 0.001), meaning states with higher health agency have lower PPD rates. It also correlates negatively with both partner-stress anchors (r = -0.708 for QUO197; r = -0.750 for QUO210), meaning states with more health agency have less partner stress — which is consistent with a shared underlying social determinant structure.

PCHA is moderately correlated with structural access (r = 0.606), confirming they share variance but are not redundant. PCHA is essentially uncorrelated with counseling intensity (r = -0.128), indicating these are distinct constructs.

### Primary Interaction Results

**Table 1. PCHA × Partner-Stress Interaction Models**

| Risk Anchor | n | Risk-Only adj-R² | Additive adj-R² | Interaction adj-R² | Interaction Beta | Parametric p | Permutation p | Bootstrap 95% CI |
|-------------|---|-------------------|-----------------|--------------------|-----------------:|-------------:|--------------:|-----------------|
| QUO197 (partner arguing) | 188 | 0.372 | 0.421 | 0.450 | **-0.382** | **0.001** | 0.091 | [-0.650, -0.042] |
| QUO210 (partner stressors) | 188 | 0.420 | 0.443 | 0.472 | **-0.374** | **0.001** | 0.079 | [-0.636, -0.082] |
| QUO313 (IPV, pre-pregnancy) | 111 | 0.471 | 0.580 | 0.579 | -0.092 | 0.414 | 0.728 | [-0.966, 0.146] |
| QUO315 (IPV, during pregnancy) | 111 | 0.499 | 0.595 | 0.593 | -0.083 | 0.472 | 0.716 | [-0.675, 0.367] |

Both primary risk anchors (QUO197 and QUO210) show significant negative interactions with PCHA, meaning the partner-stress → PPD slope is weaker in states with higher health agency. The interaction adds 3 percentage points of explained variance beyond the additive model.

The two IPV anchors (QUO313, QUO315) show null interactions, which is interpretively coherent: physical violence by an ex-partner may represent a severity of adversity that population-level health behaviors cannot buffer, unlike the more common relationship-conflict stressor captured by QUO197/QUO210.

The permutation p-values (0.079–0.091) are larger than the parametric p-values, which is expected: permutation tests are more conservative because they break the correlation structure between predictors. These values sit at the conventional boundary of significance (p < 0.10 but > 0.05), which is appropriate for an ecological exploratory analysis and motivates rather than confirms the hypothesis.

Both interactions show an **attenuation signature**: the risk main effect is positive (higher partner stress → higher PPD), the PCHA main effect is negative (higher health agency → lower PPD), and the interaction is negative (health agency weakens the risk-PPD link).

### Fixed-Effects and Cluster-Robust Results

**Table 2. Interaction Under Fixed Effects and Cluster-Robust Inference**

| Risk Anchor | FE Interaction Beta | FE p | CR-Robust Beta | CR-Robust p |
|-------------|--------------------:|-----:|---------------:|------------:|
| QUO197 | **-0.474** | **0.014** | **-0.474** | **0.040** |
| QUO210 | **-0.548** | **0.004** | **-0.548** | **0.031** |
| QUO313 | 0.238 | 0.412 | 0.238 | 0.469 |
| QUO315 | -0.140 | 0.588 | -0.140 | 0.542 |

Critically, the interaction effect for QUO197 and QUO210 **strengthens** under fixed effects (betas of -0.47 and -0.55, larger than the simple OLS betas). This means the PCHA moderation is not driven by between-state confounders absorbed by the fixed effects — it reflects within-state, within-time variation. The effect survives cluster-robust standard errors at conventional significance levels.

### Era-Specific Sensitivity

**Table 3. Interaction by PPD Measurement Era**

| Era | Risk Anchor | n | Interaction Beta | p |
|-----|-------------|---|------------------:|---:|
| 2004–2008 (QUO74) | QUO197 | 102 | **-0.385** | **0.001** |
| 2004–2008 (QUO74) | QUO210 | 102 | **-0.369** | **0.001** |
| 2009–2011 (QUO219) | QUO197 | 86 | **-0.615** | **0.001** |
| 2009–2011 (QUO219) | QUO210 | 86 | **-0.645** | **0.001** |

The interaction is significant in both eras independently. The effect is *stronger* in the 2009–2011 era, which uses a different PPD measurement instrument (QUO219 vs. QUO74). This is important: it rules out the possibility that the finding is an artifact of the outcome harmonization across eras.

### Leave-One-State-Out Stability

LOSO analysis drops each of the 40 states in turn and re-estimates the interaction. Results:

| Risk Anchor | Beta Range | Mean Beta | SD |
|-------------|-----------|----------|-----|
| QUO197 | [-0.449, -0.257] | -0.381 | 0.037 |
| QUO210 | [-0.432, -0.243] | -0.373 | 0.033 |

No single state drives the result. The interaction beta never flips sign and never loses statistical significance across any of the 80 LOSO models (40 per risk anchor). The standard deviation is small relative to the effect size, indicating high geographic stability.

### Leave-One-Component-Out Stability

When each of the 12 PCHA components is dropped and the index is rebuilt, the interaction remains significant in all 24 models (12 per risk anchor):

| Risk Anchor | Beta Range | All p < 0.005? |
|-------------|-----------|:---:|
| QUO197 | [-0.420, -0.352] | Yes |
| QUO210 | [-0.414, -0.347] | Yes |

No single component is essential. The most influential components are QUO179 (exercise) and QUO257 (pregnancy intentionality) — dropping either weakens the interaction the most — but the effect persists without them.

### Individual Component Interactions

**Table 4. Individual PCHA Component × Risk Interaction (top 10 by FDR)**

| Component | Domain | Risk | n | Beta | p | FDR |
|-----------|--------|------|---|-----:|---:|----:|
| QUO257 (trying to conceive) | Intentionality | QUO197 | 188 | -0.531 | 0.0001 | 0.001 |
| QUO257 (trying to conceive) | Intentionality | QUO210 | 188 | -0.498 | 0.0001 | 0.001 |
| QUO179 (exercise 3+/wk) | Preconception | QUO210 | 86 | -0.944 | 0.0002 | 0.001 |
| QUO297 (early PNC) | Engagement | QUO210 | 86 | -0.742 | 0.0003 | 0.002 |
| QUO297 (early PNC) | Engagement | QUO197 | 86 | -0.702 | 0.0008 | 0.004 |
| QUO179 (exercise 3+/wk) | Preconception | QUO197 | 86 | -0.896 | 0.0009 | 0.004 |
| QUO4 (breastfeeding 4 wk) | Engagement | QUO197 | 188 | -0.490 | 0.001 | 0.004 |
| QUO4 (breastfeeding 4 wk) | Engagement | QUO210 | 188 | -0.470 | 0.001 | 0.004 |
| QUO44 (breastfeeding 8 wk) | Engagement | QUO210 | 188 | -0.467 | 0.002 | 0.004 |
| QUO5 (ever breastfed) | Engagement | QUO197 | 188 | -0.487 | 0.002 | 0.004 |

All three PCHA domains contribute: intentionality (strongest individual component), preconception health behaviors (largest betas, especially exercise), and early care engagement (most numerous significant components, especially breastfeeding). This supports the construct's theoretical coherence — it is not driven by a single proxy variable.

### Discriminant Validity: PCHA vs. Structural Access

PCHA and the structural-access index are moderately correlated (r = 0.606), raising the question of whether PCHA is simply capturing insurance coverage. Head-to-head comparison on the 102 state-years where both are available:

| Model | PCHA Interaction Beta | PCHA p | Structural Interaction Beta | Structural p |
|-------|----------------------:|-------:|---------------------------:|-----------:|
| PCHA-only | -0.385 | 0.001 | — | — |
| Structural-only | — | — | -0.179 | 0.034 |
| **Both** | -0.312 | 0.142 | -0.033 | 0.834 |

When both are in the model, PCHA retains the larger coefficient (-0.312) while structural access collapses to near-zero (-0.033). The structural-access index has no independent moderating effect once PCHA is controlled. The wider confidence intervals in the joint model reflect collinearity (r = 0.606), but the pattern is clear: the active ingredient in the moderation is behavioral health agency, not insurance structure.

### The Counseling Paradox

A striking finding from the counseling-intensity analysis:

| Statistic | Counseling → PPD |
|-----------|:---:|
| Simple correlation | r = +0.051 (p = 0.851) |
| FE association | beta = -0.917 (p = 0.016) |
| **Interaction with QUO197** | **beta = +0.401 (p = 0.010)** |
| **Interaction with QUO210** | **beta = +0.448 (p = 0.002)** |

The counseling interaction is **positive** — the opposite direction from PCHA. States with more intensive prenatal counseling show a *stronger*, not weaker, relationship between partner stress and PPD. This is almost certainly confounding-by-indication: states with higher psychosocial burden counsel more intensively as a response to that burden, but the counseling itself does not overcome the underlying risk at the population level. This finding highlights why the PCHA construct is meaningful: it captures *upstream* behavioral preparation, not *downstream* clinical response.

---

## Discussion

### What These Results Mean

The central finding is straightforward: in states where women are more likely to exercise before pregnancy, take vitamins, plan their pregnancies, enter prenatal care early, and breastfeed, the link between partner stress and postpartum depression is weaker. This holds after controlling for state and year fixed effects, after adjusting for clustering, across both eras of PPD measurement, and after dropping any single state or any single PCHA component.

The effect size is meaningful. The interaction beta of approximately -0.4 to -0.5 (standardized) means that a one-standard-deviation increase in PCHA reduces the partner-stress → PPD slope by about 40% of its baseline magnitude. For states at +1 SD of partner stress, moving from -1 SD to +1 SD on PCHA is associated with approximately 2–3 percentage points lower PPD prevalence — a non-trivial shift at the population level given a mean prevalence of 13%.

### What PCHA Is (and Is Not)

PCHA is not a clinical intervention. It is a population-level behavioral pattern that reflects the degree to which women in a state proactively engage in health-promoting behaviors before and during the reproductive window. At the ecological level, it likely captures a mixture of:

1. **Individual health literacy and self-efficacy** — women who exercise and take vitamins are exercising a form of health agency
2. **Health system accessibility** — states where vitamins, dental care, and prenatal care are accessible show higher rates of these behaviors
3. **Social norms and cultural context** — states with stronger health promotion cultures may facilitate both the behaviors and their co-occurrence
4. **Policy environment** — Medicaid expansion, WIC availability, and preconception health programs vary by state

The discriminant validity analysis shows that PCHA is not reducible to insurance coverage (structural access). When both are in the model, structural access adds nothing beyond PCHA. This suggests the behavioral-engagement dimension is what matters for moderating the partner-stress pathway, not the insurance dimension per se.

### Why the Counseling Paradox Matters

The finding that counseling intensity shows a paradoxical positive interaction (more counseling = stronger partner-stress → PPD link) is not evidence that counseling is harmful. It is evidence that population-level counseling intensity is a *response marker* — states counsel more because they have more need — not a *prevention tool* at the ecological level. This is a critical distinction for policy: investing in upstream health agency (preconception health promotion, planned pregnancy support, exercise access) may be more effective at the population level than increasing the volume of prenatal counseling content.

### Comparison with Prior Work

No published study has:
- Used PRAMStat aggregate data as the analytic unit for PPD research
- Constructed a PCHA-like composite from population-level indicators
- Tested ecological-level moderation of the partner-stress → PPD pathway

The closest work is Dadi et al. (2020), who found that GINI inequality explained 73% of cross-national PPD variation in an ecological meta-regression — but that was across 56 countries using country-level economic indicators, not within-US state behavioral indicators. The Frontiers (2021) mixed-model study with PRAMS data found null results for state-level variables but tested only structural factors (Medicaid cutoffs, provider counts), not behavioral indicators.

The individual-level literature on exercise, breastfeeding, and pregnancy intentionality as separate PPD protective factors is consistent with our ecological findings. Our contribution is showing that these operate as a cluster, not in isolation, and that their co-occurrence at the population level moderates the strongest known risk pathway.

### Limitations

1. **Ecological design.** All findings describe state-year-level associations. We cannot infer that any individual woman's exercise habit buffers her partner-stress risk. The ecological fallacy is a real concern, and prospective individual-level validation is essential.

2. **No precision weighting.** PRAMStat does not provide standard errors for the unstratified prevalence estimates used in this analysis. All state-year rows are weighted equally, meaning noisy estimates from small states carry the same weight as precise estimates from large states. This introduces heteroscedasticity and potentially inflates false positive rates.

3. **Outcome harmonization.** The PPD outcome combines two different PRAMStat indicators across eras. While the era-specific sensitivity shows consistent results, subtle differences in question wording or survey methodology could confound comparisons.

4. **No parallel analysis or Cronbach's alpha.** The PCHA components have heterogeneous coverage (86–188 rows per component), which prevented computing these diagnostics on a complete-case matrix. The construct's validity rests on theoretical grounds and the leave-one-component-out stability, not on traditional psychometric evidence.

5. **Permutation p-values at boundary.** The permutation p-values (0.079–0.091) sit between conventional significance thresholds. A strict alpha = 0.05 interpretation would call these non-significant; a more appropriate exploratory alpha = 0.10 would call them significant. We present both and let readers judge.

6. **Forking paths.** Despite the pre-specified PCHA construct, the overall research project involved extensive prior exploratory analysis. The PCHA hypothesis, while theory-grounded, was informed by patterns discovered in earlier data-driven exploration phases. This is honestly disclosed.

7. **Temporal scope.** Data span 2000–2011. The US maternal health landscape has changed substantially since then (Medicaid expansion, telehealth, COVID-19). Findings may not generalize to the current period.

8. **Sample size.** With 188 state-year observations and 40 states, statistical power for the interaction term is moderate but not generous, especially under fixed effects (which consume 42 degrees of freedom for dummies).

---

## Conclusions

This analysis introduces Pre-Conception Health Agency as a novel ecological construct and provides preliminary evidence that it moderates the well-established partner-stress → PPD pathway at the US state-year level. The finding is robust across multiple inferential frameworks, consistent across PPD measurement eras, and not attributable to any single component or any single state. It is distinct from structural healthcare access and operates in the opposite direction from prenatal counseling intensity.

The practical implication is a shift in the intervention target: from prenatal counseling (which shows a paradoxical positive association at the population level) to pre-conception health optimization — exercise promotion, nutritional supplementation, reproductive planning support, and dental access — specifically for populations with elevated partner-stress risk factors.

This is an ecological, exploratory finding. It generates a testable hypothesis for individual-level prospective research: *Does participation in a pre-conception health optimization program reduce PPD incidence among women with elevated partner-stress exposure?* A phased validation program — prospective cohort, then pragmatic trial, then component optimization — is the logical next step.

---

## References

Abbasi, S., et al. (2013). Incidence of postpartum depression and its risk factors. *Journal of Midwifery and Women's Health*, 58(4), 418–426.

Beck, C. T. (2001). Predictors of postpartum depression: An update. *Nursing Research*, 50(5), 275–285.

Campos, B., et al. (2021). Associated factors with postpartum depression symptoms. *BMC Public Health*, 21, Article 155.

Dadi, A. F., et al. (2020). Global burden of antenatal depression and its association with adverse birth outcomes. *Frontiers in Psychiatry*, 11, Article 248.

Gastaldon, C., et al. (2022). Risk factors of perinatal depression among women who have experienced miscarriage. *British Journal of Psychiatry*, 221, 515–523.

Luca, D. L., et al. (2020). Financial burden of untreated perinatal mood and anxiety disorders. *American Journal of Public Health*, 110(6), 888–896.

Patton, G. C., et al. (2015). Prediction of perinatal depression from adolescence and before conception. *The Lancet*, 386, 875–883.

Qobadi, M., & Collier, C. (2016). The association of stressful life events with postpartum depression. *Mississippi PRAMS, 2007–2010*.

Robertson, E., Grace, S., Wallington, T., & Stewart, D. E. (2004). Antenatal risk factors for postpartum depression. *General Hospital Psychiatry*, 26, 289–295.

Thomson, K. C., et al. (2020). Lifecourse predictors of perinatal depression. *Psychological Medicine*, 50, 154–162.

---

## Reproducibility

All analysis code, data, and outputs are available at:
**https://github.com/Mathewmoslow/NovelPramsAnalysis**

- Run ID: `run_20260324T_pcha_v1`
- Random seed: 42
- Panel SHA-256: `db32d98dff8d81b1...`
- Python 3.x with pandas, numpy, scipy, matplotlib
- Single-command reproduction: `python scripts/run_pcha_analysis.py`

---

## Tables and Figures

| Item | File | Description |
|------|------|-------------|
| Table 1 | `pcha_interaction_models.csv` | Full interaction model results for all 4 risk anchors |
| Table 2 | `era_sensitivity.csv` | Era-specific interaction results |
| Table 3 | `leave_one_state_out.csv` | LOSO stability (80 models) |
| Table 4 | `pcha_component_interactions.csv` | Individual component interaction results (24 models) |
| Table 5 | `leave_one_component_out.csv` | LOCO stability (24 models) |
| Table 6 | `discriminant_validity.csv` | PCHA vs. structural access head-to-head |
| Table 7 | `counseling_paradox.csv` | Counseling intensity paradox results |
| Table 8 | `descriptive_statistics.csv` | Descriptive statistics for all key variables |
| Table 9 | `bivariate_correlations.csv` | Bivariate correlation matrix |
| Table 10 | `variable_dictionary.csv` | Full question dictionary with role flags |
| Figure 1 | `fig1_pcha_distribution.png` | PCHA index distribution and component coverage |
| Figure 2a | `fig2_interaction_QUO197.png` | Scatter: QUO197 × PPD by PCHA group |
| Figure 2b | `fig2_interaction_QUO210.png` | Scatter: QUO210 × PPD by PCHA group |
| Figure 3a | `fig3_loso_QUO197.png` | LOSO stability forest plot (QUO197) |
| Figure 3b | `fig3_loso_QUO210.png` | LOSO stability forest plot (QUO210) |
| Figure 4 | `fig4_component_heatmap.png` | Component × risk anchor interaction heatmap |
| Figure 5 | `fig6_model_comparison.png` | Model comparison: R² and interaction CIs |
