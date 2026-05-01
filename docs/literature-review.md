# Literature Review: Matched-Tuple Randomization for Factorial
# Designs with Application to Trial-Ready Cohorts

## 1. Introduction

This review synthesizes three bodies of literature that converge on
the design problem addressed in this project: (1) matched and
covariate-adaptive randomization methods, (2) factorial experimental
designs in clinical trials, and (3) trial-ready cohorts (TRCs) in
Alzheimer's disease research. The goal is to identify the
methodological gap -- the absence of matched-tuple randomization
procedures specifically developed for factorial treatment structures
-- and to motivate the proposed research.

## 2. Matched and Covariate-Adaptive Randomization

### 2.1 Foundations

The tension between randomization (for validity) and balance (for
precision) has been recognized since the earliest experimental
designs. Fisher (1935) introduced randomization as the basis for
inference but also advocated blocking to reduce variance.
Cox (1958) formalized the analysis of blocked designs and
established the principle that 'block what you can, randomize what
you cannot.'

Simple randomization assigns treatments with fixed probabilities
regardless of covariates. While unbiased, it can produce substantial
covariate imbalance in moderate sample sizes. Efron (1971) proposed
the biased coin design to force sequential balance. Pocock and
Simon (1975) introduced minimization, a dynamic balancing procedure
that assigns each incoming subject to the arm that minimizes a
global imbalance function across multiple stratification factors.
Minimization remains widely used in clinical trials but its
inferential properties require care (Rosenberger and Lachin, 2016).

### 2.2 Matched-Pair Randomization

In matched-pair designs, subjects are grouped into pairs based on
baseline covariates, and within each pair, one subject is randomized
to treatment and the other to control. This eliminates all
between-pair variability from the treatment comparison.

Imai (2008) provided a rigorous randomization-based analysis of
matched-pair designs, deriving the variance of the
difference-in-means estimator without modeling assumptions. The key
result is that the variance under matched-pair randomization depends
on within-pair differences in potential outcomes, not on
between-pair differences. When pairs are well matched, this variance
is substantially smaller than under simple randomization.

Bai (2022) proved a notable optimality result: among all stratified
randomization schemes that treat each unit with probability one-half,
a specific matched-pair design achieves the minimum asymptotic
variance for the average treatment effect estimator. In simulations
based on ten published RCTs, this optimal design reduced standard
errors by 10% on average and up to 34% relative to the original
designs.

Fogarty (2018) showed that regression adjustment within matched-pair
designs provides further precision gains, and that the two approaches
(matching and regression) are complementary rather than redundant.

### 2.3 Rerandomization

Morgan and Rubin (2012) formalized rerandomization: repeatedly
generate random allocations and accept only those satisfying a
covariate balance criterion (e.g., the Mahalanobis distance between
treatment group covariate means falls below a threshold).
Rerandomization occupies a middle ground between simple
randomization and deterministic designs: it improves balance while
preserving the stochastic basis for inference.

Li, Ding, and Rubin (2018) developed the asymptotic theory of
rerandomization, showing that the distribution of the treatment
effect estimator under rerandomization is a truncated normal.
Morgan and Rubin (2015) extended rerandomization to balance tiers of
covariates of varying importance.

### 2.4 Optimal Matching and Algorithmic Approaches

Greevy et al. (2004) developed methods for optimal multivariate
matching before randomization, formulating the problem as a minimum-
cost network flow. Their approach produces globally optimal pair
matches and extends to more complex structures via integer
programming.

Balzer et al. (2015) introduced adaptive pair-matching for
community-randomized trials, combining matching with targeted
maximum likelihood estimation (TMLE) for effect estimation. Their
approach adaptively selects matching variables to maximize precision.

Bertsimas, Johnson, and Kallus (2015) demonstrated that optimization-
based designs can dramatically outperform randomization in small
samples, providing a theoretical framework for when deterministic
balance should be preferred over random balance.

Kallus (2018) characterized the optimal a priori balance criterion,
showing that the Mahalanobis distance used in rerandomization is
optimal under a specific semiparametric efficiency criterion.

### 2.5 Covariate-Adaptive Randomization Inference

Bai, Tabord-Meehan, and Romano (2024) recently addressed inference
under covariate-adaptive randomization in matched designs. Their
framework modifies permutation probabilities to account for the
restricted randomization, avoiding the need to model outcomes or
exclude matched sets.

### 2.6 Extension to k > 2

The literature on matched-pair randomization (k = 2) is
well-developed. Extension to matched tuples (k > 2) is considerably
less studied. Martin et al. (1993) examined matching in community
intervention studies and noted that with small numbers of units,
the loss of degrees of freedom from matching can outweigh precision
gains -- a concern that is amplified when k is large. This
tradeoff is central to the proposed research.

The algorithmic challenge also increases: optimal matching into
tuples of size k > 2 is a partitioning problem that is NP-hard in
general. Practical approaches require greedy heuristics, local
search, or integer programming formulations.

## 3. Factorial Designs in Clinical Trials

### 3.1 Classical Theory

Fisher (1935) introduced factorial designs to study multiple factors
simultaneously. The fundamental efficiency of factorial designs lies
in the hidden replication principle: each observation contributes to
the estimation of every main effect and interaction. A 2^K factorial
provides the same precision for each main effect as K separate
one-factor experiments of the same total size.

Montgomery (2017) and Wu and Hamada (2009) provide thorough
modern treatments of factorial design theory, including fractional
factorials and response surface methods.

### 3.2 Causal Inference Framework

Dasgupta, Pillai, and Rubin (2015) recast 2^K factorial designs
within the potential outcomes framework. They defined factorial
effects as contrasts of potential outcomes and derived exact
randomization-based inference for main effects and interactions.
Their framework reveals that interaction effects require stronger
assumptions for identification than main effects.

Lu (2016) extended this framework to incorporate covariate
adjustment, showing that regression adjustment is valid and
precision-enhancing for factorial effects estimated under
randomization-based inference.

### 3.3 Multi-Component Interventions

Collins et al. (2007) introduced the Multiphase Optimization Strategy
(MOST), which uses factorial experiments to identify the active
components of complex interventions before assembling the optimized
package. This approach has been widely adopted in behavioral science
and public health, where interventions often combine multiple
components (e.g., counseling modality, reminder frequency, incentive
level).

In Alzheimer's disease, factorial structures arise naturally when
evaluating combination therapies (e.g., anti-amyloid + anti-tau +
anti-inflammatory) or when optimizing intervention parameters within
a prevention program.

### 3.4 Gap: Matching in Factorial Designs

The intersection of matched randomization and factorial design is
essentially unexplored. Existing matched designs handle k = 2;
existing factorial theory assumes simple or stratified randomization.
No published method provides:

- Algorithms for forming matched k-tuples where k is the number of
  factorial treatment combinations.
- Exact or asymptotic variance formulas for factorial effects under
  matched-tuple randomization.
- Randomization-based inference procedures that exploit the tuple
  structure.
- Software implementing these methods.

## 4. Trial-Ready Cohorts in Alzheimer's Disease

### 4.1 The Trial-Ready Cohort Paradigm

The trial-ready cohort concept emerged from the recognition that
Alzheimer's disease prevention trials face severe recruitment
challenges. Screening failure rates exceed 80% in preclinical AD
trials, primarily because amyloid positivity (required for
inclusion) is only identifiable through PET imaging or lumbar
puncture. By pre-screening and characterizing large pools of
potential participants, TRCs reduce the time and cost of enrollment.

### 4.2 TRC-PAD

The Trial-Ready Cohort for Preclinical/Prodromal Alzheimer's Disease
(TRC-PAD) was established as a collaborative effort across the
Alzheimer's Clinical Trials Consortium (ACTC). Aisen et al. (2020)
described the project architecture:

1. **APT Webstudy:** An online longitudinal study that screens
   participants quarterly with cognitive assessments and subjective
   complaint measures. By 2020, over 30,000 individuals had
   consented.

2. **Risk prediction:** Machine learning models predict brain amyloid
   elevation from demographic, genetic (APOE), and cognitive data.
   High-risk individuals are invited for in-person screening.

3. **In-person screening:** Cognitive battery, APOE genotyping, and
   plasma biomarker assessment (p-tau217, A-beta 42/40 ratio).
   Recent advances in plasma biomarkers (McDonough et al., 2024)
   have dramatically improved prescreening efficiency.

4. **TRC enrollment:** Biomarker-confirmed eligible participants
   enter the trial-ready cohort and are followed longitudinally
   until a matching clinical trial opens.

The TRC-PAD design targets enrollment of approximately 2,000
participants: 1,000 preclinical and 1,000 prodromal AD.

### 4.3 AHEAD 3-45

The AHEAD 3-45 study (Rafii et al., 2023) illustrates how trial-
ready cohorts feed directly into clinical trials. This study
evaluated lecanemab in preclinical AD, using a sister-trial design:

- **A3 trial:** Phase 2, intermediate amyloid (20--40 Centiloids),
  PET imaging endpoints.
- **A45 trial:** Phase 3, elevated amyloid (>40 Centiloids),
  cognitive composite primary endpoint.

Both trials share a common screening process, drawing from the
TRC-PAD infrastructure and direct site recruitment. The screening
includes a blood-based biomarker step to 'screen out' individuals
unlikely to have elevated amyloid, followed by PET confirmation.

### 4.4 DIAN-TU

The Dominantly Inherited Alzheimer Network (Bateman et al., 2012)
maintains a registry of individuals carrying autosomal dominant AD
mutations (PSEN1, PSEN2, APP). These individuals have near-certain
disease onset, typically in their 30s--50s, making them an ideal
population for prevention trials. The DIAN-TU platform
(Bateman et al., 2017) implemented adaptive designs with a shared
placebo arm and cognitive run-in period.

### 4.5 EPAD

The European Prevention of Alzheimer's Dementia (EPAD) programme
(Solomon et al., 2022) recruited over 2,000 research participants
into a longitudinal cohort study across 39 partner institutions.
EPAD's design explicitly envisioned an adaptive proof-of-concept
trial platform drawing from the longitudinal cohort, though IMI
funding ended before these trials launched.

### 4.6 Relevance to Matched-Tuple Factorial Designs

Trial-ready cohorts create the ideal conditions for matched-tuple
randomization:

1. **Large characterized pools.** TRC-PAD's 30,000+ web
   participants and 2,000 in-person participants provide ample
   excess for high-quality matching.

2. **Rich covariate profiles.** Amyloid PET, tau PET, APOE genotype,
   cognitive composites (PACC; Donohue et al., 2014), demographics,
   and comorbidities provide strong prognostic covariates for
   matching.

3. **Factorial treatment opportunities.** Combination therapy trials
   (anti-amyloid + anti-tau, or drug + lifestyle intervention) create
   natural factorial structures. The next generation of AD
   prevention trials is expected to test combinations
   (Cummings et al., 2024).

4. **High per-subject cost.** Prevention trials run 3--4 years per
   subject with PET imaging and cognitive testing. Even modest
   efficiency gains from matching translate to substantial cost
   savings or power improvements.

5. **Biomarker-defined subgroups.** The A3/A45 distinction in AHEAD
   illustrates that subjects are already stratified by biomarker
   level. Within each stratum, further matching on additional
   covariates is feasible and beneficial.

## 5. Summary and Gap Analysis

The three literatures reviewed -- matched randomization, factorial
design, and trial-ready cohorts -- converge on a clear
methodological need:

| Literature             | Contribution          | Limitation            |
|------------------------|-----------------------|-----------------------|
| Matched randomization  | Optimal k=2 designs   | No k>2 theory         |
| Rerandomization        | Balance criteria      | Not tuple-structured  |
| Factorial design       | Efficient multi-factor| Assumes simple rand.  |
| Causal inference       | Potential outcomes    | No matched factorial  |
| Trial-ready cohorts    | Large characterized   | Under-exploited for   |
|                        | pools                 | design optimization   |

No existing method combines matched-tuple formation (k > 2) with
factorial treatment structure and provides valid inference under the
resulting constrained randomization distribution. This project
fills that gap.

## 6. References

See `analysis/report/references.bib` for the complete bibliography.
Key references by topic:

**Matched randomization:** Bai (2022), Imai (2008), Morgan and
Rubin (2012, 2015), Greevy et al. (2004), Balzer et al. (2015),
Bai et al. (2024).

**Factorial design:** Fisher (1935), Dasgupta et al. (2015),
Lu (2016), Collins et al. (2007), Montgomery (2017).

**Trial-ready cohorts:** Aisen et al. (2020), Rafii et al. (2023),
Bateman et al. (2012, 2017), Solomon et al. (2022).

**General:** Rosenberger and Lachin (2016), Rubin (2008),
Lin (2013), Kernan et al. (1999).
