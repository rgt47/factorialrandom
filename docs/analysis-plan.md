# Analysis Plan: Matched-Tuple Randomization for Factorial Designs

## 1. Problem Statement

Factorial clinical trials assign subjects to combinations of two or
more treatment factors. A 2x2 design has four arms; a 2x2x2 has
eight. Simple randomization in these settings can produce substantial
covariate imbalance across arms, particularly when sample sizes are
moderate or baseline risk profiles are heterogeneous.

When a trial-ready cohort (TRC) is available -- a large pool of
pre-screened, characterized subjects -- the investigator has the
opportunity to exploit baseline covariate information at the design
stage. Specifically, subjects can be grouped into matched tuples of
size k (where k equals the number of factorial treatment
combinations), and treatment assignment can be constrained so that
each member of a tuple receives a distinct treatment combination.

This project develops the statistical theory, algorithms, and
software for matched-tuple randomization in factorial designs, with
a motivating application to Alzheimer's disease prevention trials
that draw from trial-ready cohorts.

## 2. Specific Aims

### Aim 1: Matching Algorithm

Develop and evaluate algorithms for forming matched k-tuples from a
pool of N subjects based on multivariate baseline covariates. Key
design choices include:

- Distance metric (Mahalanobis, propensity score, prognostic score)
- Matching strategy (optimal vs. greedy, with and without calipers)
- Handling of residual unmatched subjects
- Computational scalability for k > 2

### Aim 2: Randomization and Inference

Define the randomization distribution for matched-tuple factorial
designs and derive:

- Exact and asymptotic distributions of treatment effect estimators
- Variance estimators that account for the matched structure
- Randomization-based p-values and confidence intervals
- Connections to rerandomization (Morgan and Rubin, 2012)

### Aim 3: Efficiency Analysis

Quantify the precision gain of matched-tuple randomization relative
to simple randomization and stratified randomization, as a function
of:

- Within-tuple covariate homogeneity
- Number of factors and levels
- Outcome-covariate correlation (prognostic strength)

### Aim 4: Simulation Study

Evaluate operating characteristics under realistic Alzheimer's
disease trial scenarios, comparing matched-tuple factorial designs
against simple randomization, stratified randomization, minimization,
and rerandomization.

### Aim 5: Software Implementation

Deliver an R package (`factorialrandom`) that implements:

- Matched-tuple formation from baseline covariate data
- Constrained randomization within tuples
- Treatment effect estimation with valid inference
- Diagnostic and visualization tools

## 3. Background and Rationale

### 3.1 Why Matching Before Randomization

Randomization guarantees unbiased estimation of treatment effects on
average, but any single randomization may produce poor covariate
balance. In small to moderate sample sizes (50 to 500 per arm),
chance imbalance is common and inflates variance. Blocking,
stratification, minimization, and rerandomization all address this
problem by restricting the set of acceptable allocations.

Matched-pair designs represent the most aggressive form of
restriction: by pairing subjects on covariates and randomizing
within pairs, all between-pair variation is eliminated from the
treatment comparison. Bai (2022) showed that among all stratified
designs with 1:1 allocation, a specific matched-pair design achieves
maximum precision for the average treatment effect.

Extending matching from pairs (k = 2) to tuples (k > 2) is
necessary when the number of treatment arms exceeds two, as in
factorial designs.

### 3.2 Factorial Designs

Factorial experiments simultaneously evaluate multiple treatment
factors and their interactions (Fisher, 1935). They are efficient
because each observation contributes information about every factor.
In clinical trials, factorial designs are used to evaluate drug
combinations, dose-schedule interactions, and multi-component
interventions (Collins et al., 2007).

The challenge is that factorial designs multiply the number of arms.
A 2x2 design needs k = 4; a 2x3 needs k = 6; a 2x2x2 needs k = 8.
Covariate balance must be maintained not only across main effect
contrasts but also across interaction contrasts. Matched tuples of
size k naturally ensure perfect balance for all estimable contrasts.

### 3.3 Trial-Ready Cohorts in Alzheimer's Disease

The trial-ready cohort paradigm, developed for Alzheimer's disease
prevention trials, pre-screens and characterizes large pools of
potential participants. The TRC-PAD project (Aisen et al., 2020)
enrolled over 30,000 individuals in its web-based screening
component. EPAD (Solomon et al., 2022) recruited over 2,000
participants across Europe. DIAN-TU (Bateman et al., 2012)
maintains a registry of individuals carrying autosomal dominant AD
mutations.

These cohorts are ideal settings for matched-tuple randomization
because:

- Large pools (N >> k * number of tuples) permit high-quality
  matching.
- Rich baseline characterization (amyloid PET, tau PET, cognitive
  composites, APOE genotype, demographics) provides strong
  prognostic covariates.
- Prevention trials have long durations (3--4 years) and high per-
  subject costs, making efficiency gains from matching particularly
  valuable.
- Combination therapies (e.g., anti-amyloid + anti-tau) create
  natural factorial structures.

## 4. Proposed Methods

### 4.1 Matching Algorithm

#### 4.1.1 Distance Metric

For each pair of subjects (i, j), compute a distance d(i, j) based
on their baseline covariate vectors X_i and X_j. Candidates:

- **Mahalanobis distance:** d(i,j) = sqrt((X_i - X_j)' S^{-1}
  (X_i - X_j)), where S is the sample covariance matrix of X. This
  is affinely invariant and appropriate when covariates are
  continuous.

- **Prognostic score:** Fit a model predicting the outcome from
  covariates (using historical data or the screening cohort), then
  match on the estimated prognostic score. This focuses matching on
  the dimension most relevant to precision.

- **Rank-based Mahalanobis distance:** Replace each covariate with
  its rank before computing Mahalanobis distance. Robust to outliers
  and nonlinear covariate effects.

#### 4.1.2 Tuple Formation

Given N subjects and k arms, form floor(N/k) tuples by solving:

  minimize sum_{tuples} (within-tuple dispersion)

where within-tuple dispersion can be defined as the sum of pairwise
distances, the maximum pairwise distance (minimax), or the trace of
the within-tuple covariance matrix.

For k = 2, this reduces to optimal pair matching, solvable in O(N^3)
by the Hungarian algorithm (or faster approximations). For k > 2,
the problem is NP-hard in general, but good greedy and local-search
heuristics exist:

- **Sequential greedy:** Select the first subject, find the k-1
  nearest neighbors to form a tuple, remove all k from the pool,
  repeat.
- **Optimal transport:** Formulate as a partitioning problem and
  solve via integer programming (feasible for N < 1000).
- **Hierarchical clustering:** Perform hierarchical clustering and
  cut at the level producing clusters of size k.

#### 4.1.3 Handling Remainders

When N is not divisible by k, the remaining N mod k subjects can be:

- Excluded (wasteful but clean)
- Assigned to a smaller incomplete tuple with restricted
  randomization
- Assigned by simple randomization (hybrid design)

### 4.2 Randomization Within Tuples

Within each matched tuple of size k, the k treatment combinations
are assigned uniformly at random to the k subjects. This produces
k! equally likely allocations per tuple, and the total randomization
space has (k!)^M elements, where M is the number of tuples.

This constrained randomization guarantees:

- Each treatment combination is assigned to exactly one subject per
  tuple.
- Perfect covariate balance across treatment arms at the tuple level.
- Balance for all main-effect and interaction contrasts.

### 4.3 Estimation and Inference

#### 4.3.1 Point Estimation

The difference-in-means estimator, computed as the average outcome
in one arm minus the average in another, is unbiased under the
matched-tuple randomization distribution.

For factorial effects, define contrasts using the standard 2^K
notation. Each main effect and interaction is a signed linear
combination of the cell means.

#### 4.3.2 Variance Estimation

Under matched-tuple randomization, the variance of the treatment
effect estimator differs from the variance under simple
randomization. Two approaches:

- **Neyman-type conservative variance:** Generalize the Neyman
  (1923) variance formula to the matched-tuple setting, accounting
  for within-tuple correlation of potential outcomes.
- **Randomization-based variance:** Compute the exact variance
  under the constrained randomization distribution. This requires
  enumerating or sampling from (k!)^M allocations.

#### 4.3.3 Regression Adjustment

Following Lin (2013) and Fogarty (2018), supplement the
difference-in-means with regression adjustment on baseline
covariates. The matched design and regression adjustment provide
complementary precision gains: matching removes between-tuple
variation, and regression adjustment captures residual within-tuple
covariate-outcome association.

### 4.4 Simulation Study Design

#### 4.4.1 Data Generating Mechanisms

Simulate data mimicking Alzheimer's disease prevention trials:

| Factor                  | Levels                       |
|-------------------------|------------------------------|
| Sample size (total)     | 100, 200, 500, 1000         |
| Design                  | 2x2, 2x3, 2x2x2            |
| Covariate dimension     | 3, 8, 20                    |
| Prognostic strength     | Weak (R^2=0.1), moderate    |
|                         |   (R^2=0.3), strong (R^2=0.5)|
| Treatment effect        | Null, small, moderate        |
| Pool/sample ratio       | 1.0, 1.5, 2.0, 3.0          |

The pool/sample ratio controls how much excess the TRC provides
for matching. A ratio of 1.0 means all subjects are used; a ratio
of 3.0 means 3N subjects are available and the best N are selected
into tuples.

#### 4.4.2 Comparator Methods

1. **Simple randomization** -- no covariate balancing.
2. **Stratified randomization** -- stratify on 1--2 key covariates.
3. **Minimization** -- Pocock-Simon dynamic balancing.
4. **Rerandomization** -- Morgan-Rubin Mahalanobis criterion.
5. **Matched-tuple randomization** (proposed).

#### 4.4.3 Evaluation Metrics

- Covariate balance (standardized mean differences)
- Bias of treatment effect estimators
- Variance and MSE of treatment effect estimators
- 95% CI coverage probability
- Power under non-null treatment effects
- Type I error rate under the null

### 4.5 Software Design

```r
# Form matched tuples
tuples <- match_tuples(
  data,
  k = 4,
  covariates = c('age', 'amyloid', 'apoe4', 'pacc'),
  distance = 'mahalanobis',
  method = 'optimal'
)

# Randomize within tuples
allocation <- randomize_tuples(tuples)

# Estimate factorial effects
fit <- estimate_effects(
  outcome ~ A * B,
  data = trial_data,
  tuples = tuples,
  method = 'neyman'
)

# Randomization inference
pval <- ri_test(fit, n_perm = 10000)
```

## 5. Timeline

| Phase | Duration | Deliverable |
|-------|----------|-------------|
| 1. Literature review   | 2 weeks | Final review doc |
| 2. Matching algorithms | 3 weeks | match_tuples()   |
| 3. Randomization       | 2 weeks | randomize_tuples() |
| 4. Inference           | 3 weeks | estimate_effects(), ri_test() |
| 5. Simulation          | 3 weeks | Full simulation results |
| 6. Package polish      | 2 weeks | Docs, vignettes, tests |
| 7. Manuscript          | 4 weeks | Draft for Biometrics |

Total: approximately 19 weeks.

## 6. Key References

- Bai (2022). Optimality of matched-pair designs. AER.
- Morgan and Rubin (2012). Rerandomization. Annals of Statistics.
- Imai (2008). Variance identification, matched-pair. Stat Med.
- Dasgupta, Pillai, Rubin (2015). Causal inference from 2^K
  factorials. JRSSB.
- Aisen et al. (2020). TRC-PAD overview. JPAD.
- Greevy et al. (2004). Optimal multivariate matching. Biostatistics.
- Fisher (1935). The Design of Experiments.
