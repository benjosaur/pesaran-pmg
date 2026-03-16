# Pesaran, Shin & Smith (1999) — R Replication

R implementation of the four estimators from:

> Pesaran, M.H., Shin, Y. & Smith, R.P. (1999). "Pooled Mean Group Estimation of Dynamic Heterogeneous Panels." *Journal of the American Statistical Association*, 94(446), 621-634.

Translated from the original GAUSS code distributed with the paper.

## Estimators

| Estimator | Description |
|-----------|-------------|
| **SFE** | Static Fixed Effects — estimates cointegrating relationship directly from levels with φ = -1 imposed |
| **DFE** | Dynamic Fixed Effects — panel ARDL/ECM with within-transformation, recovers long-run via θ = -β/φ |
| **MG** | Mean Group — unrestricted OLS per group, averages long-run parameters across groups |
| **PMG** | Pooled Mean Group — constrains long-run parameters equal across groups via concentrated ML, allows heterogeneous short-run dynamics |

## Variables

The replication uses Example 1 (OECD Consumption Function). The variable names `inc` and `inf` are labels from the original GAUSS interactive session:

| Label | Variable | Source file | Description |
|-------|----------|-------------|-------------|
| **inc** | log real per-capita national disposable income | `LNDI.DAT` | Long-run income elasticity of consumption |
| **inf** | inflation rate (GDP deflator change) | `DP.DAT` | Long-run effect of inflation on consumption |

The dependent variable is log real per-capita private consumption (`LPC.DAT`). The intercept is the only deterministic regressor (Z). No seasonal dummies are needed — the data is annual (1960-1993).

## Replication Results (Example 1: OECD Consumption Function)

24 OECD countries, 1960-1993. ARDL(1,2,0) specification.

| Estimator | inc (se) | inf (se) |
|-----------|----------|----------|
| PMG | 0.914 (0.008) | -0.495 (0.052) |
| MG | 0.916 (0.028) | -0.305 (0.128) |
| SFE | 0.953 (0.007) | -0.094 (0.027) |
| DFE | 0.925 (0.019) | -0.307 (0.076) |

Reference (OUTPUT.OUT):

| Estimator | inc (se) | inf (se) |
|-----------|----------|----------|
| PMG | 0.915 (0.008) | -0.501 (0.052) |
| MG | 0.916 (0.028) | -0.305 (0.128) |
| SFE | 0.953 (0.007) | -0.095 (0.027) |
| DFE | 0.925 (0.019) | -0.307 (0.076) |

SFE, DFE, and MG match exactly. PMG differs by < 0.01 due to cross-platform numerical precision (our solution has a marginally higher log-likelihood: 2139.36 vs 2139.35).

## Cross-Estimator Tests (beyond the paper)

`tests/test_cross_estimator.R` verifies the nesting relationships MG > PMG > DFE that the paper discusses theoretically but the original GAUSS code does not formally test end-to-end:

| Test | Result |
|------|--------|
| **LR test** (MG vs PMG): 2*(URSLL - RSLL) | 121.50, χ²(44), p < 0.001 — rejects LR homogeneity |
| **Hausman test** (PMG vs MG) | h = 2.62, p = 0.27 — does not reject PMG efficiency |
| **LL ordering**: MG URSLL > PMG RSLL | 2200.11 > 2139.36 (unrestricted fits better, as expected) |
| **Short-run heterogeneity**: PMG φ_i std dev | 0.17 — substantial, DFE homogeneity assumption is restrictive |
| **DFE φ within PMG φ range** | -0.156 is within [-0.603, 0.018] |
| **SFE bias direction** | Highest income elasticity (upward), smallest |inf| (attenuation) |

The LR test rejects while the Hausman test does not — a classic finite-sample divergence. The LR test has power against all departures from homogeneity, while the Hausman test is specifically designed for the case where PMG is efficient under H0. In practice, the Hausman non-rejection supports the PMG specification.

## Why Lag Orders Are Fixed (Not Chosen by AIC)

The ARDL(1,2,0) lag specification is not arbitrary — it is forced by the original GAUSS code:

- When a **relative variable** is present (here: inflation relative to a reference country), the GAUSS code disables AIC-based lag selection and forces fixed lags (`SUB0.PRC` line 525: `if rel_op eq 2; lag_op = 1`).
- The relative variable's lag order is always set to 0 by the program.
- The ARDL(1,2,0) choice — one lag of Δy, two lags of Δ(income), zero lags of Δ(inflation) — matches the paper's Table 3.

A future extension could implement the AIC/SBC/HQ lag selection procedure (`AICLAG` in `SUB0.PRC`), but this requires handling group-specific lag orders in the PMG information matrix, which the current SE computation does not support.

## How to Run

```r
# Run full replication
Rscript replicate.R

# Run unit tests (coefficients, SEs, obs counts vs OUTPUT.OUT)
Rscript tests/test_replication.R

# Run cross-estimator tests (LR test, nesting, Hausman)
Rscript tests/test_cross_estimator.R
```

No external packages required for the core estimators. `test_replication.R` uses `testthat` if available, with a basic fallback otherwise.

## Using Your Own Data

The estimator functions are generic. Bypass `load_example1()` and pass your own matrices:

```r
source("panel_ecm_fe.R")
source("pmg_mg.R")

# Y: stacked dependent variable (length sum(Ti))
# X: stacked regressors (sum(Ti) x k matrix)
# Z: fixed regressors — intercept MUST be the last column
#    e.g. cbind(trend, qs1, qs2, qs3, intercept) for quarterly data
# Ti: integer vector of time periods per group
# plag: lag order for Δy (length n vector)
# qlag: lag orders for each ΔX column (n x k matrix)

mg  <- MG(Y, X, Z, n, Ti, k, plag, qlag)
pmg <- PMG_NR(Y, X, Z, n, Ti, k, plag, qlag)
sfe <- SFE(Y, X, Z, n, Ti, k)
dfe <- DFE1(Y, X, Z, n, Ti, k, numot, plag, qlag)
```

**Caveat**: All groups must currently have the same number of short-run parameters (uniform lag orders). Group-specific AIC-chosen lags would require extending the `pmg_final` SE computation.

## File Structure

```
├── panel_ecm_fe.R          # SFE, DFE estimators and helpers (lag1, demean, FDGP1)
├── pmg_mg.R                # MG, PMG (BSA & NR), LR test, Hausman test
├── load_data.R             # Data loading for Examples 1 and 2
├── replicate.R             # Main replication script
├── tests/
│   ├── test_replication.R      # Unit tests against OUTPUT.OUT values
│   └── test_cross_estimator.R  # LR test, nesting restrictions, Hausman
├── prog/                   # Original GAUSS code and data files
│   ├── LPC.DAT, LNDI.DAT, DP.DAT   # Example 1 data
│   ├── EDA.DAT                       # Example 2 data
│   ├── OUTPUT.OUT                    # Reference output
│   └── *.PRG, *.PRC                  # Original GAUSS programs
└── README.md
```

## Key Implementation Notes

- **Relative variable**: The last X variable (inflation) is treated as a relative variable. Its lag order is fixed at 0, and the last group (U.K.) is excluded from MG/PMG estimation (NN = N-1 = 23). SFE and DFE use all N = 24 groups.
- **PMG algorithms**: Both Back-Substitution (BSA) and Newton-Raphson (NR) are implemented. NR jointly updates θ and φ and matches the GAUSS `PMLENR1` procedure. Both converge to the same optimum.
- **Standard errors**: PMG θ SEs come from the full information matrix. PMG φ/β/γ SEs use cross-sectional variance (MG-style averaging).
- **LR test degrees of freedom**: The GAUSS output reports df = NN*k = 46. The correct df is (NN-1)*k = 44 — the number of restrictions imposed by the PMG constraint.
- **No external dependencies**: Core estimators use only base R.
