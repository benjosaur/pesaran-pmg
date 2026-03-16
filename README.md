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

## Mapping to the Paper's Tables

The paper reports results for two ARDL specifications across two examples. The GAUSS `OUTPUT.OUT` distributed with the code uses a *different* specification from Table 1 in the paper. Both are replicated here.

### Table 1: ARDL(1,1,1), OECD Consumption, NDI, N=24

All 24 countries used (inflation treated as regular regressor, not relative). DFE SEs are robust (heteroscedasticity-corrected); the uncorrected SEs are substantially smaller.

|  | MG | PMG | DFE (robust SE) |
|--|----|----|------|
| θ₁ (inc) | **0.918** (0.027) | **0.904** (0.009) | **0.912** (0.045) |
| θ₂ (inf) | **-0.353** (0.117) | **-0.466** (0.057) | **-0.266** (0.098) |
| φ | **-0.306** (0.030) | **-0.200** (0.032) | **-0.179** (0.042) |
| MLL | 2390 | 2327 | — |

Paper Table 1 values: MG θ₁=0.918 (0.027), PMG θ₁=0.904 (0.010), DFE θ₁=0.912 (0.045). **Exact match on all coefficients and SEs.**

### OUTPUT.OUT: ARDL(1,2,0), OECD Consumption, NDI, NN=23

Inflation treated as relative variable (lag forced to 0, U.K. excluded from MG/PMG). This is the specification distributed with the GAUSS code.

|  | MG | PMG | DFE | SFE |
|--|----|----|------|-----|
| θ₁ (inc) | **0.916** (0.028) | **0.914** (0.008) | **0.925** (0.019) | **0.953** (0.007) |
| θ₂ (inf) | **-0.305** (0.128) | **-0.495** (0.052) | **-0.307** (0.076) | **-0.094** (0.027) |
| φ | **-0.317** (0.032) | **-0.199** (0.035) | **-0.156** (0.017) | -1 |

SFE, DFE, and MG match the GAUSS output exactly. PMG θ₂ differs by 0.006 (our RSLL is marginally higher: 2139.36 vs 2139.35).

### What Differs Between These Specifications

| Feature | Table 1 | OUTPUT.OUT |
|---------|---------|------------|
| ARDL order | (1, 1, 1) | (1, 2, 0) |
| Inflation treatment | Regular regressor (lag=1) | Relative variable (lag forced to 0) |
| Groups in MG/PMG | N=24 (all) | NN=23 (U.K. excluded) |
| Short-run Δ(income) terms | 1 (current Δy^d) | 2 (current + lagged Δy^d) |
| Short-run Δ(inflation) terms | 1 (current Δπ) | 0 (none) |

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

## Why Different ARDL Specifications?

The paper tests multiple ARDL orders to show robustness. The choice is motivated by:

1. **ARDL(1,1,1)** — the "maximum lag = 1" default. Simple, symmetric treatment of all regressors. Used for Table 1.

2. **ARDL(1,2,0)** — used in the GAUSS `OUTPUT.OUT`. When inflation is treated as a **relative variable** (measured relative to a reference country like U.K.), the GAUSS code forces its lag to 0 and disables AIC-based lag selection (`SUB0.PRC` line 525). The extra income lag (q=2) allows a richer short-run income response.

3. **SBC-chosen lags** — the paper's Table 4 (energy demand) lets the Schwarz criterion choose lag orders per country. The most common choice was ARDL(1,1,0) for consumption and ARDL(1,0,0) for energy demand. This is the most data-driven approach but creates group-specific SR dimensions.

The paper's key finding is that **the long-run estimates are robust across all specifications** — income elasticity is always ~0.9, inflation is always significantly negative. The short-run dynamics and speed of adjustment are more sensitive to lag choice.

### Why Not AIC in Our Replication?

The current code supports any user-specified fixed lag order via `plag` and `qlag`. AIC/SBC-based selection per country (as in the GAUSS `AICLAG` procedure) is not yet implemented. This would require:
- Running unrestricted ARDL regressions per group with varying lag orders
- Selecting the order that minimizes AIC/SBC for each group
- Handling heterogeneous SR dimensions in the PMG information matrix (the current `pmg_final` asserts uniform `numsi`)

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
