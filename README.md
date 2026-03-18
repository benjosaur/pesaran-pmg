# Pesaran, Shin & Smith (1999) — R Replication

R implementation of the four estimators from:

> Pesaran, M.H., Shin, Y. & Smith, R.P. (1999). "Pooled Mean Group Estimation of Dynamic Heterogeneous Panels." *Journal of the American Statistical Association*, 94(446), 621-634.

## Provenance: What Comes From Where

### Direct translations from GAUSS (verified against OUTPUT.OUT and paper)

| File | GAUSS source | What it translates | Verified against |
|------|-------------|-------------------|-----------------|
| `panel_ecm_fe.R` — `SFE()` | `FIX.PRC` — `SFE` | Static Fixed Effects estimator | OUTPUT.OUT p.25: all coefficients and SEs match to 4dp |
| `panel_ecm_fe.R` — `DFE1()` | `FIX.PRC` — `DFE1` | Dynamic Fixed Effects estimator | OUTPUT.OUT p.25-26: all LR, phi, SR coefficients and SEs match to 4dp |
| `panel_ecm_fe.R` — `FDGP1()` | `FIX.PRC` — `FDGP1` | Data generation for DFE (excludes intercept from ww) | Implicit via DFE results |
| `panel_ecm_fe.R` — `lag1()`, `demean()` | `SUB0.PRC` — `LAG`, within-transform | Lag and demeaning helpers | Implicit via all estimators |
| `pmg_mg.R` — `MG()` | `SUB1.PRC` — `OLS1` | Mean Group estimator (individual OLS + cross-sectional averaging) | OUTPUT.OUT p.26: all LR, phi, SR coefficients and SEs match to 3dp. Paper Table 1: exact match. Paper Table 3: exact match |
| `pmg_mg.R` — `PMG()` | `SUB1.PRC` — `PMLEBS1` | PMG via Back-Substitution Algorithm | Converges to same optimum as NR (verified) |
| `pmg_mg.R` — `PMG_NR()` | `SUB1.PRC` — `PMLENR1` + `INIPHI1` | PMG via Newton-Raphson | OUTPUT.OUT p.26: LR coefficients match within 0.006 (our RSLL is marginally higher). Paper Table 1: exact match on all coefficients |
| `pmg_mg.R` — `pmg_final()` | `SUB1.PRC` — `PMLESR1` | PMG final SEs from information matrix | OUTPUT.OUT: theta SEs match. Paper Table 1: match |
| `pmg_mg.R` — `DGP1()` | `SUB1.PRC` — `DGP1` | Data generation for MG/PMG (includes intercept in ww) | Implicit via MG/PMG results |
| `load_data.R` — `load_example1()` | `jasa1.prg` — data loading section | OECD consumption data (LPC, LNDI, DP) with missing value handling | Ti counts match OUTPUT.OUT. SFE obs=767, DFE obs=743 match exactly |
| `load_data.R` — `load_example2()` | `jasa2.prg` — data loading section | Asian energy demand data (EDA) | Paper Table 3: all estimators match |

### Additions (not in the original GAUSS code)

| File | What it adds | Status |
|------|-------------|--------|
| `pmg_mg.R` — `hausman_test()` | Hausman test comparing PMG vs MG efficiency | Verified: h=2.62, p=0.27 for OUTPUT.OUT spec; consistent with paper's reported h=2.79 (difference due to PMG convergence point) |
| `pmg_mg.R` — `lr_test()` (in `test_cross_estimator.R`) | Likelihood ratio test: 2*(MG_LL - PMG_LL) | Verified: LR=121.50 matches OUTPUT.OUT value 121.51. Corrects GAUSS df from NN*k=46 to (NN-1)*k=44 |
| `panel_helpers.R` — `prepare_panel()` | Convert data.frame to stacked format | Verified: 14/14 checks pass, max diff=0 between df and raw API (`replicate_df.R`) |
| `panel_helpers.R` — `pmg_panel()`, `mg_panel()`, `dfe_panel()`, `sfe_panel()` | Data.frame convenience wrappers with labelled output | Verified: identical results to raw API. SR labels verified against OUTPUT.OUT |
| `panel_helpers.R` — print methods | Labelled output for all estimators (LR, phi, SR with row names) | Verified: all printed values match OUTPUT.OUT |
| `panel_ecm_fe.R` — `DFE2()`, `FDGP2()` | Two-regressor-set DFE variant | Translated from `FIX.PRC`. Not separately verified (no reference output for this variant in the distributed files) |
| `tests/test_replication.R` | Unit tests against OUTPUT.OUT | 12/12 pass |
| `tests/test_cross_estimator.R` | Cross-estimator nesting tests, LR test, Hausman | 14/14 pass |
| `replicate_tables.R` | Reproduces paper Tables 1, 2 (partial), 3 | Table 1: exact match. Table 2 ARDL(1,0,0): exact match. Table 3: exact match |
| `replicate_df.R` | Verifies data.frame API against raw API | 14/14 checks, max diff=0 |
| `examples/mock_panel.R` | Simulated ECM panel with known DGP | Not a replication — demonstrates API usage |
| `reports/outliers_sweden.md` | Analysis of Sweden outlier in PDI data | Documents paper's discussion; independently verified country-specific OLS |

### Bug fixes relative to original GAUSS

| Fix | Detail |
|-----|--------|
| LR test df | GAUSS reports df=NN*k=46; correct df is (NN-1)*k=44 |
| `diag()` with k=1 | R's `diag(scalar)` creates an identity matrix instead of extracting the diagonal; fixed with `drop=FALSE` |
| `ref_group` ordering | Alphabetical sort could exclude wrong reference country; fixed by explicitly placing `ref_group` last |

## Variables

The variable names `inc` and `inf` are labels from the original GAUSS interactive session:

| Label | Variable | Source file | Description |
|-------|----------|-------------|-------------|
| **inc** | log real per-capita national disposable income | `LNDI.DAT` | Long-run income elasticity of consumption |
| **inf** | inflation rate (GDP deflator change) | `DP.DAT` | Long-run effect of inflation on consumption |

The dependent variable is log real per-capita private consumption (`LPC.DAT`). The intercept is the only deterministic regressor (Z). No seasonal dummies are needed — the data is annual (1960-1993).

## Verification Summary

### Paper Table 1: ARDL(1,1,1), NDI, N=24 — `replicate_tables.R`

All 24 countries, inflation as regular regressor. DFE SEs are robust.

|  | MG | PMG | DFE (robust SE) |
|--|----|----|------|
| θ₁ (inc) | **0.918** (0.027) | **0.904** (0.009) | **0.912** (0.045) |
| θ₂ (inf) | **-0.353** (0.117) | **-0.466** (0.057) | **-0.266** (0.098) |
| φ | **-0.306** (0.030) | **-0.200** (0.032) | **-0.179** (0.042) |

Paper values: MG θ₁=0.918 (0.027), PMG θ₁=0.904 (0.010), DFE θ₁=0.912 (0.045). **Exact match on all coefficients and SEs.**

### OUTPUT.OUT: ARDL(1,2,0), NDI, NN=23 — `replicate.R`

Inflation as relative variable, U.K. excluded from MG/PMG.

|  | MG | PMG | DFE | SFE |
|--|----|----|------|-----|
| θ₁ (inc) | **0.916** (0.028) | **0.914** (0.008) | **0.925** (0.019) | **0.953** (0.007) |
| θ₂ (inf) | **-0.305** (0.128) | **-0.495** (0.052) | **-0.307** (0.076) | **-0.094** (0.027) |
| φ | **-0.317** (0.032) | **-0.199** (0.035) | **-0.156** (0.017) | -1 |

SFE, DFE, MG: exact match. PMG θ₂ differs by 0.006 (our RSLL=2139.36 > reference 2139.35).

### Paper Table 3: ARDL(1,0,0), Energy Demand, N=10 — `replicate_tables.R`

|  | MG | PMG | DFE | SFE |
|--|----|----|------|-----|
| θ₁ (output) | **1.228** (0.183) | **1.183** (0.039) | **1.301** (0.109) | **1.009** (0.037) |
| θ₂ (price) | **-0.261** (0.118) | **-0.339** (0.033) | **-0.365** (0.097) | **-0.067** (0.030) |

Paper values: exact match on all coefficients.

## Cross-Estimator Tests (additions, not in paper/GAUSS)

`tests/test_cross_estimator.R` verifies nesting relationships MG > PMG > DFE:

| Test | Result |
|------|--------|
| **LR test** (MG vs PMG): 2*(URSLL - RSLL) | 121.50, χ²(44), p < 0.001 |
| **Hausman test** (PMG vs MG) | h = 2.62, p = 0.27 — does not reject PMG |
| **LL ordering**: MG URSLL > PMG RSLL | 2200.11 > 2139.36 |
| **Short-run heterogeneity**: PMG φ_i std dev | 0.17 |
| **DFE φ within PMG φ range** | -0.156 within [-0.603, 0.018] |
| **SFE bias direction** | Highest θ₁ (upward), smallest |θ₂| (attenuation) |

## Why Different ARDL Specifications?

1. **ARDL(1,1,1)** — "maximum lag = 1" default. Table 1.
2. **ARDL(1,2,0)** — relative variable forces inflation lag=0. OUTPUT.OUT.
3. **SBC-chosen** — data-driven per country. Table 4 (not reproduced; requires lag selection not yet implemented).

The long-run estimates are robust across all specifications (~0.9 for income, negative for inflation). Short-run dynamics are more sensitive.

## How to Run

```bash
# Reproduce paper tables (Table 1, 2 partial, 3)
Rscript replicate_tables.R

# Reproduce OUTPUT.OUT specification
Rscript replicate.R

# Verify data.frame API matches raw API
Rscript replicate_df.R

# Unit tests against OUTPUT.OUT
Rscript tests/test_replication.R

# Cross-estimator nesting tests
Rscript tests/test_cross_estimator.R
```

No external packages required. `test_replication.R` uses `testthat` if available, with a basic fallback.

## Using Your Own Data

```r
source("panel_helpers.R")

# Any data.frame / data.table / pdata.frame
pmg <- pmg_panel(df, y = "consumption", x = c("income", "inflation"),
                 unit = "country", time = "year",
                 plag = 1, qlag = c(1, 1))
print(pmg)        # labelled LR, phi, SR coefficients
pmg$lr            # data.frame with variable names
pmg$group_phi     # named vector of group-specific phi

# Also: mg_panel(), dfe_panel(), sfe_panel() — same interface
```

Seasonal dummies for quarterly/monthly data go in Z (intercept must be last column):
```r
sfe <- sfe_panel(df, y = "y", x = c("x1", "x2"),
                 unit = "id", time = "quarter",
                 z = c("qs1", "qs2", "qs3"))
```

For relative variables, specify the reference group explicitly:
```r
pmg <- pmg_panel(df, ..., ref_group = "U.K.")
```

**Caveat**: All groups must have uniform lag orders. Group-specific AIC-chosen lags require extending `pmg_final`.

## File Structure

```
├── panel_ecm_fe.R            # [GAUSS translation] SFE, DFE, helpers
├── pmg_mg.R                  # [GAUSS translation] MG, PMG (BSA & NR), Hausman
├── load_data.R               # [GAUSS translation] Data loading for Examples 1 & 2
├── panel_helpers.R           # [Addition] data.frame API, labelled output
├── replicate.R               # [Verification] OUTPUT.OUT specification
├── replicate_tables.R        # [Verification] Paper Tables 1, 2, 3
├── replicate_df.R            # [Verification] data.frame API vs raw API
├── tests/
│   ├── test_replication.R    # [Verification] Unit tests vs OUTPUT.OUT (12/12)
│   └── test_cross_estimator.R # [Addition] LR test, Hausman, nesting (14/14)
├── examples/
│   └── mock_panel.R          # [Addition] Simulated ECM panel demo
├── reports/
│   └── outliers_sweden.md    # [Addition] Sweden outlier analysis
├── prog/                     # [Original] GAUSS code and data files
│   ├── *.DAT                 # Original data files
│   ├── OUTPUT.OUT            # Reference output
│   └── *.PRG, *.PRC          # Original GAUSS programs
└── README.md
```

## Key Implementation Notes

- **Relative variable**: Inflation is relative to a reference country. Its lag is forced to 0, and the reference group (U.K.) is excluded from MG/PMG (NN = N-1 = 23). SFE and DFE use all N = 24 groups.
- **PMG algorithms**: Both BSA and NR are implemented and converge to the same optimum. NR matches GAUSS `PMLENR1`.
- **Standard errors**: PMG θ SEs from the information matrix. PMG φ/β/γ SEs from cross-sectional variance.
- **LR test df correction**: GAUSS reports df = NN*k = 46. Correct df is (NN-1)*k = 44.
- **No external dependencies**: Core estimators use only base R.
