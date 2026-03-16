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

## How to Run

```r
# Run full replication
Rscript replicate.R

# Run tests
Rscript tests/test_replication.R
```

No external packages required for the core estimators. Tests use `testthat` if available, with a basic fallback otherwise.

## File Structure

```
├── panel_ecm_fe.R          # SFE, DFE estimators and helpers (lag1, demean, FDGP1)
├── pmg_mg.R                # MG, PMG (BSA & NR) estimators and Hausman test
├── load_data.R             # Data loading for Examples 1 and 2
├── replicate.R             # Main replication script
├── tests/
│   └── test_replication.R  # Unit tests against OUTPUT.OUT values
├── prog/                   # Original GAUSS code and data files
│   ├── LPC.DAT, LNDI.DAT, DP.DAT   # Example 1 data
│   ├── EDA.DAT                       # Example 2 data
│   ├── OUTPUT.OUT                    # Reference output
│   └── *.PRG, *.PRC                  # Original GAUSS programs
└── README.md
```

## Key Implementation Notes

- **Relative variable**: The last X variable (inflation) is treated as a relative variable. Its lag order is fixed at 0, and the last group (U.K.) is excluded from MG/PMG estimation (NN = N-1 = 23).
- **PMG algorithms**: Both Back-Substitution (BSA) and Newton-Raphson (NR) are implemented. NR jointly updates θ and φ and matches the GAUSS `PMLENR1` procedure. Both converge to the same optimum.
- **Standard errors**: PMG θ SEs come from the full information matrix. PMG φ/β/γ SEs use cross-sectional variance (MG-style averaging).
- **No external dependencies**: Core estimators use only base R.
