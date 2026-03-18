# Outliers in Heterogeneous Panel Estimation: The Sweden Case

## The Problem in One Table

Country-specific OLS regressions for the OECD consumption function using private disposable labor income (PDI), ARDL(1,1,1):

| Country | phi | theta_inc | theta_inf | sigma | T |
|---------|----:|----------:|----------:|------:|--:|
| Australia | -0.218 | 1.032 | -0.712 | 0.007 | 31 |
| Belgium | -0.293 | 0.805 | -1.188 | 0.009 | 22 |
| Canada | -0.337 | 0.951 | -1.120 | 0.011 | 32 |
| Finland | -0.345 | 0.949 | -0.258 | 0.015 | 32 |
| France | -0.075 | 0.231 | -1.601 | 0.007 | 23 |
| Italy | -0.443 | 0.814 | -1.001 | 0.004 | 13 |
| Japan | -0.292 | 0.902 | -1.402 | 0.007 | 32 |
| Netherlands | -0.666 | 0.801 | 0.889 | 0.003 | 8 |
| New Zealand | -0.496 | 1.055 | 0.086 | 0.007 | 10 |
| Norway | -0.491 | 1.366 | -0.171 | 0.025 | 16 |
| Portugal | -0.382 | 1.101 | -0.833 | 0.007 | 12 |
| Spain | -0.721 | 0.835 | -1.051 | 0.010 | 13 |
| **Sweden** | **-0.001** | **-206.064** | **-239.938** | **0.018** | **23** |
| U.S. | -0.218 | 0.965 | -0.449 | 0.013 | 32 |
| U.K. | -0.487 | 1.041 | -0.789 | 0.008 | 32 |

Sweden's phi is -0.001. Every other country has phi between -0.07 and -0.72.

## Why phi near Zero Causes Extreme Long-Run Estimates

The long-run parameter is computed as:

```
theta = -beta / phi
```

When phi is close to zero, this ratio explodes. For Sweden:
- beta_inc = 0.206, phi = -0.001 --> theta_inc = -0.206 / -0.001 = -206
- beta_inf = 0.240, phi = -0.001 --> theta_inf = -0.240 / -0.001 = -240

These are nonsensical values. An income elasticity of -206 means a 1% increase in income causes a 206% decrease in consumption. The true long-run relationship exists but is unidentifiable from this short time series because the error correction mechanism is essentially zero — consumption and income are near-random-walk cointegrated, and the ECM regression cannot distinguish the adjustment speed from zero noise.

## The Intuition: What phi near Zero Means Economically

phi measures how fast consumption returns to its long-run equilibrium after a shock. If phi = -0.3, then 30% of the gap closes each year. If phi = -0.001, the gap closes at 0.1% per year — effectively no adjustment at all. With only 23 annual observations, a speed of adjustment this slow is statistically indistinguishable from a unit root (no cointegration). The long-run relationship becomes unidentified: you cannot estimate where the system is heading if it never gets there.

## How It Propagates Through Estimators

### Mean Group (MG)

MG simply averages the country-specific estimates:

```
theta_bar = (1/N) * sum(theta_i)
```

One extreme value dominates:
- **With Sweden (N=15)**: theta_inc = -12.88 (se=13.80)
- **Without Sweden (N=14)**: theta_inc = 0.918 (se=0.066)

A single outlier swings the average by 13.8 and inflates the standard error by a factor of 200. The MG estimate is meaningless.

### Pooled Mean Group (PMG)

PMG constrains theta to be common across groups. Because it estimates via concentrated ML (weighting by 1/sigma^2), groups with large residual variance get downweighted. However, Sweden's sigma (0.018) is not unusually large — the problem is phi, not sigma. The PMG iteration attempts to find a common theta such that every group's phi is well-defined, but Sweden's near-zero phi makes the information matrix singular (the phi_Sweden row/column has near-zero entries). The Newton-Raphson update and the final SE computation both fail:

```
Error: system is computationally singular: reciprocal condition number = 2.27e-17
```

**Without Sweden (N=14)**: PMG converges normally, theta_inc = 0.919 (se=0.014).

### Dynamic Fixed Effects (DFE)

DFE pools all groups into a single within-regression, estimating a common phi. Sweden's near-zero phi is averaged into the pool and has modest impact because it's one group of 15. The DFE estimate (theta_inc = 0.935, phi = -0.179) is reasonable.

### Static Fixed Effects (SFE)

SFE imposes phi = -1 by construction, so the near-zero phi problem never arises. SFE is immune to this class of outlier but pays for it with potential bias from the wrong dynamics.

## Summary: Estimator Sensitivity to Outliers

| Estimator | Mechanism | Outlier Sensitivity |
|-----------|-----------|---------------------|
| **MG** | Simple average of group theta | **Catastrophic** — one extreme theta dominates |
| **PMG** | Concentrated ML, weighted by 1/sigma^2 | **Numerically singular** when any phi near 0 |
| **DFE** | Pooled within-regression | **Modest** — outlier averaged into common phi |
| **SFE** | Imposes phi = -1 | **Immune** — phi never estimated |

## General Principles for Outlier Handling

### 1. Diagnose Before Dropping

Always run MG first and inspect the group-specific estimates. Look for:
- **phi near zero**: signals no cointegration / unidentified LR relationship
- **theta implausibly large**: symptom of phi near zero
- **sigma much larger than other groups**: genuinely noisy data
- **T very small**: few observations make all estimates unreliable

The MG is the diagnostic tool — its whole purpose is to reveal heterogeneity.

### 2. Understand Why the Outlier Exists

Sweden's outlier is not "bad data." The paper explains: the SBC chose ARDL(1,1,0) for Sweden rather than ARDL(1,1,1). Under ARDL(1,1,0), Sweden's estimates are reasonable (phi = -0.35, theta_inc = 0.94). The outlier is a **specification problem**, not a data problem. The ARDL(1,1,1) model is wrong for Sweden because the additional lagged inflation term absorbs all the adjustment dynamics, leaving phi ≈ 0.

This matters because:
- Dropping Sweden "fixes" the numbers but hides the fact that the common lag specification is inappropriate
- Using SBC-chosen lags per country avoids the problem entirely (Table 4 in the paper)
- The PMG's robustness to specification comes precisely from pooling — the paper shows PMG barely changes when Sweden is dropped

### 3. Quantify the Impact

Run the analysis with and without the suspected outlier:

| | MG (N=15) | MG (N=14) | PMG (N=14) |
|--|-----------|-----------|------------|
| theta_inc | -12.88 (13.80) | 0.918 (0.066) | 0.919 (0.014) |
| theta_inf | -16.64 (15.95) | -0.686 (0.176) | -0.994 (0.077) |
| phi | -0.364 (0.051) | -0.390 (0.047) | -0.266 (0.057) |

Key observations:
- MG is completely destroyed; dropping one group changes theta by 13.8
- PMG (N=14) barely differs from PMG (N=24, NDI, Table 1: 0.904) — the pooling provides robustness
- phi averages are stable because Sweden's phi (-0.001) is close to zero and has small absolute effect on the average

### 4. Decision Framework

```
Is phi_i near zero for group i?
├── Yes: Is it due to specification (wrong lag order)?
│   ├── Yes: Use SBC to choose lags, or use a simpler ARDL spec
│   └── No (genuine unit root): Exclude group from MG/PMG, note it
└── No: Keep group in estimation

Is the MG average driven by one group?
├── Yes: Report both with and without, prefer PMG
└── No: MG is reliable
```

### 5. PMG as a Natural Robustness Check

The paper's main argument for PMG over MG is precisely this: PMG constrains theta to be common, so one group's extreme OLS estimate cannot hijack the pooled estimate. The concentrated ML weights each group's contribution by phi^2/sigma^2, which naturally downweights groups with small phi. However, when phi is *exactly* zero, the information matrix becomes singular rather than just downweighted — this is a numerical rather than statistical failure.

A practical fix (not implemented in the original GAUSS code) would be to detect and exclude groups with |phi| below a threshold (e.g., 0.01) from the PMG estimation, similar to how the code already handles the phi = -1 case.

## Reproducing This Analysis

```r
source("panel_helpers.R")

# Build PDI data.frame
# (see replicate_tables.R for the full data loading code)

# Run MG to diagnose
mg <- mg_panel(df_pdi, y="lpc", x=c("lpdi","dp"), unit="country", time="year",
               plag=1, qlag=c(1,1))
print(mg$group_phi)  # Sweden will have phi near 0

# Compare with and without
mg_no_sweden <- mg_panel(df_no_sweden, ...)
pmg_no_sweden <- pmg_panel(df_no_sweden, ...)
```
