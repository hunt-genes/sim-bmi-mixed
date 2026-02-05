# sim-bmi-mixed
Simulation of BMI by time with repeated measures in a mixed model

This repository contains R code to simulate a longitudinal dataset and fit a mixed-effects model for BMI using age splines, categorical survey-year periods, and a categorical genetic risk score (GRS). The simulation and model are intended as a reproducible example for exploring time-varying genetic effects and age trajectories in a repeated-measures setting.

## Overview

- We simulate data for n_id individuals, each measured n_obs times across an age window.
- Age is represented using piecewise linear splines to flexibly model non-linear age-BMI relationships.
- A categorical GRS with ordered levels is assigned at the individual level.
- Year categories (survey periods) capture secular/time effects and interact with GRS to produce deterministic, monotonic genetic effects that increase over later survey periods.
- A linear mixed model (lmer) is fit with random intercepts and random slopes for age (grouped by `id`) plus multiple interaction terms between splines, year categories, and GRS.

## Files

- `simulate_and_fit.R` — creates the simulated longitudinal dataset and fits the mixed-effects model.
- `README.md` — this file.

## Requirements

- R >= 4.0 
- Packages:
  - lme4
  - lspline

Install packages (one-time):
```r
install.packages("lme4")
install.packages("lspline")
```

Load in script:
```r
library(lme4)
library(lspline)
```

## How to run

1. Save the provided R code into a file, e.g. `simulate_and_fit.R`.
2. Open R or RStudio and run:
```r
source("simulate_and_fit.R")
```
3. The script will:
   - create a simulated dataset `df` with columns described below,
   - fit a mixed-effects model with `lmer()`,
   - print `summary(model)` and aggregated means by `yearcat` and `GRS`.

Set a seed for reproducibility.

## Simulated data layout

The dataset `df` created by the script includes (not exhaustive):

- `id` — integer subject identifier (1..n_id)
- `age` — numeric age at observation
- `yearcat` — factor with levels: `66_69`, `84_86`, `95_97`, `00_01`, `06_08`
- `GRS` — factor genetic risk score categories: `ref`, `GRS1`, `GRS2`, `GRS3`, `GRS4`
- `as1`..`as7` — piecewise linear spline basis columns for age (from `lspline`)
- `u0`, `u1` — simulated per-subject random intercept and slope (used to generate BMI)
- `BMI` — simulated outcome (body mass index)

## Model specification used

The model fit in the script is:

BMI ~ yearcat + GRS +
  as1 + as2 + as3 + as4 + as5 + as6 + as7 +
  yearcat:GRS +
  (as1+as2+as3+as4+as5+as6+as7):GRS +
  (as1+as2+as3+as4+as5+as6+as7):yearcat +
  (1 + age | id)

Key points:
- Fixed effects include main effects of survey `yearcat`, `GRS`, and all age-spline components.
- Interaction terms allow the effect of GRS to vary by year period and the spline shapes to interact with GRS and year period.
- Random effects include a random intercept and random slope on continuous `age` for each subject (`id`).

The model is fit with REML = FALSE (maximum likelihood) in the example script:
```r
model <- lmer(..., data = df, REML = FALSE)
```

## Interpretation guidance

- Coefficients for `yearcat` represent average shifts in BMI across survey periods (reference level is `66_69`).
- `GRS` coefficients represent average differences in BMI across genetic-risk categories (reference is `ref`) when other variables are at reference.
- `yearcat:GRS` interaction coefficients show how the GRS effect differs by survey period. In the simulated data these interactions were constructed to increase the genetic effect with later periods.
- Spline coefficients (`as1..as7`) together define the age–BMI relationship; interpreting single spline coefficients in isolation is not intuitive — visualize predicted curves instead.
- Random-effect variance components indicate between-subject variability in intercepts and age slopes.

## Recommended diagnostics and checks

- Check model convergence and warnings printed by `lmer()`. If convergence warnings occur, consider:
  - simplifying random-effects structure,
  - rescaling predictors (e.g., center age),
  - refitting using different optimizers (e.g., `control = lmerControl(optimizer = "bobyqa")`).
- Residual diagnostics:
  - plot residuals vs fitted,
  - check QQ-plot of residuals,
  - check heteroscedasticity across `age`, `yearcat`, and `GRS`.
- Random effects:
  - plot subject-specific fitted lines or random intercepts/slopes.
- Multicollinearity:
  - when using many spline bases and interactions, check correlations/VIF-like diagnostics for fixed effects (VIF is less well-defined with lmer; consider inspecting the fixed-effects design matrix).

## Visualizations



