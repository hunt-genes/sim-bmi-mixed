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

1. Predicted BMI by time point
This repository includes an example visualization that uses emmeans to compute marginal (population-level) predicted BMI for two GRS groups (bottom fifth and top fifth) across survey time points, holding age fixed at a representative value (the median age in the sample). The plot shows estimated means and confidence intervals for each time point and GRS group.

Purpose
- To illustrate how predicted (population-level) BMI differs between high- and low-risk GRS groups across survey periods.
- To show time-varying genetic susceptibility while controlling for age via spline terms (held constant at `age_ref`).

Key points about the code
- Use the same spline knot locations used when building `df` (this is critical so that the spline basis for predictions matches the model input).
- Predictions use `re.form = NA` in `emmeans()` so the plot shows the marginal (fixed-effects) population-level predictions, not subject-specific predictions.
- Age is held fixed at `age_ref` (here chosen as the median age). This isolates the predicted difference in BMI by GRS and year at a representative age.
- The code requests predictions for two GRS levels (`ref` and `GRS4`) and then labels them as "Bottom fifth" vs "Top fifth".
- The code handles possible differences in the names for CI columns produced by emmeans (e.g., `lower.CL` vs `asymp.LCL`).
- `emm_options(lmerTest.limit = Inf, disable.tests = TRUE)` disables lmerTest-related tests that can otherwise cause spurious messages; the exact option usage may depend on your emmeans and lmerTest versions.

How to use 
```r
## ---- Predicted BMI by time point using emmeans ----
## Run after your model is fitted

# Packages
library(emmeans)
library(ggplot2)
library(lspline)

# 1) Representative age and matching spline values
age_ref <- median(df$age)

# IMPORTANT: use the same knots you used when building df
spline_knots <- c(20, 30, 40, 50, 60, 70)

Xref <- lspline(age_ref, knots = spline_knots)
colnames(Xref) <- paste0("as", 1:ncol(Xref))
Xref <- as.list(as.data.frame(Xref))  # list: as1=..., as2=..., etc.

# 2) Get predictions for the two GRS groups across all yearcat levels
#    This holds age/splines fixed at age_ref; averages over random effects
grs_levels <- c("ref", "GRS4")

# Disable tests
emm_options(lmerTest.limit = Inf, disable.tests = TRUE)

# Predictions
emm <- emmeans(
  model,
  specs = ~ yearcat | GRS,
  at = c(list(GRS = grs_levels), Xref),
  re.form = NA
)

# Ensure CI columns are included
emm_df <- as.data.frame(summary(emm, infer = c(TRUE, TRUE)))

# Label groups nicely
emm_df$GRS_group <- factor(
  ifelse(emm_df$GRS == "GRS4",
         "Top fifth (most susceptible)",
         "Bottom fifth (least susceptible)"),
  levels = c("Top fifth (most susceptible)", "Bottom fifth (least susceptible)")
)

# Standardize CI column names
if ("lower.CL" %in% names(emm_df)) {
  emm_df$LCL <- emm_df$lower.CL
  emm_df$UCL <- emm_df$upper.CL
} else if ("asymp.LCL" %in% names(emm_df)) {
  emm_df$LCL <- emm_df$asymp.LCL
  emm_df$UCL <- emm_df$asymp.UCL
} else {
  stop("No CI columns found in emm_df — cannot plot error bars.")
}

# Plot using LCL/UCL
ggplot(emm_df, aes(x = yearcat, y = emmean, color = GRS_group, group = GRS_group)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.15, alpha = 0.7) +
  scale_color_manual(
    values = c("Top fifth (most susceptible)" = "blue",
               "Bottom fifth (least susceptible)" = "red")
  ) +
  labs(
    title = sprintf("Predicted BMI by time point for top vs bottom fifth of GRS (age = %.1f)", age_ref),
    x = "Time point (year category)",
    y = "Predicted BMI",
    color = "GRS group"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")
```
<img width="1342" height="735" alt="image" src="https://github.com/user-attachments/assets/2823a0c4-3e81-4e4c-a086-f0b0a73fd5b3" />

Interpretation of the plotted results
- Each point shows the model-predicted mean BMI for a GRS group at a survey time point, at the representative age (`age_ref`).
- Error bars represent the confidence interval for the predicted mean.
- Because predictions are computed with `re.form = NA`, these are marginal (population-average) estimates, not conditional on particular subjects.
- The plot therefore communicates how the fixed-effects portion of the model predicts BMI across time for low vs high genetic risk, at the chosen age.



