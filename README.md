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

2. Predicted BMI by time point, stratified by age group and GRS

Purpose
- To illustrate how predicted (population-level) BMI differs between high- and low-risk GRS groups across survey periods.
- To show time-varying genetic susceptibility while controlling for age via spline terms (held constant at `age_ref`).

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

<img width="652" height="456" alt="image" src="https://github.com/user-attachments/assets/ca951f9f-cd49-4e8e-9b49-2430cb4c1603" />

Each point shows the model-predicted mean BMI for a GRS group at a survey time point, at the representative age (age_ref).
Error bars represent the confidence interval for the predicted mean.

- This code extends the single-age visualization to compare predicted BMI across multiple representative ages (age groups). For each age-group midpoint it computes the spline basis (so predictions match the model input), requests marginal model-predictions (emmeans with re.form = NA) for two GRS levels across all year categories, combines results, and plots predicted BMI with CI. The final plot uses color for age groups and linetype for the two GRS groups so you can simultaneously compare (1) how predictions change across time, (2) how they differ by GRS, and (3) how those differences vary by age.

Walkthrough:
  1. Define age groups and representative ages
     - CreateS 4 age bins (25–35, 35–50, 50–65, 65–80) and computeS an age_ref for each as the midpoint. These midpoints serve as representative ages for predictions.
  2. Generate spline values for each age reference
     - For each age_ref computes the same lspline basis used when fitting the model (knots = c(20,30,40,50,60,70).) Converting each spline result to a named list (as1..as7) lets pass them into emmeans via the at = argument so the model receives the correct predictor values.
  3. Specify GRS groups to compare
     - Choose the two categorical levels to compare: "ref" and "GRS4" (lowest vs highest risk).
  4. Create predictions for each age group × GRS × yearcat
     - For each age group call emmeans(model, specs = ~ yearcat | GRS, at = c(list(GRS = grs_levels), spline_vals), re.form = NA). This returns marginal predicted means (population-level) for each yearcat and GRS, holding the spline terms fixed at the representative age.
     - Results for each age group are collected in a list and then combined with bind_rows into one data.frame (emm_df).
  5. Standardize CI column names
     - emmeans output can name CI columns differently depending on package versions (lower.CL vs asymp.LCL).
  6. Create descriptive GRS labels
  7. Final plot: color = age_group, linetype = GRS_group



3. Predicted BMI trajectories across Age by Birth Cohort (population-level, averaged over GRS)

-Build a dense age grid (20–80).
-Create a prediction data.frame with all combinations of age × yearcat × GRS, preserving factor levels.
-Recreate the exact spline basis used to fit the model for each age in the grid.
-Predict the fixed-effect (population-level) BMI for each row using predict(..., re.form = NA).
-Compute correct standard errors for the fixed-effect predictions using Xβ and the fixed-effect covariance matrix (vcov).
-Average the predicted BMI and CIs over GRS (so the resulting curves are marginal with respect to GRS).
-Plot predicted BMI vs age for each yearcat with ribbons showing 95% CIs.

```r
# Step 1 — Make prediction grid with correct factor structure
age_grid <- seq(20, 80, by = 1)

pred_grid <- expand.grid(
  age = age_grid,
  yearcat = levels(df$yearcat),
  GRS = levels(df$GRS)
)

# Force same factor structure as model data
pred_grid$yearcat <- factor(pred_grid$yearcat,
                            levels = levels(df$yearcat))

pred_grid$GRS <- factor(pred_grid$GRS,
                        levels = levels(df$GRS))

# Step 2 — Recreate spline variables EXACTLY like training data
library(lspline)

make_splines <- function(age_vec){
  X <- lspline(age_vec, knots = c(30,40,50,60,70))
  colnames(X) <- paste0("as",1:ncol(X))
  
  # Ensure 7 spline columns
  if(ncol(X) < 7){
    X <- cbind(X, matrix(0, nrow(X), 7-ncol(X)))
    colnames(X) <- paste0("as",1:7)
  }
  
  as.data.frame(X[,1:7])
}

pred_grid <- cbind(pred_grid, make_splines(pred_grid$age))

# Step 3 — Population predictions (fixed effects only)
pred_grid$BMI <- predict(
  model,
  newdata = pred_grid,
  re.form = NA
)

# Step 4 — Correct standard errors
Xmat <- model.matrix(
  lme4::nobars(formula(model)),
  model.frame(lme4::nobars(formula(model)), pred_grid)
)

beta <- fixef(model)
vc <- vcov(model)

# Keep only columns that exist in the model
Xmat <- Xmat[, names(beta), drop = FALSE]

# Ensure correct column order
Xmat <- Xmat[, names(beta)]

pred_grid$SE <- sqrt(diag(Xmat %*% vc %*% t(Xmat)))

pred_grid$lower <- pred_grid$BMI - 1.96 * pred_grid$SE
pred_grid$upper <- pred_grid$BMI + 1.96 * pred_grid$SE

# Safety check
stopifnot(identical(colnames(Xmat), names(beta)))

# Step 5 — Average over GRS
library(dplyr)

plot_df <- pred_grid %>%
  group_by(age, yearcat) %>%
  summarise(
    BMI = mean(BMI),
    lower = mean(lower),
    upper = mean(upper),
    .groups = "drop"
  )

# Step 6 — Plot trajectories
library(ggplot2)

ggplot(plot_df,
       aes(age, BMI, colour = yearcat, fill = yearcat)) +
  
  geom_line(linewidth = 1.2) +
  
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.2,
    colour = NA
  ) +
  
  labs(
    x = "Age (years)",
    y = "Predicted BMI",
    colour = "Birth cohort",
    fill = "Birth cohort"
  ) +
  
  theme_minimal(base_size = 14)
```
   <img width="610" height="440" alt="image" src="https://github.com/user-attachments/assets/f530e04e-a1db-424a-9274-744dfdf48d8f" />

Interpretation

Each curve shows the model-predicted average BMI across age for one yearcat.
Ribbons show the approximate 95% confidence interval around the fixed-effect prediction (averaged over GRS).
This plot highlights how age trajectories differ by survey period / birth cohort while integrating over genetic risk.
