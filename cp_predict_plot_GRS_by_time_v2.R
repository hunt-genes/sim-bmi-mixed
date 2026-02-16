## ---- Predicted BMI by time point using emmeans ----

# Packages
library(emmeans)
library(ggplot2)
library(lspline)

# 1) Representative age and matching spline values
age_ref <- median(df$age)

# IMPORTANT: use the same knots I used when building df
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

# Label groups 
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
