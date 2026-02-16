library(emmeans)
library(ggplot2)
library(lspline)
library(dplyr)

# ---------------------------------------------------------
# 1. Define age groups and representative ages
# ---------------------------------------------------------
age_groups <- data.frame(
  age_group = c("25-35", "35-50", "50-65", "65-80"),
  age_min   = c(25, 35, 50, 65),
  age_max   = c(35, 50, 65, 80)
)

age_groups$age_ref <- (age_groups$age_min + age_groups$age_max) / 2

# ---------------------------------------------------------
# 2. Generate spline values for each age reference
# ---------------------------------------------------------
knots <- c(20, 30, 40, 50, 60, 70)

age_groups$splines <- lapply(age_groups$age_ref, function(a) {
  S <- lspline(a, knots = knots)
  colnames(S) <- paste0("as", 1:ncol(S))
  as.list(as.data.frame(S))
})

# ---------------------------------------------------------
# 3. Specify GRS groups to compare
# ---------------------------------------------------------
grs_levels <- c("ref", "GRS4")

# Turn off hypothesis tests 
emm_options(disable.tests = TRUE, lmerTest.limit = Inf)

# ---------------------------------------------------------
# 4. Create predictions for each age group × GRS × yearcat
# ---------------------------------------------------------
emm_list <- list()

for (i in seq_len(nrow(age_groups))) {
  
  spline_vals <- age_groups$splines[[i]]
  
  emm_i <- emmeans(
    model,
    specs = ~ yearcat | GRS,
    at = c(list(GRS = grs_levels), spline_vals),
    re.form = NA
  )
  
  tmp <- as.data.frame(summary(emm_i, infer = c(TRUE, TRUE)))
  tmp$age_group <- age_groups$age_group[i]
  
  emm_list[[i]] <- tmp
}

emm_df <- bind_rows(emm_list)

# ---------------------------------------------------------
# 5. Standardize CI column names 
# ---------------------------------------------------------
if ("lower.CL" %in% names(emm_df)) {
  emm_df$LCL <- emm_df$lower.CL
  emm_df$UCL <- emm_df$upper.CL
} else if ("asymp.LCL" %in% names(emm_df)) {
  emm_df$LCL <- emm_df$asymp.LCL
  emm_df$UCL <- emm_df$asymp.UCL
} else {
  stop("No confidence interval columns found in emmeans output.")
}

# ---------------------------------------------------------
# 6. Create descriptive GRS labels
# ---------------------------------------------------------
emm_df$GRS_group <- factor(
  ifelse(emm_df$GRS == "GRS4",
         "Top fifth (most susceptible)",
         "Bottom fifth (least susceptible)"),
  levels = c("Bottom fifth (least susceptible)",
             "Top fifth (most susceptible)")
)

# ---------------------------------------------------------
# 7. Final plot
# ---------------------------------------------------------
ggplot(
  emm_df,
  aes(
    x = yearcat,
    y = emmean,
    color = age_group,
    linetype = GRS_group,
    group = interaction(age_group, GRS_group)
  )) +
  
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = LCL, ymax = UCL), width = 0.12, alpha = 0.6) +
  
  scale_color_brewer(palette = "Dark2") +
  scale_linetype_manual(values = c("solid", "dashed")) +
  
  labs(
    title = "Predicted BMI by time point, stratified by age group and GRS",
    x = "Time point (year category)",
    y = "Predicted BMI",
    color = "Age group",
    linetype = "GRS group"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )
