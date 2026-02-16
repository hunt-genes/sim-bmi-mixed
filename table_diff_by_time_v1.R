library(emmeans)

# Turn off tests that require pbkrtest/lmerTest; still gives p-values for contrasts
emm_options(disable.tests = FALSE, lmerTest.limit = Inf)

# ------------------------------------------------------------
# 1. Compute marginal means for each yearcat × GRS combination
#    (Model-fixed effects only, re.form = NA)
#    Average over age and spline terms.
# ------------------------------------------------------------
emm_grs <- emmeans(
  model,
  specs = ~ GRS | yearcat,
  re.form = NA
)

# ------------------------------------------------------------
# 2. Compute contrasts: GRS4 - ref
#    Gives the difference at each time period
# ------------------------------------------------------------
contrast_results <- contrast(
  emm_grs,
  method = list("GRS4 - ref" = c(-1, 0, 0, 0, 1)),   # subtract ref from GRS4
  by = "yearcat",                                    # separate estimate per time point
  adjust = "none"
)

# ------------------------------------------------------------
# 3. Get full table with SE, CI, and p-values
# ------------------------------------------------------------
results_table <- summary(
  contrast_results,
  infer = c(TRUE, TRUE),   # include CI
  level = 0.95
)

# Print 
print(results_table)

# Tibble for export

library(dplyr)

table_clean <- results_table %>%
  as.data.frame() %>%
  select(
    yearcat,
    estimate,
    SE,
    df,
    asymp.LCL,
    asymp.UCL,
    p.value
  )

write.csv(table_clean, "GRS_difference_table.csv", row.names = FALSE)
