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

# Step 2 — Recreate spline variables 
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

