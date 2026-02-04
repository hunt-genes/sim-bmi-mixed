# Creates a longitudinal dummy data set
rm(list = ls())

# Packages
library(lme4)
#library(splines)
install.packages("lspline")   # once
library(lspline)

set.seed(1234)

# Study design
n_id  <- 400     # individuals
n_obs <- 5       # repeated measures per person

N <- n_id * n_obs

# Create longitudinal structure
df <- data.frame(
  id = rep(1:n_id, each = n_obs)
)

# Generate age
df$age <- rep(runif(n_id, 25, 75), each = n_obs) +
  rep(seq(0, 8, length.out = n_obs), times = n_id)

# Year categories
year_levels <- c("66_69", "84_86", "95_97", "00_01", "06_08")

df$yearcat <- sample(year_levels, N, replace = TRUE,
                     prob = c(0.25, 0.25, 0.15, 0.20, 0.15))

df$yearcat <- factor(df$yearcat, levels = year_levels)

##### Genetic risk score categories
# Individual-level GRS
GRS_id <- sample(c("ref", "GRS1", "GRS2", "GRS3", "GRS4"),
                 n_id, replace = TRUE)

# Assign to observations
df$GRS <- factor(GRS_id[df$id],
                 levels = c("ref", "GRS1", "GRS2", "GRS3", "GRS4"))


# Linear age splines (20 as reference, knots every 10 years)
knots <- c(30, 40, 50, 60, 70)

X <- lspline(df$age, knots = c(20,30,40,50,60,70))

colnames(X) <- paste0("as", 1:ncol(X))

df <- cbind(df, X)

# Random effects
u0 <- rnorm(n_id, 0, 2.0)   # random intercept SD
u1 <- rnorm(n_id, 0, 0.15)  # random slope SD

df$u0 <- u0[df$id]
df$u1 <- u1[df$id]

# Fixed-effect coefficients (arbitrary but realistic)
beta0 <- 23

beta_year <- c(0, 1.0, 1.5, 2.0, 2.5)
names(beta_year) <- year_levels

#beta_grs <- c(0, 0.8, 1.4, 2.0, 2.8)
#names(beta_grs) <- levels(df$GRS)

beta_grs <- c(
  ref  = 0,
  GRS1 = 2,
  GRS2 = 4,
  GRS3 = 6,
  GRS4 = 8
)

beta_spline <- runif(7, -0.3, 0.5)

#### Year-specific multipliers that increase over time
year_levels <- c("66_69", "84_86", "95_97", "00_01", "06_08")

year_multiplier <- c(
  "66_69" = 1.0,
  "84_86" = 1.4,
  "95_97" = 1.8,
  "00_01" = 2.3,
  "06_08" = 3.0
)

# Interaction coefficients
#beta_year_grs <- matrix(runif(4*4, -0.4, 0.4),
#                        nrow = 4, ncol = 4,
#                        dimnames = list(
#                          year = c("84_86","95_97","00_01","06_08"),
#                          grs  = c("GRS1","GRS2","GRS3","GRS4")
#                        ))

beta_year_grs <- matrix(0, nrow = 4, ncol = 4,
                        dimnames = list(
                          year = year_levels[-1],
                          grs  = c("GRS1","GRS2","GRS3","GRS4")
                        ))

for (y in year_levels[-1]) {
  for (g in c("GRS1","GRS2","GRS3","GRS4")) {
    beta_year_grs[y, g] <- year_multiplier[y] * beta_grs[g]
  }
}

#beta_spline_grs <- matrix(runif(7*4, -0.3, 0.3),
#                          nrow = 7, ncol = 4)

beta_spline_grs <- matrix(0, 7, 4)

beta_spline_year <- matrix(runif(7*4, -0.3, 0.3),
                           nrow = 7, ncol = 4)

# Construct BMI
BMI <- beta0

# year main effects (keep small)
beta_year <- c(
  "66_69" = 0,
  "84_86" = 0.5,
  "95_97" = 1.0,
  "00_01" = 1.5,
  "06_08" = 2.0
)

BMI <- BMI + beta_year[df$yearcat]

# age splines (keep moderate)
for (k in 1:7) {
  BMI <- BMI + 0.2 * df[[paste0("as", k)]]
}

# deterministic monotonic genetic effect
genetic_effect <- beta_grs[df$GRS] * year_multiplier[df$yearcat]
BMI <- BMI + genetic_effect

# random effects (small)
BMI <- BMI + 0.5 * df$u0 + 0.05 * df$u1 * df$age

# noise (small)
BMI <- BMI + rnorm(nrow(df), 0, 0.5)

df$BMI <- BMI

# Fit your model
model <- lmer(
  BMI ~ yearcat + GRS +
    as1 + as2 + as3 + as4 + as5 + as6 + as7 +
    yearcat:GRS +
    (as1+as2+as3+as4+as5+as6+as7):GRS +
    (as1+as2+as3+as4+as5+as6+as7):yearcat +
    (1 + age | id),
  data = df,
  REML = FALSE
)

# Verify
summary(model)

# Verify
aggregate(BMI ~ yearcat + GRS, df, mean)





