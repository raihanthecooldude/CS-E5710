########### Installing and loading packages #########
if (!require(tidybayes)) {
  install.packages("tidybayes")
  library(tidybayes)
}
if (!require(brms)) {
  install.packages("brms")
  library(brms)
}
if (!require(metadat)) {
  install.packages("metadat")
  library(metadat)
}
if(!require(cmdstanr)){
  install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
  library(cmdstanr)
}
cmdstan_installed <- function(){
  res <- try(out <- cmdstanr::cmdstan_path(), silent = TRUE)
  !inherits(res, "try-error")
}
if(!cmdstan_installed()){
  install_cmdstan()
}
library(dplyr)



########### Loading data ###########
data <- read.csv("maternal-mortality-2021-all.csv")
View(data)

# Add a Region_Country column
data <- data |>
  dplyr::mutate(
    Country = factor(Country),
    Region = factor(Region),
    Region_Country = interaction(Region, Country, drop = TRUE, sep = "_")
  )

##### Specify priors #####
test_hierarchical_formula  <- bf(
  DeathRate ~ 1 + (1 | Region_Country),
  family = "gaussian",
  center = FALSE
)

test_priors <- get_prior(test_hierarchical_formula, data = data)
View(test_priors)

##### Calculate priors #####
# b-country
country_means <- data |> group_by(Country) |>
  summarize(MeanRate = mean(DeathRate, na.rm = TRUE))
group_variance <- var(country_means$MeanRate)
group_sd <- sqrt(group_variance)
group_sd # 1.19

######## Separate model ########
separate_formula <- bf(
  formula = DeathRate ~ 0 + Country,
  family = "gaussian"
)  

separate_fit <- brm(
  formula = separate_formula,
  data = data,
  iter = 2000,      # Number of iterations
  warmup = 1000,    # Number of warmup iterations
  chains = 4,        # Number of chains
  seed = 2024,
  prior = prior(normal(3.639, 3.06))
)

summary(separate_fit)


######## Pooled model #######

pooled_formula <- bf(
  formula = DeathRate ~ 1 + Country,
  family = "gaussian",
  center = FALSE
)

pooled_fit <- brm(
  formula = pooled_formula,
  data = data,
  iter = 2000,      # Number of iterations
  warmup = 1000,    # Number of warmup iterations
  chains = 4,        # Number of chains
  seed = 2024,
  prior = prior(normal(3.64, 3.06))
)

summary(pooled_fit)

# pp_check create the default density overlay plot
pp_check(pooled_fit)
# add the leave-one-out CV criterion to the fit object
pooled_fit <- add_criterion(
  pooled_fit,
  criterion = "loo"
)
# view the output
loo(pooled_fit)


######## Partially pooled model within countries ########
pooled_countries_formula <- bf(
  formula = DeathRate ~ 1 + (1 | Country),
  family = "gaussian",
  center = FALSE
)
pooled_countries_priors <- c(
  prior(normal(3.64, 3.06), class = "b", coef = "Intercept"),  # Population-level intercept
  prior(normal(0, 1.19), class = "sd", group = "Country"),    # Standard deviation of varying intercepts for Country
  prior(normal(0, 3), class = "sigma") # variability within countries
)
pooled_countries_fit <- brm(
  formula = pooled_countries_formula,
  data = data,
  iter = 2000,      # Number of iterations
  warmup = 1000,    # Number of warmup iterations
  chains = 4,        # Number of chains
  seed = 2024,
  prior = pooled_countries_priors
)
summary(pooled_countries_fit)

pp_check(pooled_countries_fit)
pooled_countries_fit <- add_criterion(
  pooled_countries_fit,
  criterion = "loo", moment_match = TRUE
)
loo(pooled_countries_fit)

# Predictions for a new country in existing regions:
new_data_seychelles <- data.frame(
  Region = factor("Africa"),
  Country = factor("Seychelles"),
  Region_Country = factor("Africa_Seychelles")
)
data_post_epred1 <- posterior_epred(
  pooled_countries_fit,
  newdata = new_data_seychelles,
  allow_new_levels = TRUE
)
# Calculate the mean prediction
mean_prediction1 <- apply(data_post_epred1, 2, mean)
# Print the mean prediction
mean_prediction1

# Predictions for a new country in non-existing regions:
new_data_mars <- data.frame(
  Region = factor("Mars"),
  Country = factor("Space"),
  Region_Country = factor("Space_Mars")
)
data_post_epred0 <- posterior_epred(
  pooled_countries_fit,
  newdata = new_data_mars,
  allow_new_levels = TRUE
)
mean_prediction0 <- apply(data_post_epred0, 2, mean)
mean_prediction0

######## Partially pooled model for countries within region ########

# class sd - group region
data2020 <- data <- read.csv("maternal-mortality-2020.csv")
region_means <- data2020 %>%
  group_by(Region) %>%
  summarize(mean_rate = mean(DeathRate, na.rm = TRUE))
sd_region <- sd(region_means$mean_rate)
sd_region # approximately 0.6

# Model
pooled_countries_regions_formula <- bf(
  formula = DeathRate ~ 1 + (1 | Region) + (1 | Country),
  family = "gaussian",
  center = FALSE
)
pooled_region_countries_priors <- c(
  prior(normal(3.64, 3.06), class = "b", coef = "Intercept"),  # Population-level intercept
  prior(normal(0, 1.19), class = "sd", group = "Country"),    # Standard deviation of varying intercepts for Country
  prior(normal(0, 3), class = "sigma") # variability within countries
  #prior(normal(0, 0.6), class = "sd", group = "Region")
)
pooled_countries_regions_fit <- brm(
  formula = pooled_countries_regions_formula,
  data = data,
  iter = 2000,      # Number of iterations
  warmup = 1000,    # Number of warmup iterations
  chains = 4,        # Number of chains
  seed = 2024,
  prior = pooled_region_countries_priors
)
summary(pooled_countries_regions_fit)

# pp_check create the default density overlay plot
pp_check(pooled_countries_regions_fit)

# add the leave-one-out CV criterion to the fit object
pooled_countries_regions_fit <- add_criterion(
  pooled_countries_regions_fit,
  criterion = "loo",
  moment_match = TRUE
)

# view the output
loo(pooled_countries_regions_fit)



# Predictions for a new country in existing regions:
data_post_epred3 <- posterior_epred(
  pooled_countries_regions_fit,
  newdata = new_data_seychelles,
  allow_new_levels = TRUE
)
# Calculate the mean prediction
mean_prediction3 <- apply(data_post_epred3, 2, mean)
# Print the mean prediction
mean_prediction3

# Predictions for a new country in non-existing regions:
data_post_epred4 <- posterior_epred(
  pooled_countries_regions_fit,
  newdata = new_data_mars,
  allow_new_levels = TRUE
)
mean_prediction4 <- apply(data_post_epred4, 2, mean)
mean_prediction4

########## Model Comparison #########

loo_compare(pooled_fit, pooled_countries_fit, pooled_countries_regions_fit)

