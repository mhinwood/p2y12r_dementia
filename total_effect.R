# Total Effect Analysis Script
library(tidyverse)
library(survival)
source("utils.R")  # Load utility functions

# Calculate IPTW
trt_den <- glm(
  treatment ~ aspirin + TIA + stroke_treat + Depression + sex + inflammation + 
    ACS + incomeq5 + bs(NIHSS) + aware_arrival + statin + 
    education_collapsed + year_fct + bs(age),
  data = wide,
  family = binomial(link = "logit")
)

wide <- wide %>%
  mutate(
    p_score = predict(trt_den, type = "response"),
    iptw = ifelse(treatment == 1, 1/p_score, 1/(1 - p_score))
  )

# Define time points (in days)
time_points <- c(365, 730, 1826)

# Fit crude models
outcome_crude <- survfit(
  Surv(time_outcome, outcome_event) ~ treatment,
  data = wide
)

competing_crude <- survfit(
  Surv(time_competing, competing_event) ~ treatment,
  data = wide
)

# Fit weighted models
outcome_weighted <- survfit(
  Surv(time_outcome, outcome_event) ~ treatment,
  data = wide,
  weights = iptw
)

competing_weighted <- survfit(
  Surv(time_competing, competing_event) ~ treatment,
  data = wide,
  weights = iptw
)

# Extract risk estimates
crude_outcome_risks <- extract_cif_metrics(outcome_crude, time_points)
crude_competing_risks <- extract_cif_metrics(competing_crude, time_points)
weighted_outcome_risks <- extract_cif_metrics(outcome_weighted, time_points)
weighted_competing_risks <- extract_cif_metrics(competing_weighted, time_points)

# Bootstrap confidence intervals
crude_outcome_ci <- bootstrap_cif_extraction(
  data = wide,
  time_var = "time_outcome",
  status_var = "outcome_event",
  trt_var = "treatment",
  time_points = time_points
)

weighted_outcome_ci <- bootstrap_cif_extraction(
  data = wide,
  time_var = "time_outcome",
  status_var = "outcome_event",
  trt_var = "treatment",
  weight_var = "iptw",
  time_points = time_points
)

# Print results
print("Crude Outcome Results:")
print(crude_outcome_risks)
print(crude_outcome_ci)

print("Weighted Outcome Results:")
print(weighted_outcome_risks)
print(weighted_outcome_ci)