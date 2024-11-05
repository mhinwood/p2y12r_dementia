# Controlled Direct Effect Analysis Script
library(tidyverse)
library(survival)
source("utils.R")  # Load utility functions

# Define time points (in quarters)
quarters_of_interest <- c(4, 8, 20)

# Treatment weights
trt_den <- glm(
  trt_3m ~ aspirin + TIA + stroke_treat + Depression + sex + inflammation + 
    ACS + incomeq5 + bs(NIHSS) + aware_arrival + statin + 
    education_collapsed + year_fct + bs(age),
  data = long,
  family = binomial
)

trt_num <- glm(trt_3m ~ 1, family = binomial(link = "logit"), data = long)

# Calculate treatment weights
long <- long %>% 
  mutate(
    p_denom_trt = predict(trt_den, type = "response"),
    trt_num_pred = predict(trt_num, type = "response"),
    w_trt_cde = ifelse(trt_3m == 1, 
                       trt_num_pred / p_denom_trt,
                       (1 - trt_num_pred) / (1 - p_denom_trt))
  )


# Baseline competing event (death) weights
death_baseline <- glm(
  competing ~ trt_3m + statin + stroke_treat + smoking + sex + ACS + 
    incomeq5 + bs(NIHSS) + aware_arrival + education_collapsed + 
    year_fct + bs(age),
  data = long,
  family = quasibinomial
)

# Time-varying competing event weights
death_tv <- glm(
  competing ~ trt_3m * bs(Treatment_block, 3) + statin + stroke_treat + 
    smoking + sex + ACS + incomeq5 + bs(NIHSS) + aware_arrival + 
    education_collapsed + year_fct + bs(age) + BP_pre,
  data = long,
  family = quasibinomial
)
# Add time-varying weights
long <- long %>%
  mutate(p_denom_tv = predict(death_tv, type = "response")) %>%
  group_by(id) %>%
  mutate(
    w_tv = 1 / cumprod(1 - p_denom_tv),
    # Cap at 99th percentile
    w_tv = pmin(w_tv, quantile(w_tv, 0.99, na.rm = TRUE))
  ) %>%
  ungroup()

# Combine all weights
long <- long %>%
  mutate(
    combined_weights = w_trt * w_death * w_tv,
    weights_trimmed = pmin(combined_weights, quantile(combined_weights, 0.99))
  )

# Fit models
# Crude model
km_crude <- survfit(
  Surv(time = tstart, time2 = fuptime, event = outcome) ~ treatment,
  data = long,
  cluster = id
)

# Treatment-weighted only
km_trt <- survfit(
  Surv(time = tstart, time2 = fuptime, event = outcome) ~ treatment,
  data = long,
  cluster = id,
  weights = w_trt
)

# Treatment and baseline competing event weights
km_baseline <- survfit(
  Surv(time = tstart, time2 = fuptime, event = outcome) ~ treatment,
  data = long,
  cluster = id,
  weights = w_trt * w_death
)

# Fully adjusted model
km_full <- survfit(
  Surv(time = tstart, time2 = fuptime, event = outcome) ~ treatment,
  data = long,
  cluster = id,
  weights = weights_trimmed
)

# Extract risks
risks_crude <- extract_risks(km_crude)
risks_trt <- extract_risks(km_trt)
risks_baseline <- extract_risks(km_baseline)
risks_full <- extract_risks(km_full)

# Bootstrap confidence intervals
ci_crude <- bootstrap_km(
  data = long,
  n_samples = 200,
  time_points = quarters_of_interest,
  seed = 126
)

ci_full <- bootstrap_km(
  data = long,
  n_samples = 200,
  time_points = quarters_of_interest,
  weight_var = "weights_trimmed",
  seed = 127
)

# Print results
print("Crude Model Results:")
print(risks_crude)
print(ci_crude)

print("Fully Adjusted Model Results:")
print(risks_full)
print(ci_full)