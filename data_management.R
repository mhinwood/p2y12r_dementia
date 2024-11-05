# Data preparation script for competing risks survival analysis
# Creates both wide and long format datasets

library(tidyverse)
library(survival)
library(splines)

# Set administrative censoring date
admin_censor_date <- as.Date("2018-12-31")

# Prepare wide format dataset
wide <- raw_data %>%
  # Calculate time-to-event variables
  mutate(
    # Time to outcome
    time_outcome = case_when(
      !is.na(outcome_date) ~ as.numeric(outcome_date - index_date),
      !is.na(competing_date) ~ as.numeric(competing_date - index_date),
      TRUE ~ as.numeric(admin_censor_date - index_date)
    ),
    # Time to competing event
    time_competing = case_when(
      !is.na(competing_date) ~ as.numeric(competing_date - index_date),
      TRUE ~ as.numeric(admin_censor_date - index_date)
    ),
    # Event indicators
    outcome_event = ifelse(!is.na(outcome_date), 1, 0),
    competing_event = ifelse(!is.na(competing_date), 1, 0)
  ) %>%
  # Remove negative follow-up times
  filter(time_outcome >= 0, time_competing >= 0)

# Prepare long format data for CDE
long <- wide %>%
  group_by(id) %>%
  mutate(
    # Calculate quarters
    time_period = case_when(
      !is.na(outcome_date) ~ pmin(ceiling(time_outcome/90), 20),
      !is.na(competing_date) ~ pmin(ceiling(time_competing/90), 20),
      TRUE ~ 20
    ),
    period = row_number()
  ) %>%
  filter(period <= time_period) %>%
  mutate(
    outcome = ifelse(outcome_event == 1 & period == time_period, 1, 0),
    competing = ifelse(competing_event == 1 & period == time_period & outcome == 0, 1, 0),
    tstart = lag(period, default = 0),
    fuptime = ifelse(outcome == 1 | competing == 1, period,
                     ifelse(period <= time_period, period, time_period))
  ) %>%
  ungroup()

# Save prepared datasets
saveRDS(wide, "wide_data.rds")
saveRDS(long, "long_data.rds")