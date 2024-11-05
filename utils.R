# Functions for extracting risks and bootstrapping estimates

# Extract cumulative incidence function metrics
extract_cif_metrics <- function(model, time_points) {
  summary_data <- summary(model, times = time_points)
  
  result <- data.frame(
    time_point = time_points,
    risk_control = NA_real_,
    risk_treatment = NA_real_,
    rd = NA_real_,
    rr = NA_real_
  )
  
  n_strata <- length(unique(summary_data$strata))
  if (n_strata != 2) {
    warning("Expected 2 strata, but found ", n_strata)
    return(result)
  }
  
  for (i in seq_along(time_points)) {
    indices <- which(summary_data$time == time_points[i])
    if (length(indices) == 2) {
      control_est <- 1 - summary_data$surv[indices[1]]
      treatment_est <- 1 - summary_data$surv[indices[2]]
      
      result$risk_control[i] <- control_est * 100
      result$risk_treatment[i] <- treatment_est * 100
      result$rd[i] <- (treatment_est - control_est) * 100
      result$rr[i] <- ifelse(control_est == 0, NA, treatment_est / control_est)
    }
  }
  
  return(result)
}

# Extract risks for time-varying analysis
extract_risks <- function(fit) {
  # Get times and survival estimates
  times <- fit$time
  surv <- fit$surv
  
  # Split survival estimates for treatment groups
  n_times <- length(times)
  surv_0 <- surv[1:n_times]
  surv_1 <- surv[(n_times + 1):(2 * n_times)]
  
  # Calculate risks
  risk_0 <- 1 - surv_0
  risk_1 <- 1 - surv_1
  
  # Calculate differences and ratios
  risk_diff <- risk_1 - risk_0
  risk_ratio <- risk_1 / risk_0
  
  # Return data frame
  data.frame(
    quarter = times,
    risk_0 = risk_0,
    risk_1 = risk_1,
    risk_diff = risk_diff,
    risk_ratio = risk_ratio
  )
}

# Bootstrap function for total effects
bootstrap_cif_extraction <- function(data, time_var, status_var, trt_var, weight_var = NULL, 
                                     n_boot = 200, time_points, seed = 123) {
  set.seed(seed)
  
  boot_results <- replicate(n_boot, {
    boot_indices <- sample(nrow(data), replace = TRUE)
    boot_data <- data[boot_indices, ]
    
    # Fit model with or without weights
    if(is.null(weight_var)) {
      fit <- survfit(as.formula(paste("Surv(", time_var, ",", status_var, ") ~", trt_var)),
                     data = boot_data)
    } else {
      fit <- survfit(as.formula(paste("Surv(", time_var, ",", status_var, ") ~", trt_var)),
                     data = boot_data, weights = boot_data[[weight_var]])
    }
    
    extract_cif_metrics(fit, time_points)
  }, simplify = FALSE)
  
  boot_combined <- do.call(rbind, boot_results)
  
  ci_results <- lapply(time_points, function(tp) {
    subset_data <- boot_combined[boot_combined$time_point == tp, ]
    list(
      time_point = tp,
      risk_control_ci = quantile(subset_data$risk_control, probs = c(0.025, 0.975), na.rm = TRUE),
      risk_treatment_ci = quantile(subset_data$risk_treatment, probs = c(0.025, 0.975), na.rm = TRUE),
      rd_ci = quantile(subset_data$rd, probs = c(0.025, 0.975), na.rm = TRUE),
      rr_ci = quantile(subset_data$rr, probs = c(0.025, 0.975), na.rm = TRUE)
    )
  })
  
  return(ci_results)
}

# Bootstrap function for controlled direct effect
bootstrap_km <- function(data, n_samples = 200, time_points, weight_var = NULL, seed = 123) {
  set.seed(seed)
  
  # Store results
  results <- vector("list", n_samples)
  
  # Get unique clusters
  clusters <- unique(data$id)
  
  # Progress bar
  pb <- txtProgressBar(min = 0, max = n_samples, style = 3)
  
  for(i in 1:n_samples) {
    # Sample clusters with replacement
    sampled_clusters <- sample(clusters, replace = TRUE)
    boot_data <- data[data$id %in% sampled_clusters, ]
    
    # Fit model
    if(is.null(weight_var)) {
      fit <- survfit(
        Surv(time = tstart, time2 = fuptime, event = outcome) ~ treatment,
        data = boot_data,
        cluster = id
      )
    } else {
      fit <- survfit(
        Surv(time = tstart, time2 = fuptime, event = outcome) ~ treatment,
        data = boot_data,
        cluster = id,
        weights = boot_data[[weight_var]]
      )
    }
    
    # Extract risks
    risks <- extract_risks(fit)
    results[[i]] <- risks[risks$quarter %in% time_points, ]
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  # Calculate CIs
  ci_results <- lapply(time_points, function(tp) {
    subset_data <- do.call(rbind, lapply(results, function(x) x[x$quarter == tp, ]))
    list(
      quarter = tp,
      risk_0_ci = quantile(subset_data$risk_0, probs = c(0.025, 0.975), na.rm = TRUE),
      risk_1_ci = quantile(subset_data$risk_1, probs = c(0.025, 0.975), na.rm = TRUE),
      risk_diff_ci = quantile(subset_data$risk_diff, probs = c(0.025, 0.975), na.rm = TRUE),
      risk_ratio_ci = quantile(subset_data$risk_ratio, probs = c(0.025, 0.975), na.rm = TRUE)
    )
  })
  
  return(ci_results)
}