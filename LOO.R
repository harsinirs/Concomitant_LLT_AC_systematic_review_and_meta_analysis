#------------------------------------------------------------------------------#
# Leave one out analysis for primary outcomes
#------------------------------------------------------------------------------#

# Function to fit brms model on data subset and extract summaries
fit_and_extract <- function(data_subset, full_data_size) {
  cat("Fitting model on data with", nrow(data_subset), "rows...\n")
  
  model <- brm(
    formula = log_IRR | se(SE_log_IRR) ~ 1 + (1 | Study),
    data = data_subset,
    family = gaussian(),
    prior = c(
      prior(normal(0, 1), class = "Intercept"),
      prior(cauchy(0, 0.5), class = "sd")
    ),
    chains = 4,
    iter = 2000,
    warmup = 1000,
    control = list(adapt_delta = 0.995),
    save_pars = save_pars(all = TRUE),
    refresh = 0  # suppress sampling progress output
  )
  
  # Extract posterior draws and rename for clarity
  post.samples <- as_draws_df(model) %>%
    rename(
      smd = b_Intercept,
      tau = sd_Study__Intercept
    )
  
  cat("Model fitting done. Extracting summaries...\n")
  
  # Posterior probability pooled effect < 0
  p_less_0 <- mean(post.samples$smd < 0)
  
  # Calculate posterior I²
  post.samples$I2 <- (post.samples$tau^2) / 
    (post.samples$tau^2 + mean(data_subset$SE_log_IRR^2)) * 100
  
  I2_q <- quantile(post.samples$I2, probs = c(0.025, 0.5, 0.975))
  
  # Skewness of the effect size posterior
  skew_val <- skewness(post.samples$smd)
  
  # Summaries for intercept including 95% credible intervals
  smd_mean <- mean(post.samples$smd)
  smd_lower <- quantile(post.samples$smd, 0.025)
  smd_upper <- quantile(post.samples$smd, 0.975)
  
  # Mean tau (between-study heterogeneity)
  tau_mean <- mean(post.samples$tau)
  
  # Study label: "All" if full data, else LOO name of excluded study
  study_label <- ifelse(
    nrow(data_subset) == full_data_size,
    "All",
    paste0("LOO: ", setdiff(unique(major$Study), unique(data_subset$Study)))
  )
  
  cat("Summaries for", study_label, "computed.\n\n")
  
  tibble(
    Study = study_label,
    Intercept = smd_mean,
    Lower_CI = smd_lower,
    Upper_CI = smd_upper,
    Skewness = skew_val,
    P_RR_less_0 = p_less_0,
    Tau = tau_mean,
    I2_lower = I2_q[1],
    I2_median = I2_q[2],
    I2_upper = I2_q[3]
  )
}

# Get study names for LOO 
#use dataset of outcome of interest
studies <- unique(vte_data$Study)
full_data_size <- nrow(vte_data)

cat("Starting full model fit with all data...\n")

# Fit full model once (no study left out)
full_result <- fit_and_extract(vte_data, full_data_size)

# Initialize list to store LOO results
loo_results <- vector("list", length(studies))

# Loop over studies
for (i in seq_along(studies)) {
  cat("Running LOO iteration leaving out study:", studies[i], "\n")
  
  # Subset data excluding current study
  loo_data <- subset(vte_data, Study != studies[i])
  
  # Fit model and extract summaries
  loo_results[[i]] <- fit_and_extract(loo_data, full_data_size)
}

cat("LOO fitting completed. Combining results...\n")

# Combine all results: full data + LOO iterations
all_results <- bind_rows(full_result, bind_rows(loo_results))

# Optional: reorder columns
all_results <- all_results %>%
  select(Study, Intercept, Lower_CI, Upper_CI, Skewness, P_RR_less_0, Tau, I2_lower, I2_median, I2_upper)

cat("Final results:\n")
print(all_results)

# Optionally, you may want to save or export the results table here
# write.csv(all_results, "loo_results_summary.csv", row.names = FALSE)

