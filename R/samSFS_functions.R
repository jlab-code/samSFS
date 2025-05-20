
safe_likelihood_evaluation_joint <- function(
    base_grid,
    af_data,
    joint_mode = TRUE,
    top_n = 1,
    return_summary = TRUE,
    n_cores = 1
) {
  bin_cols <- grep("^bin_", names(base_grid), value = TRUE)
  epsilon <- 1e-10
  
  for (col in bin_cols) {
    af_data[[col]] <- as.numeric(af_data[[col]])
  }
  
  has_batch_id <- "batch_id" %in% names(af_data)
  if (!has_batch_id && joint_mode) {
    warning("Batch ID column is missing. Falling back to single-sample analysis.")
    joint_mode <- FALSE
  }
  
  if (joint_mode) {
    af_data$batch_id[is.na(af_data$batch_id) | af_data$batch_id == ""] <- NA
    if (all(is.na(af_data$batch_id))) {
      warning("Batch ID column is empty. Falling back to single-sample analysis.")
      joint_mode <- FALSE
    }
  }
  
  batches <- if (joint_mode) unique(na.omit(af_data$batch_id)) else unique(af_data$branch_name)
  
  evaluate_one_batch <- function(batch_id) {
    sample_ids <- if (joint_mode) {
      unique(af_data$branch_name[af_data$batch_id == batch_id])
    } else {
      batch_id
    }
    
    total_mutations <- sum(sapply(sample_ids, function(sid) {
      sum(as.numeric(af_data[af_data$branch_name == sid, bin_cols, drop = FALSE]))
    }))
    
    if (total_mutations == 0) {
      warning(sprintf("⚠ Skipping batch '%s': no mutations.", batch_id))
      return(NULL)
    }
    
    log_liks <- sapply(seq_len(nrow(base_grid)), function(i) {
      params <- base_grid[i, ]
      joint_loglik <- 0
      for (sid in sample_ids) {
        p_row <- as.numeric(params[bin_cols]) + epsilon
        p_row <- p_row / sum(p_row)
        obs_counts <- as.numeric(af_data[af_data$branch_name == sid, bin_cols, drop = FALSE])
        alpha_s <- sum(obs_counts) / total_mutations
        joint_loglik <- joint_loglik + alpha_s * sum(obs_counts * log(p_row))
      }
      joint_loglik
    })
    
    result_df <- base_grid[, c("rho", "m", "omega")]
    result_df$joint_log_likelihood <- log_liks
    result_df$batch_id <- batch_id
    
    top_n_df <- result_df[order(-result_df$joint_log_likelihood), ]
    head(top_n_df, top_n)
  }
  
  evaluate_one_sample <- function(branch_id) {
    total_mutations <- sum(as.numeric(af_data[af_data$branch_name == branch_id, bin_cols, drop = FALSE]))
    
    if (total_mutations == 0) {
      warning(sprintf("⚠ Skipping sample '%s': no mutations.", branch_id))
      return(NULL)
    }
    
    log_liks <- sapply(seq_len(nrow(base_grid)), function(i) {
      params <- base_grid[i, ]
      p_row <- as.numeric(params[bin_cols]) + epsilon
      p_row <- p_row / sum(p_row)
      obs_counts <- as.numeric(af_data[af_data$branch_name == branch_id, bin_cols, drop = FALSE])
      sum(obs_counts * log(p_row))
    })
    
    result_df <- base_grid[, c("rho", "m", "omega")]
    result_df$single_log_likelihood <- log_liks
    result_df$branch_name <- branch_id
    
    top_n_df <- result_df[order(-result_df$single_log_likelihood), ]
    head(top_n_df, top_n)
  }
  
  joint_results_list <- parallel::mclapply(
    batches,
    evaluate_one_batch,
    mc.cores = n_cores
  )
  joint_results_list <- Filter(Negate(is.null), joint_results_list)
  joint_best_models <- do.call(rbind, joint_results_list)
  rownames(joint_best_models) <- NULL
  
  single_results_list <- parallel::mclapply(
    unique(af_data$branch_name),
    evaluate_one_sample,
    mc.cores = n_cores
  )
  single_results_list <- Filter(Negate(is.null), single_results_list)
  single_best_models <- do.call(rbind, single_results_list)
  rownames(single_best_models) <- NULL
  
  if (!return_summary) {
    return(list(
      joint_best_models = joint_best_models,
      single_best_models = single_best_models
    ))
  }
  
  list(
    joint_best_models = joint_best_models,
    single_best_models = single_best_models
  )
}





# ----------------------------
# Final Unified AF Mixture Analysis Script (Fully Fixed)
# ----------------------------

library(parallel)
library(dplyr)

# Requires simulate_process_to_fixed_iteration() from simulate_processes.cpp
# Rcpp::sourceCpp("simulate_processes.cpp")

# Run simulation for one process
run_n_processes <- function(m, rho, n, small_arc_threshold, mu_behavior) {
  simulate_process_to_fixed_iteration(
    rho = rho,
    m = m,
    n_start = 1,
    n = n,
    small_arc_threshold = small_arc_threshold,
    mu_behavior = mu_behavior
  )
}

# Bin helper
count_true_partitions <- function(af_vector, u = 10, lower_bound = 0.01) {
  af_vector <- af_vector[!is.na(af_vector) & af_vector >= lower_bound & af_vector <= 1]
  bin_edges <- seq(lower_bound, 1, length.out = u + 1)
  bin_counts <- hist(af_vector, breaks = bin_edges, plot = FALSE, include.lowest = TRUE)$counts
  range_info <- data.frame(lower = head(bin_edges, -1), upper = tail(bin_edges, -1))
  return(list(counts = bin_counts, range_info = range_info))
}

# Final fixed bin label function
generate_average_mixture_vector <- function(rho, m, omega, time_vec,
                                            n_replicates = 100, G = 240e6,
                                            u = 10, lower_bound = 0.01,
                                            mu_behavior = "random_each_iter",
                                            small_arc_threshold = 1e-11,
                                            digits = 5) {
  n_processes <- length(time_vec)
  lambda_vec <- 2 * G * omega * time_vec
  weights <- lambda_vec / sum(lambda_vec)
  
  mixture_matrix <- matrix(0, nrow = n_replicates, ncol = u)
  
  for (rep in 1:n_replicates) {
    q_matrix <- matrix(0, nrow = u, ncol = n_processes)
    for (k_idx in seq_len(n_processes)) {
      sim_matrix <- run_n_processes(m, rho, n_processes, small_arc_threshold, mu_behavior)
      af_k <- sim_matrix[, k_idx]
      bin_counts <- count_true_partitions(af_k, u = u, lower_bound = lower_bound)$counts
      bin_freqs <- bin_counts / sum(bin_counts)
      q_matrix[, k_idx] <- bin_freqs
    }
    p_vec <- as.vector(q_matrix %*% weights)
    p_vec <- p_vec / sum(p_vec + 1e-12)
    mixture_matrix[rep, ] <- p_vec
  }
  
  avg_mixture <- colMeans(mixture_matrix)
  avg_mixture <- avg_mixture / sum(avg_mixture)
  avg_mixture <- round(avg_mixture, digits = digits)
  
  bin_edges <- seq(lower_bound, 1, length.out = u + 1)
  bin_labels <- paste0("bin_", sprintf("%.3f", head(bin_edges, -1)), "-", sprintf("%.3f", tail(bin_edges, -1)))
  names(avg_mixture) <- bin_labels
  
  return(avg_mixture)
}



# Rebuilds grid
compute_af_grid_across_samples <- function(rho_vals, m_vals, omega_vals,
                                           time_data, u, lower_bound,
                                           G, n_replicates,
                                           mu_behavior = "random_each_iter",
                                           small_arc_threshold = 1e-11,
                                           digits = 5,
                                           n_cores = 1) {
  samples <- time_data$branch_name
  param_grid <- expand.grid(rho = rho_vals, m = m_vals, omega = omega_vals,
                            sample_id = samples, stringsAsFactors = FALSE)
  
  cat("🔄 Starting simulation for", nrow(param_grid), "parameter combinations\n")
  
  # --- Cache to avoid redundant grid simulation ---
  grid_cache <- new.env(hash = TRUE)
  
  grid_results <- parallel::mclapply(seq_len(nrow(param_grid)), function(i) {
    row <- param_grid[i, ]
    time_vec <- as.numeric(time_data[time_data$branch_name == row$sample_id, grep("^t", names(time_data))])
    time_vec <- time_vec[!is.na(time_vec)]
    
    if (length(time_vec) == 0) {
      warning(sprintf("Skipping %s: empty or invalid time vector", row$sample_id))
      return(NULL)
    }
    
    # Build cache key from time_vec and parameter combo
    time_key <- paste0(round(time_vec, 6), collapse = "_")
    cache_key <- paste(time_key, row$rho, row$m, row$omega, sep = "_")
    
    # Use cached grid if available
    if (exists(cache_key, envir = grid_cache)) {
      mixture <- get(cache_key, envir = grid_cache)
    } else {
      mixture <- generate_average_mixture_vector(row$rho, row$m, row$omega, time_vec,
                                                 n_replicates = getOption("samSFS.n_replicates", 100),
                                                 G = G, u = u, lower_bound = lower_bound,
                                                 mu_behavior = mu_behavior,
                                                 small_arc_threshold = small_arc_threshold,
                                                 digits = digits)
      assign(cache_key, mixture, envir = grid_cache)
    }
    
    if (!is.null(mixture)) {
      return(data.frame(row, as.list(mixture), stringsAsFactors = FALSE, check.names = FALSE))
    } else {
      warning(sprintf("Simulation failed for %s", row$sample_id))
      return(NULL)
    }
  }, mc.cores = n_cores)
  
  df <- do.call(rbind, grid_results)
  
  if (is.null(df) || nrow(df) == 0) {
    warning("Grid simulation failed: no rows returned.")
    return(NULL)
  }
  
  return(df)
}


analyze_real_data_wrapper <- function(time_data, af_data, output_dir,
                                      AM_placement = "random_each_iter",
                                      rho_vals, m_vals, omega_vals,
                                      n_replicates = 100, u = 10, lower_bound = 0.01,
                                      G = 240e6, n_cores_samples = 2, n_cores_grid = 4,
                                      n_cores_eval = 4, joint_mode = TRUE, top_n = 100,
                                      save_outputs = TRUE, return_summary = TRUE,
                                      verbose = TRUE,
                                      grid_override = NULL,
                                      total_mutations_override = NULL,
                                      apply_omega_penalty = TRUE) {
  
  start_time <- Sys.time()
  
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    log_file_path <- file.path(output_dir, "run_log.txt")
    log_connection <- file(log_file_path, open = "wt")
    sink(log_connection, type = "output")
    sink(log_connection, type = "message")
    cat(sprintf("📅 Run started at: %s\n", start_time))
    cat(sprintf("Output folder: %s\n\n", normalizePath(output_dir)))
  }
  
  if (verbose) cat("Starting real data analysis...\n")
  if (verbose) cat("Computing theoretical allele frequency grid...\n")
  
  if (!is.null(grid_override)) {
    grid_output <- grid_override
  } else {
    grid_output <- compute_af_grid_across_samples(...)
  }
  
  if (!is.null(output_dir) && save_outputs) {
    saveRDS(grid_output, file = file.path(output_dir, "grid_output.rds"))
    cat("✔️ Saved grid_output.rds\n")
  }
  
  if (verbose) cat("Binning real allele frequency data...\n")
  bin_edges <- seq(lower_bound, 1, length.out = u + 1)
  bin_labels <- paste0("bin_", sprintf("%.3f", head(bin_edges, -1)), "-", sprintf("%.3f", tail(bin_edges, -1)))
  
  binned_list <- lapply(split(af_data, af_data$branch_name), function(df_branch) {
    af_vec <- df_branch$allele_frequency
    result <- count_true_partitions(af_vec, u = u, lower_bound = lower_bound)
    full_bin_vec <- setNames(rep(0, u), bin_labels)
    names(result$counts) <- paste0("bin_", sprintf("%.3f", result$range_info$lower), "-", sprintf("%.3f", result$range_info$upper))
    full_bin_vec[names(result$counts)] <- result$counts
    out <- as.data.frame(t(full_bin_vec))
    out$branch_name <- unique(df_branch$branch_name)
    out
  })
  
  af_data_binned <- do.call(rbind, binned_list)
  af_data_binned <- af_data_binned[, c("branch_name", bin_labels)]
  rownames(af_data_binned) <- NULL
  
  if ("batch_id" %in% names(time_data)) {
    if (verbose) cat("Merging batch_id information from time_data...\n")
    af_data_binned <- merge(af_data_binned, time_data[, c("branch_name", "batch_id")],
                            by = "branch_name", all.x = TRUE)
  } else {
    warning("No 'batch_id' column found in time_data. Proceeding without batch grouping.")
  }
  
  if (verbose) cat("Evaluating likelihoods for real AF data...\n")
  lik_output <- safe_likelihood_evaluation_joint(grid_output, af_data_binned,
                                                 joint_mode, top_n, return_summary,
                                                 n_cores = n_cores_eval)
  
  # 🔧 L2 Loss on Mutation Count with dynamic scaling
  if (apply_omega_penalty && !is.null(total_mutations_override) && total_mutations_override > 0) {
    observed <- total_mutations_override
    correction <- 1 / (1 - lower_bound)
    
    # Use all time intervals from all samples
    time_vec_all <- as.numeric(unlist(time_data[, grep("^t", names(time_data))]))
    time_vec <- time_vec_all[!is.na(time_vec_all)]
    
    # Define scaling parameter per mutation
    penalty_strength_per_mutation <- 1.0  # adjust this if needed
    
    if (!is.null(lik_output$single_best_models)) {
      lambda_single <- 2*G * lik_output$single_best_models$omega * sum(time_vec) * correction
      l2_loss_single <- (observed - lambda_single)^2
      penalty_weight_single <- penalty_strength_per_mutation * observed
      
      lik_output$single_best_models$single_log_likelihood <-
        lik_output$single_best_models$single_log_likelihood - penalty_weight_single * l2_loss_single
    }
    
    if (!is.null(lik_output$joint_best_models)) {
      lambda_joint <- 2*G * lik_output$joint_best_models$omega * sum(time_vec) * correction
      l2_loss_joint <- (observed - lambda_joint)^2
      penalty_weight_joint <- penalty_strength_per_mutation * observed
      
      # 🔎 Diagnostics
      cat("\n🔎 L2 Penalty Diagnostics (Dynamic Scaling):\n")
      #cat("  Observed mutations:", observed, "\n")
      #cat("  Omega (first):", lik_output$joint_best_models$omega[1], "\n")
      #cat("  Expected lambda (at omega):", lambda_joint[1], "\n")
      #cat("  L2 penalty weight:", penalty_weight_joint, "\n")
      #cat("  Total penalty:", penalty_weight_joint * l2_loss_joint[1], "\n")
      
      lik_output$joint_best_models$joint_log_likelihood <-
        lik_output$joint_best_models$joint_log_likelihood - penalty_weight_joint * l2_loss_joint
    }
    
  } else if (isTRUE(apply_omega_penalty)) {
    warning("No L2 count loss added (zero or missing total_mutations_override).")
  }
  
  if (!is.null(output_dir) && save_outputs) {
    saveRDS(lik_output, file = file.path(output_dir, "likelihood_results.rds"))
    cat("✔️ Saved likelihood_results.rds\n")
  }
  
  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "mins")
  cat(sprintf("\n Run finished at: %s\n", end_time))
  cat(sprintf("Total runtime: %.2f minutes\n", as.numeric(runtime)))
  
  if (!is.null(output_dir)) {
    sink(type = "message")
    sink(type = "output")
    close(log_connection)
  }
  
  list(grid = grid_output, likelihood = lik_output, af_data_binned = af_data_binned)
}





run_estimation_and_summarize <- function(time_data,
                                         af_data,
                                         true_params,
                                         grid_cache,
                                         total_mutations = NULL,
                                         joint_mode = TRUE,
                                         output_dir = "simulation_output",
                                         ...) {
  if (nrow(af_data) == 0) {
    stop("Simulation aborted: no allele frequencies were generated.")
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  time_vec <- as.numeric(time_data[1, grep("^t", names(time_data))])
  grid_key <- paste0(round(time_vec, 6), collapse = "_")
  
  if (!grid_key %in% names(grid_cache)) {
    stop("Grid cache does not contain entry for time vector: ", grid_key)
  }
  
  grid_output <- grid_cache[[grid_key]]
  
  res <- analyze_real_data_wrapper(
    time_data = time_data,
    af_data = af_data,
    grid_override = grid_output,
    output_dir = NULL,
    joint_mode = joint_mode,
    top_n = getOption("samSFS.grid_top_n", 25),
    save_outputs = FALSE,
    return_summary = TRUE,
    verbose = FALSE,
    n_cores_grid = getOption("samSFS.grid_cores", 2),
    n_cores_eval = getOption("samSFS.eval_cores", 2),
    rho_vals = getOption("samSFS.grid_rho"),
    m_vals = getOption("samSFS.grid_m"),
    omega_vals = getOption("samSFS.grid_omega"),
    u = getOption("samSFS.grid_u", 10),
    lower_bound = getOption("samSFS.grid_lower_bound", 0.01),
    G = getOption("samSFS.grid_G", 240e6),
    total_mutations_override = total_mutations,
    ...
  )
  
  if (is.null(res$likelihood$single_best_models)) stop("Estimation failed: single_best_models is NULL")
  if (is.null(res$likelihood$joint_best_models)) stop("Estimation failed: joint_best_models is NULL")
  
  ind_df <- res$likelihood$single_best_models
  ind_df$log_likelihood <- ind_df$single_log_likelihood
  individual_summary <- summarize_top_models_weighted(
    df = ind_df,
    id_col = "branch_name",
    likelihood_col = "log_likelihood",
    top_n = getOption("samSFS.grid_top_n", 25)
  )
  
  joint_df <- res$likelihood$joint_best_models
  joint_df$log_likelihood <- joint_df$joint_log_likelihood
  batched_summary <- summarize_top_models_weighted(
    df = joint_df,
    id_col = "batch_id",
    likelihood_col = "log_likelihood",
    top_n = getOption("samSFS.grid_top_n", 25)
  )
  
  for (df in list(individual_summary, batched_summary)) {
    df$true_rho <- true_params$rho
    df$true_m <- true_params$m
    df$true_omega <- true_params$omega
  }
  
  return(list(individual = individual_summary, batched = batched_summary))
}




precompute_grid_cache <- function(time_data, rho_vals, m_vals, omega_vals,
                                  u = 10, G = 240e6, lower_bound = 0.01,
                                  n_replicates = 10, n_cores = 1) {
  time_vec_list <- extract_unique_timevecs(time_data)
  grid_cache <- list()
  
  for (tv_row in time_vec_list) {
    tv <- as.numeric(tv_row)
    key <- paste0(round(tv, 6), collapse = "_")
    
    message("⏳ Precomputing grid for time_vec: ", key)
    
    dummy_time_data <- data.frame(
      branch_name = "template_sample",
      t1 = tv[1], t2 = tv[2], t3 = tv[3]
    )
    
    grid_cache[[key]] <- compute_af_grid_across_samples(
      rho_vals = rho_vals,
      m_vals = m_vals,
      omega_vals = omega_vals,
      time_data = dummy_time_data,
      u = u,
      lower_bound = lower_bound,
      G = G,
      n_replicates = n_replicates,
      n_cores = n_cores
    )
  }
  
  return(grid_cache)
}




analyze_real_data_fresh <- function(time_data,
                                    af_data,
                                    output_dir,
                                    G,
                                    n_cores_grid,
                                    n_cores_eval,
                                    rho_vals = seq(0.05, 0.95, by = 0.05),
                                    m_vals = seq(2, 20, by = 1),
                                    omega_vals = 10^seq(-11, -7, length.out = 200),
                                    n_replicates = 100,
                                    u = 5,
                                    lower_bound = 0.01,
                                    top_n = 30,
                                    save_outputs = TRUE,
                                    return_summary = TRUE,
                                    verbose = TRUE,
                                    apply_omega_penalty = TRUE,
                                    total_mutations_override = sum(af_data$allele_frequency > lower_bound),
                                    AM_placement = "random_each_iter") {
  
  start_time <- Sys.time()
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    log_file_path <- file.path(output_dir, "run_log.txt")
    log_connection <- file(log_file_path, open = "wt")
    sink(log_connection, type = "output")
    sink(log_connection, type = "message")
    cat(sprintf("Run started at: %s\n", start_time))
    cat(sprintf("Output folder: %s\n\n", normalizePath(output_dir)))
  }
  
  if (verbose) cat("Starting real data analysis...\n")
  if (verbose) cat("Computing theoretical allele frequency grid...\n")
  
  grid_output <- compute_af_grid_across_samples(
    rho_vals = rho_vals,
    m_vals = m_vals,
    omega_vals = omega_vals,
    time_data = time_data,
    u = u,
    lower_bound = lower_bound,
    G = G,
    n_replicates = n_replicates,
    mu_behavior = AM_placement,
    small_arc_threshold = 1e-11,
    digits = 5,
    n_cores = n_cores_grid
  )
  
  if (!is.null(output_dir) && save_outputs) {
    saveRDS(grid_output, file = file.path(output_dir, "grid_output.rds"))
    cat("✔️ Saved grid_output.rds\n")
  }
  
  if (verbose) cat("Binning real allele frequency data...\n")
  bin_edges <- seq(lower_bound, 1, length.out = u + 1)
  bin_labels <- paste0("bin_", sprintf("%.3f", head(bin_edges, -1)), "-", sprintf("%.3f", tail(bin_edges, -1)))
  
  binned_list <- lapply(split(af_data, af_data$branch_name), function(df_branch) {
    af_vec <- df_branch$allele_frequency
    result <- count_true_partitions(af_vec, u = u, lower_bound = lower_bound)
    full_bin_vec <- setNames(rep(0, u), bin_labels)
    names(result$counts) <- paste0("bin_", sprintf("%.3f", result$range_info$lower), "-", sprintf("%.3f", result$range_info$upper))
    full_bin_vec[names(result$counts)] <- result$counts
    out <- as.data.frame(t(full_bin_vec))
    out$branch_name <- unique(df_branch$branch_name)
    out
  })
  
  af_data_binned <- do.call(rbind, binned_list)
  af_data_binned <- af_data_binned[, c("branch_name", bin_labels)]
  rownames(af_data_binned) <- NULL
  
  if ("batch_id" %in% names(time_data)) {
    if (verbose) cat("Merging batch_id information from time_data...\n")
    af_data_binned <- merge(af_data_binned, time_data[, c("branch_name", "batch_id")], by = "branch_name", all.x = TRUE)
  } else {
    warning("No 'batch_id' column found in time_data. Proceeding without batch grouping.")
  }
  
  if (verbose) cat("Evaluating likelihoods for real AF data...\n")
  lik_output <- safe_likelihood_evaluation_joint(grid_output, af_data_binned, TRUE, top_n, return_summary, n_cores = n_cores_eval)
  
  if (apply_omega_penalty && !is.null(total_mutations_override) && total_mutations_override > 0) {
    observed <- total_mutations_override
    correction <- 1 / (1 - lower_bound)
    time_vec_all <- as.numeric(unlist(time_data[, grep("^t", names(time_data))]))
    time_vec <- time_vec_all[!is.na(time_vec_all)]
    penalty_strength_per_mutation <- 1.0
    
    if (!is.null(lik_output$single_best_models)) {
      lambda_single <- 2*G * lik_output$single_best_models$omega * sum(time_vec) * correction
      l2_loss_single <- (observed - lambda_single)^2
      penalty_weight_single <- penalty_strength_per_mutation * observed
      lik_output$single_best_models$single_log_likelihood <-
        lik_output$single_best_models$single_log_likelihood - penalty_weight_single * l2_loss_single
    }
    
    if (!is.null(lik_output$joint_best_models)) {
      lambda_joint <- 2*G * lik_output$joint_best_models$omega * sum(time_vec) * correction
      l2_loss_joint <- (observed - lambda_joint)^2
      penalty_weight_joint <- penalty_strength_per_mutation * observed
      lik_output$joint_best_models$joint_log_likelihood <-
        lik_output$joint_best_models$joint_log_likelihood - penalty_weight_joint * l2_loss_joint
    }
  }
  
  summary_individual <- summarize_top_models_weighted(
    df = lik_output$single_best_models,
    id_col = "branch_name",
    likelihood_col = "single_log_likelihood",
    top_n = top_n
  )
  
  summary_joint <- summarize_top_models_weighted(
    df = lik_output$joint_best_models,
    id_col = "batch_id",
    likelihood_col = "joint_log_likelihood",
    top_n = top_n
  )
  
  if (!is.null(output_dir) && save_outputs) {
    saveRDS(lik_output, file = file.path(output_dir, "likelihood_results.rds"))
    write.csv(summary_individual, file = file.path(output_dir, "summary_individual.csv"), row.names = FALSE)
    write.csv(summary_joint, file = file.path(output_dir, "summary_joint.csv"), row.names = FALSE)
    cat("✔️ Saved likelihood_results.rds, summary_individual.csv, summary_joint.csv\n")
  }
  
  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "mins")
  cat(sprintf("\n Run finished at: %s\n", end_time))
  cat(sprintf("Total runtime: %.2f minutes\n", as.numeric(runtime)))
  
  if (!is.null(output_dir)) {
    sink(type = "message")
    sink(type = "output")
    close(log_connection)
  }
  
  list(grid = grid_output,
       likelihood = lik_output,
       af_data_binned = af_data_binned,
       summary_individual = summary_individual,
       summary_joint = summary_joint)
}




# samSFS_setup.R

set_samSFS_threads <- function(simulation_cores = 2, grid_cores = 2, eval_cores = 2) {
  RcppParallel::setThreadOptions(numThreads = simulation_cores)
  options(samSFS.grid_cores = grid_cores)
  options(samSFS.eval_cores = eval_cores)
}

set_samSFS_grid <- function(rho_vals, m_vals, omega_vals,
                            u = 5, G = 240e6, lower_bound = 0.01,
                            top_n = 30, n_replicates = 100) {
  options(samSFS.grid_rho = rho_vals)
  options(samSFS.grid_m = m_vals)
  options(samSFS.grid_omega = omega_vals)
  options(samSFS.grid_u = u)
  options(samSFS.grid_G = G)
  options(samSFS.grid_lower_bound = lower_bound)
  options(samSFS.grid_top_n = top_n)
  options(samSFS.n_replicates = n_replicates)
}




# samSFS_simulation_core.R

# This file contains the patched simulation core logic including mutation tracking and Poisson support

# --- simulate_af_data() and related logic in patched form already delivered in samSFS_poisson_patch

# simulate_af_data: returns af_data and total_mutations
# run_samSFS_simulation: coordinates one simulation run
# run_estimation_and_summarize: connects AF data to inference engine
# These functions are already adapted and included in samSFS_poisson_patch


simulate_af_data <- function(time_data,
                             rho, m, omega,
                             G = getOption("samSFS.grid_G", 240e6),
                             mu_behavior = "random_each_iter",
                             small_arc_threshold = 1e-11,
                             lower_bound = getOption("samSFS.grid_lower_bound", 0.01)) {
  
  n_samples <- nrow(time_data)
  long_af_rows <- list()
  total_mutations <- 0
  correction <- 1 / (1 - lower_bound)
  
  for (i in seq_len(n_samples)) {
    branch_name <- time_data$branch_name[i]
    batch_id <- time_data$batch_id[i]
    time_vec <- as.numeric(time_data[i, grep("^t", names(time_data))])
    n_intervals <- length(time_vec)
    
    sim_matrix <- simulate_process_to_fixed_iteration(
      rho = rho,
      m = m,
      n_start = 1,
      n = n_intervals,
      small_arc_threshold = small_arc_threshold,
      mu_behavior = mu_behavior
    )
    
    mutation_counts <- rpois(length(time_vec), lambda = 2*G * omega * time_vec * correction)
    total_mutations <- total_mutations + sum(mutation_counts)
    
    for (k in seq_len(n_intervals)) {
      af_k <- sim_matrix[, k]
      if (mutation_counts[k] > 0) {
        sampled_afs <- sample(af_k, size = mutation_counts[k], replace = TRUE)
        long_af_rows[[length(long_af_rows) + 1]] <- data.frame(
          branch_name = branch_name,
          allele_frequency = sampled_afs,
          batch_id = batch_id
        )
      }
    }
  }
  
  if (total_mutations == 0) {
    warning("No mutations were generated in this simulation. Using failsafe value of 1.")
    total_mutations <- 1
  }
  
  af_df <- do.call(rbind, long_af_rows)
  return(list(af_data = af_df, total_mutations = total_mutations))
}





# PATCHED run_samSFS_simulation() with progress messages
run_samSFS_simulation <- function(n_samples, branching_times, samples_per_batch,
                                  true_params, grid_cache,
                                  output_dir = "simulation_output") {
  
  if (is.list(true_params) && !is.data.frame(true_params)) {
    true_params <- as.data.frame(true_params)
  }
  
  if (!("rho" %in% names(true_params) && "m" %in% names(true_params) && "omega" %in% names(true_params))) {
    stop("true_params must contain columns: rho, m, omega")
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  results_list <- list()
  
  for (i in seq_len(nrow(true_params))) {
    param_row <- true_params[i, ]
    param_suffix <- sprintf("rho%.3f_m%d_omega%.2e", param_row$rho, param_row$m, param_row$omega)
    message(sprintf("Running simulation %d of %d: %s", i, nrow(true_params), param_suffix))
    
    param_output_dir <- file.path(output_dir, param_suffix)
    dir.create(param_output_dir, showWarnings = FALSE, recursive = TRUE)
    
    time_data <- generate_time_data(
      n_samples = n_samples,
      branching_times = branching_times,
      samples_per_batch = samples_per_batch
    )
    
    af_result <- simulate_af_data(
      time_data = time_data,
      rho = param_row$rho,
      m = param_row$m,
      omega = param_row$omega
    )
    
    af_data <- af_result$af_data
    total_mutations <- af_result$total_mutations
    
    est_result <- run_estimation_and_summarize(
      time_data = time_data,
      af_data = af_data,
      true_params = list(
        rho = param_row$rho,
        m = param_row$m,
        omega = param_row$omega
      ),
      grid_cache = grid_cache,
      total_mutations = total_mutations,
      output_dir = param_output_dir
    )
    
    # Add parameter reference columns
    est_result$individual$rho <- param_row$rho
    est_result$individual$m <- param_row$m
    est_result$individual$omega <- param_row$omega
    
    est_result$batched$rho <- param_row$rho
    est_result$batched$m <- param_row$m
    est_result$batched$omega <- param_row$omega
    
    # Write outputs
    write.csv(est_result$individual,
              file = file.path(output_dir, paste0("individual_", param_suffix, ".csv")), row.names = FALSE)
    
    write.csv(est_result$batched,
              file = file.path(output_dir, paste0("batched_", param_suffix, ".csv")), row.names = FALSE)
    
    # Remove subfolder if it's empty
    if (length(list.files(param_output_dir, recursive = TRUE)) == 0) {
      unlink(param_output_dir, recursive = TRUE)
      message("Removed empty folder: ", param_output_dir)
    }
    
    results_list[[i]] <- list(individual = est_result$individual, batched = est_result$batched)
    
  }
  
  # If only one combination, return like before
  if (nrow(true_params) == 1) {
    return(results_list[[1]])
  }
  
  # Otherwise return all
  individual_all <- do.call(rbind, lapply(results_list, `[[`, "individual"))
  batched_all <- do.call(rbind, lapply(results_list, `[[`, "batched"))
  
  return(list(individual = individual_all, batched = batched_all))
}



# samSFS_time.R

generate_time_data <- function(n_samples = 10,
                               branching_times = c(1.0, 1.0, 1.0),
                               samples_per_batch = 1) {
  n_intervals <- length(branching_times)
  time_df <- as.data.frame(matrix(rep(branching_times, each = n_samples), nrow = n_samples, byrow = TRUE))
  names(time_df) <- paste0("t", seq_len(n_intervals))
  time_df$branch_name <- paste0("sample_", seq_len(n_samples))
  
  batch_ids <- rep(paste0("batch_", seq_len(ceiling(n_samples / samples_per_batch))),
                   each = samples_per_batch)[seq_len(n_samples)]
  time_df$batch_id <- batch_ids
  
  time_df <- time_df[, c("branch_name", paste0("t", seq_len(n_intervals)), "batch_id")]
  return(time_df)
}

extract_unique_timevecs <- function(time_data) {
  time_cols <- grep("^t", names(time_data), value = TRUE)
  unique_rows <- unique(time_data[, time_cols, drop = FALSE])
  split(unique_rows, seq_len(nrow(unique_rows)))
}



summarize_top_models_weighted <- function(
    df,
    id_col = "sample_id",
    likelihood_col = "log_likelihood",
    top_n = 10,
    round_output = FALSE,
    digits = 4,
    output_dir = NULL
) {
  library(dplyr)
  
  df <- df %>%
    mutate(across(c(rho, m, omega), as.numeric))
  
  # Check that all groups have enough data
  too_few <- df %>%
    group_by(.data[[id_col]]) %>%
    summarise(n_rows = n(), .groups = "drop") %>%
    filter(n_rows < top_n)
  
  if (nrow(too_few) > 0) {
    stop("Error: The following ", id_col, "(s) have fewer than ", top_n, " rows:\n",
         paste(too_few[[id_col]], collapse = ", "))
  }
  
  # Weighted median helper
  weighted_median <- function(values, weights) {
    ord <- order(values)
    values <- values[ord]
    weights <- weights[ord]
    cum_weights <- cumsum(weights)
    values[which(cum_weights >= 0.5)[1]]
  }
  
  result <- df %>%
    group_by(.data[[id_col]]) %>%
    slice_max(order_by = .data[[likelihood_col]], n = top_n) %>%
    group_modify(~{ 
      logliks <- .x[[likelihood_col]]
      weights <- exp(logliks - max(logliks))  # stabilize
      weights <- weights / sum(weights)
      
      wrho   <- sum(.x$rho   * weights)
      wm     <- sum(.x$m     * weights)
      womega <- sum(.x$omega * weights)
      avg_ll <- sum(logliks * weights)
      
      data.frame(
        rho_mean   = wrho,
        rho_median = weighted_median(.x$rho, weights),
        rho_sd     = sqrt(sum((.x$rho - wrho)^2 * weights)),
        
        m_mean     = wm,
        m_median   = weighted_median(.x$m, weights),
        m_sd       = sqrt(sum((.x$m - wm)^2 * weights)),
        
        omega_mean   = womega,
        omega_median = weighted_median(.x$omega, weights),
        omega_sd     = sqrt(sum((.x$omega - womega)^2 * weights)),
        
        avg_log_likelihood = avg_ll
      )
    }) %>%
    ungroup()
  
  if (round_output) {
    result <- result %>%
      mutate(across(where(is.numeric), ~ signif(.x, digits = digits)))
  }
  
  result <- as.data.frame(result)
  
  # Write to CSV if requested
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    out_name <- paste0("summary_", id_col, ".csv")
    out_path <- file.path(output_dir, out_name)
    write.csv(result, file = out_path, row.names = FALSE)
    message("Summary written to: ", out_path)
  }
  
  return(result)
}

