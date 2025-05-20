# === Simulation Demo: Parameter Inference from Simulated Data ===

# --- Clean up workspace ---
rm(list = ls()); gc()

# --- OPTIONAL: Set working directory manually if needed ---
# setwd("path/to/samSFS")

# --- 1. Load required packages ---
required_pkgs <- c("dplyr", "tidyr", "circular", "Rcpp", "RcppArmadillo", "parallel", "RcppParallel")
new_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(new_pkgs)) install.packages(new_pkgs, repos = "https://cloud.r-project.org", quiet = TRUE)
suppressPackageStartupMessages({
  lapply(required_pkgs, library, character.only = TRUE)
})

# --- 2. Load samSFS source files ---
source("R/samSFS_functions.R")
Rcpp::sourceCpp("src/simulate_processes.cpp")

# --- 3. Set threading configuration ---
set_samSFS_threads(simulation_cores = 2, grid_cores = 4, eval_cores = 1)

# --- 4. Define simulation grid parameters ---
set_samSFS_grid(
  rho_vals     = seq(0.35, 0.85, by = 0.05),
  m_vals       = seq(2, 15, by = 1),
  omega_vals   = 10^seq(-11, -7, length.out = 10),
  u            = 5,
  G            = 240e6,
  lower_bound  = 0.01,
  top_n        = 30
)

# --- 5. Simulate time data ---
n_samples         <- 100
branching_times   <- c(10, 10, 10)
samples_per_batch <- 20

time_data <- generate_time_data(
  n_samples         = n_samples,
  branching_times   = branching_times,
  samples_per_batch = samples_per_batch
)

# --- 6. Precompute grid cache ---
grid_cache <- precompute_grid_cache(
  time_data     = time_data,
  rho_vals      = getOption("samSFS.grid_rho"),
  m_vals        = getOption("samSFS.grid_m"),
  omega_vals    = getOption("samSFS.grid_omega"),
  u             = getOption("samSFS.grid_u"),
  G             = getOption("samSFS.grid_G"),
  lower_bound   = getOption("samSFS.grid_lower_bound"),
  n_replicates  = getOption("samSFS.n_replicates"),
  n_cores       = getOption("samSFS.grid_cores")
)

# --- 7. Define true parameters ---
true_params <- expand.grid(
  rho   = c(0.5),
  m     = c(3),
  omega = c(1e-8)
)

# --- 8. Run simulation and analysis ---
result <- run_samSFS_simulation(
  n_samples         = n_samples,
  branching_times   = branching_times,
  samples_per_batch = samples_per_batch,
  true_params       = true_params,
  grid_cache        = grid_cache,
  output_dir        = "results/test_simulation_output"
)

# --- 9. Notify user ---
cat("✅ Simulation complete. Results saved to: results/test_simulation_output\n")
