

# --- Set simulation options ---
#  | Setting            | Affects Stage                        | Tool Used            |
#  | ------------------ | ------------------------------------ | -------------------- |
#  | `simulation_cores` | Simulating AF data (C++ backend)     | `RcppParallel`       |
#  | `grid_cores`       | Grid generation over parameter space | `parallel::mclapply` |
#  | `eval_cores`       | Likelihood evaluation per sample     | `parallel::mclapply` |
#
# This code uses nested parallelism:
# 
# ┌────────────────────────────────────────────┐
# │ 1. R-level Parallelism (grid-level)        │
# │    - Controlled by: n_cores_grid           │
# │    - Used in: compute_af_grid_across_samples() │
# │    - Runs one forked R process per sample/grid row │
# └────────────────────────────────────────────┘
#
# ┌────────────────────────────────────────────┐
# │ 2. C++-level Parallelism (per simulation)  │
# │    - Controlled by: simulation_cores       │
# │    - Used inside each forked process       │
# │    - Applies to simulate_process_to_fixed_iteration() via RcppParallel │
# └────────────────────────────────────────────┘
#
# Total CPU usage = n_cores_grid × simulation_cores
#
# Tip: Choose values so that total thread usage ≈ your number of logical CPU cores.
# Example for a 16-core machine:
#   - simulation_cores = 2
#   - grid_cores       = 8
#



# === Set working directory and load packages ===
rm(list = ls()); gc()
setwd("~/Dropbox/jlab/Projects/Fixation_of_somatic_mutations_in_trees/FUNCTIONS_SIMULATION-2")

# === 0. Install & load required packages ===
required_pkgs <- c("dplyr", "tidyr", "circular", "Rcpp", "RcppArmadillo", "parallel", "RcppParallel")
new_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(new_pkgs)) install.packages(new_pkgs, repos = "https://cloud.r-project.org", quiet = TRUE)
suppressPackageStartupMessages({
  lapply(required_pkgs, library, character.only = TRUE)
})

# --- Load your source files ---
source("samSFS_setup.R")
source("samSFS_time.R")
source("samSFS_grid.R")
source("samSFS_simulation_core.R")
source("samSFS_analysis.R")
source("summarize_top_models_weighted.R")
source("samSFS_allinone_FIXED.R")
source("safe_likelihood_evaluation_joint_full.R")
Rcpp::sourceCpp("simulate_processes.cpp")

# Set these before running:
set_samSFS_threads(simulation_cores = 2, grid_cores = 4, eval_cores = 1)

# --- Setting grid ---
set_samSFS_grid(
  rho_vals = seq(0.35, 0.85, by=0.05),
  m_vals = seq(2,15, by=1),
  omega_vals = 10^seq(-11, -7, length.out = 20),
  u = 5, G = 240e6, lower_bound = 0.01, top_n = 30
)

# --- Generate data and grid ---
n_samples <- 100
branching_times <- c(10, 10, 10)
samples_per_batch <- 20

time_data <- generate_time_data(n_samples, branching_times, samples_per_batch)

grid_cache <- precompute_grid_cache(
  time_data = time_data,
  rho_vals = getOption("samSFS.grid_rho"),
  m_vals = getOption("samSFS.grid_m"),
  omega_vals = getOption("samSFS.grid_omega"),
  u = getOption("samSFS.grid_u"),
  G = getOption("samSFS.grid_G"),
  lower_bound = getOption("samSFS.grid_lower_bound"),
  n_replicates = getOption("samSFS.n_replicates"),
  n_cores = getOption("samSFS.grid_cores")
)

# --- Run simulation for a grid of true parameter values ---
true_params <- expand.grid(
  rho = c(0.5),
  m = c(3),
  omega = c(1e-10)
)

result <- run_samSFS_simulation(
  n_samples = n_samples,
  branching_times = branching_times,
  samples_per_batch = samples_per_batch,
  true_params = true_params,
  grid_cache = grid_cache,
  output_dir = "test_output"
)

