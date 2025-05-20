# === Real Data Analysis (Run interactively or via Rscript) ===

# --- Clean up workspace ---
rm(list = ls()); gc()

# --- OPTIONAL: Manually set working directory if not in samSFS root ---
# setwd("path/to/samSFS")

# --- 1. Load required packages ---
required_pkgs <- c("dplyr", "tidyr", "circular", "Rcpp", "RcppArmadillo", "parallel", "RcppParallel")
new_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(new_pkgs)) install.packages(new_pkgs, repos = "https://cloud.r-project.org", quiet = TRUE)
suppressPackageStartupMessages({
  lapply(required_pkgs, library, character.only = TRUE)
})

# --- 2. Load samSFS source code ---
source("R/samSFS_functions.R")
Rcpp::sourceCpp("src/simulate_processes.cpp")

# --- 3. Configure threading ---
set_samSFS_threads(simulation_cores = 2, grid_cores = 4, eval_cores = 2)

# --- 4. Load input data ---
time_data <- read.csv("data/time_data.csv", stringsAsFactors = FALSE)
af_data   <- read.csv("data/af_data.csv", stringsAsFactors = FALSE)


# --- 5. Run real data analysis (using coarse grid here for short run time)
analyze_real_data_fresh(
  time_data     = time_data,
  af_data       = af_data,
  output_dir    = "results/real_data_results",
  G             = 240e6,
  n_cores_grid  = 4,
  n_cores_eval  = 2,
  rho_vals      = seq(0.4, 0.7, by = 0.1),
  m_vals        = seq(3, 10, by = 2),
  omega_vals    = 10^seq(-10, -8, length.out = 5),
  top_n         = 20
)

# --- 6. Notify user ---
cat("✅ Real data analysis complete. Results saved to: results/real_data_results\n")
