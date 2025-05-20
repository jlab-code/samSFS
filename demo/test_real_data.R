# === Real Data Analysis Runner ===

rm(list = ls()); gc()

# --- Optional: manually set working directory if needed ---
# setwd("path/to/samSFS")

# --- Load required packages ---
required_pkgs <- c("dplyr", "tidyr", "circular", "Rcpp", "RcppArmadillo", "parallel", "RcppParallel")
new_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(new_pkgs)) install.packages(new_pkgs, repos = "https://cloud.r-project.org", quiet = TRUE)
suppressPackageStartupMessages({
  lapply(required_pkgs, library, character.only = TRUE)
})

# --- Load source files ---
source("R/samSFS_functions.R")
Rcpp::sourceCpp("src/simulate_processes.cpp")

# --- Threading configuration ---
set_samSFS_threads(simulation_cores = 2, grid_cores = 4, eval_cores = 2)

# --- Load real data ---
time_data <- read.csv("data/time_data.csv", stringsAsFactors = FALSE)
af_data <- read.csv("data/af_data.csv", stringsAsFactors = FALSE)

# --- Run real data analysis ---
analyze_real_data_fresh(
  time_data = time_data,
  af_data = af_data,
  output_dir = "results/real_data_results",
  G = 240e6,
  n_cores_grid = getOption("samSFS.grid_cores", 2),
  n_cores_eval = getOption("samSFS.eval_cores", 2)
)

cat("✅ Real data analysis complete. Results saved to: results/real_data_results\n")
