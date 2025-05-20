# === Real Data Analysis Runner (with integrated wrapper) ===
rm(list = ls()); gc()
setwd("~/Dropbox/jlab/Projects/Fixation_of_somatic_mutations_in_trees/FUNCTIONS_SIMULATION-2")

# --- Load packages ---
required_pkgs <- c("dplyr", "tidyr", "circular", "Rcpp", "RcppArmadillo", "parallel", "RcppParallel")
new_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(new_pkgs)) install.packages(new_pkgs, repos = "https://cloud.r-project.org", quiet = TRUE)
suppressPackageStartupMessages({
  lapply(required_pkgs, library, character.only = TRUE)
})

# --- Load sources ---
source("samSFS_setup.R")
source("samSFS_time.R")
source("samSFS_grid.R")
source("samSFS_simulation_core.R")
source("samSFS_analysis.R")
source("summarize_top_models_weighted.R")
source("samSFS_allinone_FIXED.R")
source("safe_likelihood_evaluation_joint_full.R")
Rcpp::sourceCpp("simulate_processes.cpp")

# --- Threading configuration ---
set_samSFS_threads(simulation_cores = 2, grid_cores = 4, eval_cores = 2)

# --- Load real data ---
time_data <- read.csv("time_data.csv", stringsAsFactors = FALSE)
af_data <- read.csv("af_data.csv", stringsAsFactors = FALSE)

# --- Output directory ---
output_dir <- "real_data_results"

# --- Integrated real-data-specific wrapper function ---
source("samSFS_realdata_wrapper.R")  # Load analyze_real_data_fresh()

# --- Run real data analysis ---
analyze_real_data_fresh(
  time_data = time_data,
  af_data = af_data,
  output_dir = "real_data_results",
  G = 240e6,
  n_cores_grid = getOption("samSFS.grid_cores", 2),
  n_cores_eval = getOption("samSFS.eval_cores", 2),
)


cat("✅ Real data analysis complete. Results saved to:", normalizePath(output_dir), "\n")
