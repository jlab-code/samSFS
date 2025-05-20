# samSFS


Tools to simulate and analyze allele frequency data from somatic evolution in tree-like structures.

## 🔧 Installation

Clone the repo and load in R:

```r
# Load all core functions
files <- list.files("R", full.names = TRUE)
sapply(files, source)
Rcpp::sourceCpp("src/simulate_processes.cpp")
```

## 🚀 Simulation Example

```r
samSFS_set_threads()
samSFS_set_grid(...)
samSFS_simulate(...)
```

## 🧬 Real Data Example

```r
samSFS_analyze_real(
  time_data = read.csv("example/time_data.csv"),
  af_data = read.csv("example/af_data.csv"),
  output_dir = "results/",
  G = 240e6,
  n_cores_grid = 4,
  n_cores_eval = 2
)
=======
This repository provides a clean, modular implementation of the samSFS simulation and inference pipeline using the original function names.

## Structure

- `R/`: All core functions (one per file)
- `run/`: Example/test scripts
- `data/`: Sample input data
- `src/`: C++ simulation engine
- `legacy/`: Unmodified original scripts

## Usage

```r
source("R/<function_name>.R")
source("run/test_simulated_data.R")
source("run/test_real_data.R")
>>>>>>> 15912e7 (Update all source files: unified R functions and simulation core)
```
