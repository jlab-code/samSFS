# samSFS

**Somatic evolution of stem cell mutations in plants**

This repository provides a clean, modular implementation of the `samSFS` simulation and inference pipeline.

---

## 🔧 Installation

To get started, clone the repository:

```bash
git clone https://github.com/jlab-code/samSFS.git
cd samSFS
```

---

### 🟦 Recommended: Use as an RStudio Project

This project is structured to work cleanly as an RStudio Project:

1. Open RStudio.
2. Go to **File → Open Project...**
3. Select the `samSFS/` folder.
4. This ensures your working directory is set correctly and scripts will run without modification.

You can also create a `.Rproj` file in the repo for convenience.

---

### 🖥️ Alternative: Manual R Console Use

If you're not using RStudio, load the code manually in R:

```r
# Load all core functions
files <- list.files("R", full.names = TRUE)
sapply(files, source)

# Load C++ backend
Rcpp::sourceCpp("src/simulate_processes.cpp")
```

---

## 📁 Project Structure

```
samSFS/
├── R/                  # Core R functions (contains all modular R core functions)
├── src/                # C++ simulation engine
├── data/               # Input data (real data examples)
├── results/            # Output folders created by scripts
├── demo/               # Demo scripts for real and simulated analyses
├── README.md           # Project documentation
```

---

## 🧬 Real Data Demo

To run the full pipeline on real observed data:

### ✅ Option 1: Open as RStudio Project (Recommended)

1. Clone the repo.
2. In RStudio, go to **File → Open Project...** and select the `samSFS/` folder.
3. Open `demo/test_real_data.R` and run the code line-by-line.
4. Output will be saved to: `results/real_data_results/`

### ✅ Option 2: In RStudio (without project)

1. Open the `samSFS` folder in RStudio.
2. Open and run the file `demo/test_real_data.R`.
3. Run line-by-line to understand each step.
4. Results will be saved to: `results/real_data_results/`

### ✅ Option 3: From Terminal

```bash
cd samSFS
Rscript demo/test_real_data.R
```

---

## 🧪 Simulated Data Demo

To run a full simulation + inference pipeline:

### ✅ Option 1: Open as RStudio Project (Recommended)

1. Clone the repo.
2. In RStudio, go to **File → Open Project...** and select the `samSFS/` folder.
3. Open `demo/test_simulations.R` and run the code line-by-line.
4. Output will be saved to: `results/test_simulation_output/`

### ✅ Option 2: In RStudio (without project)

1. Open the `samSFS` folder in RStudio.
2. Open and run the file `demo/test_simulations.R`.
3. Step through each line to understand simulation and inference steps.
4. Output will be saved to: `results/test_simulation_output/`

### ✅ Option 3: From Terminal

```bash
cd samSFS
Rscript demo/test_simulations.R
```

---

## 📌 Notes

- All demo scripts assume you run them from the **root `samSFS/` folder**.
- Each script sets thread usage, simulation parameters, and output paths explicitly.
- Results folders are automatically created during analysis.

---

## 📫 Contact

For questions or collaboration, reach out via GitHub Issues or open a pull request.
