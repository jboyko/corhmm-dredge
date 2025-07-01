# corhmm-dredge
**Automatic Discovery of Optimal Discrete Character Models**

Automatically search discrete character evolution models and regularize parameter values using a combination of regularization and simulated annealing.

## Overview
This project implements automatic model selection for discrete character evolution in phylogenetic comparative methods. The "dredge" algorithm combines regularization techniques (L1, L2) with simulated annealing to optimize across different model structures without requiring user specification of model complexity.

## Simulation study and associated code
Each simulation setting has an identifier following a '-' and has a corresponding dataset in /data/ and model fit in /fits/. I vary several parameters for each script in order to test a variety of empirically relavent scenarios. Rates of evolution and tree size are two of the varried parameters. We are also interested in detecting particular model types when the data structure allows. For example, we may be interested in testing for correlation when there are two or more characters in the dataset. Finally we need to vary the number of rate classes. 
**code/simulate-01:** nChar=1, nStates=c(2), qRates=list(c(0.01,0.01), c(0.1, 0.1), c(1,1))

## analysis-param-est and simulate scripts
The scripts runs three main simulation scenarios (1, 2, and 3-character models), compares the estimated parameters from different model fitting approaches (e.g., `l0`, `l1`, `l2`, `er`) to the known "true" parameters, and then calculates various summary statistics (bias, variance, RMSE).

### Common Column Definitions

*   `ntips`: The number of tips in the phylogeny (e.g., 50, 100, 200).
*   `type`: The model fitting approach used for estimation (e.g., `er`, `l0`, `l1`, `l2`).
*   `par`: The specific transition rate parameter being estimated (e.g., `q[1][2]`).
*   `include`: A filter applied to the dataset based on the estimated rate values before summarizing.
    *   `all`: All estimates are included.
    *   `<50`: Only estimates with a rate less than 50 are included.
    *   `<10`: Only estimates with a rate less than 10 are included.
*   `model`: The complexity of the generating model (`1Char`, `2Char`, `3Char`).

---

### Results Organization
- **Parameter Results** (`param_results/`): Estimates by regularization type (L1/L2/ER)
- **Structure Results** (`structure_results/`): Model comparison outputs (dep_model, hmm_model, ord_model)
- **Summary Tables** (`tables/`): Performance metrics and parameter comparisons
- **Figures** (`figures/`, `plots/`): Publication-ready visualizations

### Key Output Files
- `parameter-comparison.csv`: Raw parameter estimates vs. true values
- `stats-summary-by-type.csv`: Bias, variance, RMSE by regularization method
- `propor-variance-by-type.csv`: Variance ratios relative to L0 baseline
- Model-specific summaries by character complexity (1Char/2Char/3Char)

### File Descriptions

#### `tables/parameter-comparison.csv`
The raw, unabridged data for every parameter estimate from every simulation run. This is the most granular file.
*   `value`: The log of the estimated parameter.
*   `true`: The log of the true (simulated) parameter.
*   `diff`: The difference between `value` and `true`.
*   `mse`: The squared difference between `value` and `true`.

#### `tables/par-summary-by-type.csv`
Summarizes estimation performance (mean, variance, median) grouped by tree size (`ntips`) and model type (`type`).
*   `estimate`, `diff`, `rmse`: These columns contain matrices of summary statistics (e.g., `[,1]` is the mean, `[,2]` is the variance, `[,3]` is the median) for the estimated values, their difference from the true value, and the Root Mean Squared Error.

#### `tables/stats-summary-by-type.csv`
A cleaner version of the summary, presenting the key performance metrics directly. This is useful for high-level comparisons.
*   `bias`: The mean difference between the estimated and true log-rates.
*   `variance`: The variance of the estimated log-rates.
*   `rmse`: The Root Mean Squared Error of the log-rates.

#### `tables/propor-variance-by-type.csv`
Shows the variance of estimates from model types `l1`, `l2`, and `er` as a proportion of the variance from the baseline `l0` model. Values greater than 1 indicate higher variance than the baseline.

#### `tables/par-summary-by-type-par.csv`
Provides summary statistics (mean, variance, median) grouped not just by `type` and `ntips`, but also by each specific parameter (`par`). This allows for a more detailed look at which specific rates are easy or hard to estimate.

#### `tables/par-summary-by-type-model.csv`
Provides summary statistics (mean, variance, median) grouped by the complexity of the true underlying model (`1Char`, `2Char`, `3Char`), in addition to `ntips` and `type`. This helps assess how model complexity affects estimation accuracy.

