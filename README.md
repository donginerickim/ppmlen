# ppmlen

`ppmlen` provides tools for estimating **Elastic-Net–regularized Poisson Pseudo Maximum Likelihood (PPML)** regressions with **high-dimensional fixed effects (HDFE)**.

This repository provides the methodological implementation used in Kim and Steinbach (2026).

The package builds on the **penalized PPML framework of Breinlich et al. (2021)** and extends it in two main directions:

1. **Elastic Net regularization**

   In addition to Lasso penalties, `ppmlen` allows Elastic Net regularization.  
   This improves estimation stability when explanatory variables are highly correlated.

2. **Efficient estimation with high-dimensional fixed effects**

   Fixed effects are removed using the **alternating projection (AP)** algorithm, which enables estimation with large sets of exporter, importer, dyad, and sector interactions typical in structural gravity models.

These features make the package particularly suitable for **gravity-style trade models with many policy variables**, where standard estimation methods may suffer from high dimensionality or multicollinearity.

---

# Installation

You can install the development version of `ppmlen` from GitHub:

```r
# install.packages("devtools")
devtools::install_github("donginerickim/ppmlen")
```

After installation:

```r
library(ppmlen)
```

---

# Example dataset

The package includes a small example dataset located at

```r
system.file("extdata", "sample.csv", package = "ppmlen")
```

The dataset is a restricted subset of the research data used in the associated study.

Dataset characteristics:

- **Years:** 1990–2015  
- **Sectors:**  
  - 12 — Fresh fruit  
  - 13 — Fresh vegetables  

The dependent variable `trade_f` represents a **forward (lead) trade flow variable**.

The dataset is intentionally small so that the example model can be estimated quickly.

---

# Example

The following example demonstrates how to estimate a gravity model with three sets of high-dimensional fixed effects.

```r
library(ppmlen)

csv_path <- system.file("extdata","sample.csv",package="ppmlen")
df <- read.csv(csv_path)

# Construct identifiers for fixed effects
df$ik <- interaction(df$exp, df$itpd_id, drop = TRUE)
df$jk <- interaction(df$imp, df$itpd_id, drop = TRUE)
df$ij <- interaction(df$exp, df$imp, drop = TRUE)

fit <- ppml_en_plugin(
  data = df,
  dep = "trade_f",
  fixed = list(
    c("ik","time"),
    c("jk","time"),
    c("ij","itpd_id")
  ),
  cluster = "ij",
  alpha = 0.5,
  K = 1,
  phipost = TRUE,
  verbose = TRUE,
  post = FALSE,
  always_keep = c("tau"),
  lasso_engine = "en",
  fes_mode = "ap"
)

fit$selected
```

---

# Main arguments

The function `ppml_en_plugin()` includes several key arguments.

**data**

Dataset used for estimation.

**dep**

Dependent variable.

**fixed**

A list describing the sets of high-dimensional fixed effects.  
In gravity models these often correspond to:

- exporter × time × sector  
- importer × time × sector  
- dyad × sector

**cluster**

Clustering variable used for variance estimation.

**alpha**

Elastic Net mixing parameter.

- `alpha = 1` corresponds to **Lasso**
- `alpha < 1` introduces a **ridge component**

**K**

Number of Elastic Net iterations used in the estimation procedure.

**always_keep**

Variables that must always remain in the model, even when regularization is applied.

For example:

```r
always_keep = c("tau")
```

ensures that the tariff variable `tau` is always included in the regression, while other variables are subject to Elastic Net selection.

**lasso_engine**

Regularization method used in the estimation.

Currently supported:

- `"en"` — Elastic Net

**fes_mode**

Method used to remove fixed effects.

- `"ap"` uses the **alternating projection algorithm**.

---

# Methodological background

The implementation relies on several methodological contributions.

**Penalized PPML estimation**

Breinlich, H., Corradi, V., Rocha, N., Ruta, M., Santos Silva, J. M. C., & Zylkin, T. (2021).  
Machine learning in international trade research: Evaluating the impact of trade agreements.

**Alternating projection algorithm for high-dimensional fixed effects**

Gaure, S. (2013).  
OLS with multiple high dimensional category variables.  
Computational Statistics & Data Analysis.

**Coordinate descent for penalized regression**

Friedman, J., Hastie, T., & Tibshirani, R. (2010).  
Regularization paths for generalized linear models via coordinate descent.  
Journal of Statistical Software.

---

# Citation

If you use this package, please cite:

Kim, Dongin, and Sandro Steinbach.  
*"Non-Tariff PTA Provisions and Agricultural Trade: New Insights from a Structural Gravity Analysis and Machine Learning."*  
**American Journal of Agricultural Economics**, forthcoming, 2026.

---

# License

MIT License.
