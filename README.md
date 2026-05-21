**apim**
apim is an R package providing convenience functions for the Actor-Partner Interdependence Model (APIM).

## Installation

### Recommended (fast, no compilation required)

```r
install.packages("apim", repos = "https://robackerman.r-universe.dev")
```

> **Note:** This installs a patched version of `glmmTMB` that adds the
> `indisting()` covariance structure for indistinguishable dyads. This
> replaces the CRAN version of `glmmTMB` on your machine. The patched
> version is fully backward compatible — all existing `glmmTMB`
> functionality is preserved.

### If you do not want to replace your glmmTMB installation

Install from GitHub without the r-universe repo. You will still need to
install the patched `glmmTMB` separately if you want to use `indisting()`,
but all other `apim` functions will work with the CRAN version of `glmmTMB`.

```r
install.packages("remotes")
remotes::install_github("RobAckerman/apim", dependencies = FALSE)
```

### Load the packages

```r
library(glmmTMB)
library(apim)
```
```
