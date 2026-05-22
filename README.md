**apim**
apim is an R package providing convenience functions for the Actor-Partner Interdependence Model (APIM).

## Installation

### Step 1: Install the patched glmmTMB (required for longitudinal models with indistinguishable dyads)

apim requires a patched version of `glmmTMB` that adds the `indisting()` covariance structure for indistinguishable dyads. This must be installed before apim.

> **Note:** This replaces the CRAN version of `glmmTMB` on your machine. The patched version is fully backward compatible — all existing `glmmTMB` functionality is preserved.

> **Note:** Installing from GitHub requires compilation tools. On Windows, install [Rtools45](https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html) first. On Mac, install Xcode Command Line Tools (`xcode-select --install`). On Linux, install `build-essential`.

```r
install.packages("remotes")
remotes::install_github("RobAckerman/glmmTMB", subdir = "glmmTMB")
```

Restart R after this step.

### Step 2: Install apim

```r
install.packages("apim", repos = "https://robackerman.r-universe.dev")
```

### Step 3: Restart R and load

```r
library(glmmTMB)
library(apim)
```

### If you do not need indistinguishable dyad models

If you only need functions for distinguishable dyads and do not use `indisting()`, you can install apim without the patched `glmmTMB`:

```r
install.packages("remotes")
remotes::install_github("RobAckerman/apim", dependencies = FALSE)
```

All apim functions except those requiring `indisting()` will work with the CRAN version of `glmmTMB`.
