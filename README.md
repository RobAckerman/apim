**apim**
apim is an R package providing convenience functions for the Actor-Partner Interdependence Model (APIM).

**Installation**
This package requires a patched version of glmmTMB that adds the
indisting() covariance structure for indistinguishable dyads. Please
install from RobAckerman/glmmTMB rather than the CRAN version of
glmmTMB.

```r
# Step 1: Install remotes if not already installed
install.packages("remotes")

# Step 2: Install patched glmmTMB with indisting() support
remotes::install_github("RobAckerman/glmmTMB", subdir = "glmmTMB")

# Step 3: Install apim from GitHub
remotes::install_github("RobAckerman/apim")

# Step 4: Load the packages
library(glmmTMB)
library(apim)
```
