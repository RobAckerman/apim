\name{summary_glmmTMB_satterthwaite}
\alias{summary_glmmTMB_satterthwaite}
\title{Satterthwaite Summary for glmmTMB Models}
\description{
Produces a formatted summary for a fitted \code{glmmTMB} model using
Satterthwaite degrees of freedom for fixed effects. Confidence intervals
are computed using Satterthwaite df and included in the fixed effects table.
}
\usage{
summary_glmmTMB_satterthwaite(
  model,
  eps = 0.001,
  digits = 5,
  alpha = 0.05,
  verbose = FALSE
)
}
\arguments{
\item{model}{A fitted \code{glmmTMB} model object.}
\item{eps}{Numeric. Step size for central differences in the Jacobian.
Only used when \code{dispformula = ~0}. Default is \code{1e-3}.}
\item{digits}{Integer. Number of digits for rounding. Default is \code{5}.}
\item{alpha}{Numeric. Significance level for confidence intervals.
Default is \code{0.05}, giving 95\% confidence intervals.}
\item{verbose}{Logical. If \code{TRUE}, prints progress messages.
Default is \code{FALSE}.}
}
\value{
Invisibly returns a list with elements: satterthwaite (data frame with
Term, Estimate, SE, df, t, p, lower, upper), vc (VarCorr output), and
fit_stats (named vector of AIC, BIC, logLik, -2*log(L), df.resid).
}
\examples{
\dontrun{
m <- glmmTMB::glmmTMB(y ~ x, data = dat)
summary_glmmTMB_satterthwaite(m)
summary_glmmTMB_satterthwaite(m, alpha = 0.10)
}
}
