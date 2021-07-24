\name{ALR-package}
\alias{ALR-package}
\alias{ALR}
\docType{package}
\title{
  Alternating Logistic Regression
}
\description{
  Fit alternating logistic regression to 
}
\details{
  The package allow modelling longitudinal binary outcome with alternating logistic regression.
}
\author{
Yuqi Tian

}
\references{
Carey, V., Zeger, S. L., & Diggle, P. (1993). Modelling multivariate binary data with alternating logistic regressions. Biometrika, 80(3), 517-526.
}
\keyword{ package }

\examples{
  \dontrun{
     sample_size <- 200
     cluster_size <- rep(10, sample_size)
     z <- matrix(rep(1, sum(choose(cluster_size,2))), ncol=1) 
     mod <- ALR(mean.formula = y ~ time * x, 
         data = data, 
         z = z,
         n = cluster_size)
  }
}
