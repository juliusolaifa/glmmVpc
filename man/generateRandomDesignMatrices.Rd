% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/design_matrices.R
\name{generateRandomDesignMatrices}
\alias{generateRandomDesignMatrices}
\title{Generate Random Design Matrices for Mixed Models}
\usage{
generateRandomDesignMatrices(ns, X = NULL)
}
\arguments{
\item{ns}{A vector of group sizes, indicating the number of
observations per group.}

\item{X}{Optional. A matrix or vector of fixed effect predictors.
The first column is assumed to be an intercept if all values are 1.}
}
\value{
A list containing the intercept and slope matrices. If \code{X} is NULL,
only the intercept matrix is returned.
}
\description{
This function generates random design matrices for
random intercepts and random slopes. If \code{X} contains a column of 1's as the
first column (indicating an intercept), the intercept matrix is not generated
separately. If \code{X} is NULL, only the intercept matrix is generated.
}
\examples{
ns <- c(2, 3, 2)
X <- matrix(c(rep(1,7),seq(2,20,length.out=7)), nrow = 7, ncol = 2)
print(X)
generateRandomDesignMatrices(ns,X)
}
