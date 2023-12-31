% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rmfanova.R
\name{pointwise_ssa_test_statistic}
\alias{pointwise_ssa_test_statistic}
\title{Pointwise SSA test statistic}
\usage{
pointwise_ssa_test_statistic(
  x,
  plot = TRUE,
  values = FALSE,
  type = "l",
  ylab = "",
  main = "SSA(t)",
  ...
)
}
\arguments{
\item{x}{a list of length \eqn{\ell} with elements being \eqn{n\times p} matrices of data
corresponding to \eqn{n} functional observations measured in \eqn{p} design time points under given
experimental conditions.}

\item{plot}{a logical indicating of whether to draw the values of the pointwise SSA test statistic.
The default is \code{TRUE}.}

\item{values}{a logical indicating of whether to return the values of the pointwise SSA test statistic.
The default is \code{FALSE}.}

\item{type}{1-character string giving the type of plot desired, the same as in the \code{plot()} function.
The default is \code{"l"} for lines.}

\item{ylab}{a label for the \eqn{y}-axis, the same as in the \code{plot()} function.
The default is the empty sign.}

\item{main}{a main title for the plot, the same as in the \code{plot()} function. The default is \code{SSA(t)}.}

\item{...}{other graphical parameters, the same as in the \code{plot()} function.}
}
\value{
If \code{values = TRUE}, a vector of values of the pointwise SSA test statistic.
}
\description{
The function \code{pointwise_ssa_test_statistic()} calculates and draws the pointwise SSA test statistic.
}
\details{
For details, see the documentation of the \code{rmfanova()} function or
the paper Kurylo and Smaga (2023).
}
\examples{
# preparation of the DTI data set, for details see Kurylo and Smaga (2023)
library(refund)
data(DTI)
# MS patients
DTI_ms <- DTI[DTI$case == 1, ]
miss_data <- c()
for (i in 1:340) if (any(is.na(DTI_ms$cca[i, ]))) miss_data <- c(miss_data, i)
DTI_ms <- DTI_ms[-miss_data, ]
DTI_ms_2 <- DTI_ms[DTI_ms$Nscans == 4, ]
xx <- vector("list", 4)
for (i in 1:4) {
  xx[[i]] <- DTI_ms_2$cca[DTI_ms_2$visit == i, ]
}
xx[[1]] <- xx[[1]][-14, ]
xx[[3]] <- xx[[3]][-14, ]
yy <- xx
for (i in seq_len(4)) yy[[i]] <- yy[[i]][1:17, ]
# pointwise SSA test statistic
pointwise_ssa_test_statistic(yy, xlab = "t", xaxt = "n")
axis(1, c(1, 15, 30, 45, 60, 75, 93), labels = c(1, 15, 30, 45, 60, 75, 93))

}
\references{
Martinez-Camblor P., Corral N. (2011) Repeated Measures Analysis for Functional Data.
Computational Statistics & Data Analysis 55, 3244–3256.

Kurylo K., Smaga L. (2023) Functional repeated measures analysis of variance and its application.
Preprint https://arxiv.org/abs/2306.03883
}
