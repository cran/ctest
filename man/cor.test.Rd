\name{cor.test}
\alias{cor.test}
\title{Test for Zero Correlation}
\usage{
cor.test(x, y, alternative = "two.sided", method = "pearson",
         exact = NULL)
}
\arguments{
  \item{x, y}{numeric vectors of data values.  \code{x} and \code{y}
    must have the same length.}
  \item{alternative}{indicates the alternative hypothesis and must be
    one of \code{"two.sided"}, \code{"greater"} or \code{"less"}.  You
    can specify just the initial letter.}
  \item{method}{a string indicating which correlation coefficient is
    used for the test.  Must be one of \code{"pearson"},
    \code{"kendall"}, or \code{"spearman"}.  Only the first character is
    necessary.}
  \item{exact}{a logical indicating whether an exact p-value should be
    computed.}
}
\description{
  \code{cor.test} tests the null that \code{x} and \code{y} are
  uncorrelated (independent).
}
\details{
  If \code{method} is \code{"pearson"}, the test statistic is based on
  Pearson's product moment correlation coefficient \code{cor(x, y)} and
  follows a t distribution with \code{length(x)-2} degrees of freedom.

  If \code{method} is \code{"kendall"} or \code{"spearman"}, Kendall's
  tau or Spearman's rho, respectively, are used to estimate the
  correlation.  These tests should be used if the data do not
  necessarily come from a bivariate normal distribution.

  For Kendall's test, by default (if \code{exact} is not specified), an
  exact p-value is computed if both samples contain less than 50 finite
  values and there are no ties.  Otherwise, the standardized estimate is
  used as the test statistic, and is approximately normally distributed.

  For Spearman's test, p-values are computed using algorithm AS 89.
}
\value{
  A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the test statistic.}
  \item{parameter}{the degrees of freedom of the test statistic in the
    case that it follows a t distribution.}
  \item{p.value}{the p-value of the test.}
  \item{estimate}{the estimated correlation coefficient, with names
    attribute \code{"cor"}, \code{"tau"}, or \code{"rho"}, correspoding
    to the method employed.}
  \item{null.value}{the value of the correlation coefficient under the
    null hypothesis, hence \code{0}.}
  \item{alternative}{a character string describing the alternative
    hypothesis.}
  \item{method}{a string indicating how the correlation was estimated}
  \item{data.name}{a character string giving the names of the data.}
}
\references{
  D. J. Best & D. E. Roberts (1975),
  Algorithm AS 89: The Upper Tail Probabilities of Spearman's
  \eqn{\rho}{rho}.
  \emph{Applied Statistics}, \bold{24}, 377--379.
}
\keyword{htest}
