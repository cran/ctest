\name{friedman.test}
\alias{friedman.test}
\title{Friedman Rank Sum Test}
\usage{
friedman.test(y, groups, blocks)
}
\arguments{
    \item{y}{either a numeric vector of data values, or a data matrix.}
    \item{groups}{a vector giving the group for the corresponding
	elements of \code{y} if this is a vector;  ignored if \code{y}
	is a matrix.  If not a factor object, it is coerced to one.}
    \item{blocks}{a vector giving the block for the corresponding
	elements of \code{y} if this is a vector;  ignored if \code{y}
	is a matrix.  If not a factor object, it is coerced to one.}
}
\description{
    Performs a Friedman rank sum test with unreplicated blocked data.
}
\details{
    \code{friedman.test} can be used for analyzing unreplicated complete
    block designs (i.e., there is exactly one observation in \code{y}
    for each combination of levels of \code{groups} and \code{blocks})
    where the normality assumption may be violated.

    The null hypothesis is that apart from an effect of \code{blocks},
    the location parameter of \code{y} is the same in each of the
    \code{groups}. 

    If \code{y} is a matrix, \code{groups} and \code{blocks} are
    obtained from the column and row indices, respectively.  \code{NA}'s
    are not allowed in \code{groups} or \code{blocks};  if \code{y}
    contains \code{NA}'s, corresponding blocks are removed.
}
\value{
    A list with class \code{"htest"} containing the following
    components:
    \item{statistic}{the value of Friedman's chi-square statistic.}
    \item{parameter}{the degrees of freedom of the approximate
	chi-square distribution of the test statistic.}
    \item{p.value}{the p-value of the test.}
    \item{method}{the string \code{"Friedman rank sum test"}.}
    \item{data.name}{a character string giving the names of the data.}
}
\keyword{htest}
