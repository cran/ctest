\name{kruskal.test}
\alias{kruskal.test}
\title{Kruskal-Wallis Rank Sum Test}
\usage{
kruskal.test(x, g)
}
\arguments{
    \item{x}{a numeric vector of data values, or a list of numeric data
    	vectors.}
    \item{g}{a vector or factor object giving the group for the
	corresponding elements of \code{x}.  Ignored if \code{x} is a
	list.}
}
\description{
    Performs a Kruskal-Wallis rank sum test.
}
\details{
    \code{kruskal.test} performs a Kruskal-Wallis rank sum test of the
    null that the location parameters of the distribution of \code{x}
    are the same in each group (sample).  The alternative is that they
    differ in at least one.

    If \code{x} is a list, its elements are taken as the samples to be
    compared, and hence have to be numeric data vectors.  In this case,
    \code{g} is ignored, and one can simply use \code{kruskal.test(x)}
    to perform the test.  If the samples are not yet contained in a
    list, use \code{kruskal.test(list(x, ...))}.

    Otherwise, \code{x} must be a numeric data vector, and \code{g} must
    be a vector or factor object of the same length as \code{x} giving
    the group for the corresponding elements of \code{x}.
}
\value{
    A list with class \code{"htest"} containing the following
    components:
    \item{statistic}{the Kruskal-Wallis rank sum statistic.}
    \item{parameter}{the degrees of freedom of the approximate
	chi-square distribution of the test statistic.}
    \item{p.value}{the p-value of the test.}
    \item{method}{the string \code{"Kruskal-Wallis rank sum test"}.}
    \item{data.name}{a character string giving the names of the data.}
}
\keyword{htest}
