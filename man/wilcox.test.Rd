\name{wilcox.test}
\alias{wilcox.test}
\title{Wilcoxon Rank Sum and Signed Rank Tests}
\usage{
wilcox.test(x, y = NULL, alternative = "two.sided", mu = 0, 
            paired = FALSE, exact = NULL, correct = TRUE)
}
\arguments{
    \item{x}{numeric vector of data values.}
    \item{y}{an optional numeric vector of data values.}
    \item{alternative}{indicates the alternative hypothesis and must be
    	one of \code{"two.sided"}, \code{"greater"} or \code{"less"}.
	You can specify just the initial letter.}
    \item{mu}{a number specifying an optional location parameter.}
    \item{paired}{a logical indicating whether you want a paired test.}
    \item{exact}{a logical indicating whether an exact p-value should be
	computed.}
    \item{correct}{a logical indicating whether to apply continuity
    	correction in the normal approximation for the p-value.}
}
\description{
    Performs one and two sample Wilcoxon tests on vectors of data.
}
\details{
    If only \code{x} is given, or if both \code{x} and \code{y} are
    given and \code{paired} is \code{TRUE}, a Wilcoxon signed rank test
    of the null that the median of \code{x} (in the one sample case) or
    of \code{x-y} (in the paired two sample case) equals \code{mu} is
    performed.

    Otherwise, if both \code{x} and \code{y} are given and \code{paired}
    is \code{FALSE}, a Wilcoxon rank sum test (equivalent to the
    Mann-Whitney test) is carried out.  In this case, the null
    hypothesis is that the location of the distributions of \code{x} and
    \code{y} differ by \code{mu}.

    By default (if \code{exact} is not specified), an exact p-value is
    computed if the samples contain less than 50 finite values and
    there are no ties.  Otherwise, a normal approximation is used.
}
\value{
    A list with class \code{"htest"} containing the following
    components:
    \item{statistic}{the value of the test statistic with a name
    	describing it.}
    \item{parameter}{the parameter(s) for the exact distribution of the
    	test statistic. Currently, only normal approximations are used.}
    \item{p.value}{the p-value for the test.}
    \item{null.value}{the location parameter \code{mu}.}
    \item{alternative}{a character string describing the alternative
    	hypothesis.}
    \item{method}{the type of test applied.}
    \item{data.name}{a character string giving the names of the data.}
}
\keyword{htest}
