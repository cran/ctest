\name{mood.test}
\alias{mood.test}
\title{Mood Two-Sample Test of Scale}
\usage{
mood.test(x, y, alternative = "two.sided")
}
\arguments{
    \item{x, y}{numeric vectors of data values.}
    \item{alternative}{indicates the alternative hypothesis and must be
    	one of \code{"two.sided"}, \code{"greater"} or \code{"less"}.
	You can specify just the initial letter.}
}
\description{
    Performs Mood's two-sample test of scale.
}
\details{
    The underlying model is that the two samples are drawn from
    \eqn{f(x-l)} and \eqn{f((x-l)/s)/s}, respectively, where \eqn{l} is
    a common location parameter and \eqn{s} is a scale parameter.

    The null hypothesis is \eqn{s = 1}.

    There are more useful tests for this problem.
}
\value{
    A list with class \code{"htest"} containing the following
    components:
    \item{statistic}{the value of the test statistic.}
    \item{p.value}{the p-value of the test.}
    \item{alternative}{a character string describing the alternative
    	hypothesis.}
    \item{method}{the string \code{"Mood two-sample test of scale"}.}
    \item{data.name}{a character string giving the names of the data.}
}
\keyword{htest}
