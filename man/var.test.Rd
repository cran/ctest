\name{var.test}
\alias{var.test}
\title{F Test to Compare Two Variances}
\usage{
var.test(x, y, ratio = 1, alternative = "two.sided", 
         conf.level = 0.95)
}
\arguments{
    \item{x, y}{numeric vectors of data values.}
    \item{ratio}{the hypothesized ratio of the population variances of
    	\code{x} and \code{y}.}
    \item{alternative}{indicates the alternative hypothesis and must be
    	one of \code{"two.sided"}, \code{"greater"} or \code{"less"}.
	You can specify just the initial letter.}
    \item{conf.level}{confidence level for the returned confidence
	interval.}
}
\description{
    Performs an F test to compare the variances of two samples from
    normal populations.
}
\details{
    The null hypothesis is that the ratio of the variances of the
    populations from which \code{x} and \code{y} were drawn is equal to
    \code{ratio}.
}
\value{
    A list with class \code{"htest"} containing the following
    components:
    \item{statistic}{the value of the F test statistic.}
    \item{parameter}{the degrees of the freedom of the F distribtion of
	the test statistic.}
    \item{p.value}{the p-value of the test.}
    \item{conf.int}{a confidence interval for the ratio of the
	population variances.}
    \item{estimate}{the ratio of the sample variances of \code{x} and
    	\code{y}.}
    \item{null.value}{the ratio of population variances under the null.}
    \item{alternative}{a character string describing the alternative
    	hypothesis.}
    \item{method}{the string \code{"F test to compare two variances"}.}
    \item{data.name}{a character string giving the names of the data.}
}
\keyword{htest}
