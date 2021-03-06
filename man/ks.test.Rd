\name{ks.test}
\alias{ks.test}
\title{Kolmogorov-Smirnov Tests}
\usage{
ks.test(x, y, \dots, alternative = "two.sided")
}
\arguments{
    \item{x}{a numeric vector of data values.}
    \item{y}{either a numeric vector of data values, or a character
	string naming a distribution function.}
    \item{\dots}{parameters of the distribution specified by \code{y}.}
    \item{alternative}{indicates the alternative hypothesis and must be
	one of \code{"two.sided"}, \code{"greater"} or \code{"less"}.
	You can specify just the initial letter.}
}
\description{
    Performs one or two sample Kolmogorov-Smirnov tests.
}
\details{
    If \code{y} is numeric, a two sample test of the null that \code{x}
    and \code{y} were drawn from the same distribution is performed.

    Alternatively, \code{y} can be a character string naming a
    distribution function.  In this case, a one sample test of the null
    that the distribution function underlying \code{x} is \code{y} with
    parameters specified by \dots{} is carried out.
}
\value{
    A list with class \code{"htest"} containing the following
    components:
    \item{statistic}{the value of the test statistic.}
    \item{p.value}{the p-value of the test.}
    \item{alternative}{a character string describing the alternative
    	hypothesis.} 
    \item{method}{a character string indicating what type of test was
    	performed.} 
    \item{data.name}{a character string giving the name(s) of the data.}
}
\examples{
x <- rnorm(50)
y <- runif(30)
# Do x and y come from the same distribution?
ks.test(x, y)
# Does x come from a gamma distribution with shape 3 and scale 2?
ks.test(x, "pgamma", 3, 2)
}
\keyword{htest}
