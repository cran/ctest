\name{binom.test}
\alias{binom.test}
\title{Exact Binomial Test}
\usage{
binom.test(x, n, p = 0.5, alternative = "two.sided")
}
\arguments{
    \item{x}{number of successes.}
    \item{n}{number of trials.}
    \item{p}{probability of success.}
    \item{alternative}{indicates the alternative hypothesis and must be
	one of \code{"two.sided"}, \code{"greater"} or \code{"less"}.
	You can specify just the initial letter.}
}
\description{
    \code{binom.test} performs an exact test of the null that the
    probability of success in a Bernoulli experiment of length \code{n}
    is \code{p}, based on the number \code{x} of successes observed.
}
\value{
    A list with class \code{"htest"} containing the following
    components:
    \item{statistic}{the number of successes, \code{x}.}
    \item{parameter}{the number of trials, \code{n}.}
    \item{p.value}{the p-value of the test.}
    \item{null.value}{the probability of success under the null,
    	\code{p}.}
    \item{alternative}{a character string describing the alternative
    	hypothesis.}
    \item{method}{the string \code{"Exact binomial test"}.}
    \item{data.name}{a character string giving the names of the data.}
}
\keyword{htest}
