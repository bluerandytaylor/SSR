% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/codeforSSR_q2.2.1.R
\name{bayespredpower}
\alias{bayespredpower}
\title{Bayesian predictive power calculation}
\usage{
bayespredpower(a, n2, n, pt, pc, alphalevel, h0delta)
}
\arguments{
\item{a}{the allocation ratio for patients.The patients are randomized into treatment arm and control arm under the ratio a:1}

\item{n2}{the initial planning total sample size of two arms at stage 2.}

\item{n}{the total sample size of two arms at stage 1.}

\item{pt}{the calculated response rate in treatment arm at interim analysis.}

\item{pc}{the calculated response rate in control arm at interim analysis.}

\item{alphalevel}{the alpha level of the one-side test.}

\item{h0delta}{the difference of efficacy expected to be detected in the hypothesises, usually 0.}
}
\description{
This function calculates the bayesian predictive power with non-informative prior for binary endpoints as in formula (21) and (22).
This function will be called for in the ssr_bpp function.
}
\examples{
library(SSreEstimate)
bpp <- bayespredpower(1,100,100,0.55,0.4,0.025,0)

}
