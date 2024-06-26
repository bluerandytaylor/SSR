% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/codeforSSR_q2.2.1.R
\name{ssr_bpp}
\alias{ssr_bpp}
\title{Sample size re-estimation based on bayesian predictive power}
\usage{
ssr_bpp(allo, n2plan, n2max, nnow, ptnow, pcnow, alpha, delta, target_power)
}
\arguments{
\item{allo}{the allocation rate into the treatment arm.}

\item{n2plan}{the initial planned sample size of two arms for stage2.}

\item{n2max}{the allowed largest sample size considering the medical resource.}

\item{nnow}{the total sample size of two arms at stage 1.}

\item{ptnow}{the calculated response rate in treatment arm at interim analysis.}

\item{pcnow}{the calculated response rate in control arm at interim analysis.}

\item{alpha}{the alpha level of the one-side test.}

\item{delta}{the difference of efficacy expected to be detected in the hyphothesises, usually 0.}

\item{target_power}{the target power wish to achieve.}
}
\description{
This function is for sample size re-estimation based on bayesian predictive power as in formula (23).
This function is based on the bayespredpower function.
}
\examples{
 library(SSreEstimate)
 ssr_bpp(1,100,400,100,0.58,0.4,0.025,0,0.8)

}
