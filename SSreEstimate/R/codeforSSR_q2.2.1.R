#' @title Conditional power calculation
#' @description
#' This function calculates the conditional power for binary endpoints as in formula (9) and (10).
#' This function will be called for in the ssr_cp function.
#'
#' @param a the allocation ratio for patients.The patients are randomized into treatment arm and control arm under the ratio a:1
#' @param n2 the initial planning total sample size of two arms at stage 2.
#' @param n the total sample size of two arms at stage 1.
#' @param pt the calculated response rate in treatment arm at interim analysis.
#' @param pc the calculated response rate in control arm at interim analysis.
#' @param alphalevel the alpha level of the one-side test.
#' @param h0delta the difference of efficacy expected to be detected in the hypothesises, usually 0.
#'
#' @examples
#' library(SSreEstimate)
#' cp <- conditionalpower(1,100,100,0.55,0.4,0.025,0)
#'
#' @export
conditionalpower <- function(a,n2,n,pt,pc,alphalevel,h0delta){
  r <-((a+1)^2/a)^(1/2)
  sn2 <- a/(a+1)*((pt*(1-pt))/a+pc*(1-pc))
  sn <- sn2^(1/2)
  gama <- qnorm(1-alphalevel)
  value <- 1/(r*sn)*((n2+n)/n2)^(1/2)*((pt-pc-h0delta)*(n+n2)^(1/2)-r*sn*gama)
  cp <- pnorm(value)
  return(cp)
}


#' @title Sample size re-estimation based on conditional power
#' @description
#' This function is for sample size re-estimation based on conditional power as in formula (11).
#' This function is based on the conditionalpower function.
#'
#' @param allo the allocation rate into the treatment arm.
#' @param n2plan the initial planned sample size of two arms for stage2.
#' @param n2max the allowed largest sample size considering the medical resource.
#' @param nnow the total sample size of two arms at stage 1.
#' @param ptnow the calculated response rate in treatment arm at interim analysis.
#' @param pcnow the calculated response rate in control arm at interim analysis.
#' @param alpha the alpha level of the one-side test.
#' @param delta the difference of efficacy expected to be detected in the hyphothesises, usually 0.
#' @param target_power the target power wish to achieve.
#'
#' @examples
#'  library(SSreEstimate)
#'  ssr_cp(1,100,300,100,0.55,0.4,0.025,0,0.8)
#'
#' @export
ssr_cp <- function(allo,
                   n2plan,
                   n2max,
                   nnow,
                   ptnow,
                   pcnow,
                   alpha,
                   delta,
                   target_power){
  size <- n2plan
  repeat{
    cp <- conditionalpower(a=allo,
                           n2=size,
                           n=nnow,
                           pt=ptnow,
                           pc=pcnow,
                           alphalevel = alpha,
                           h0delta = delta)
    cat("size=",size," conditional power=",cp,"\n")
    if (target_power>cp){size <- size+2}
    if (cp>=target_power){break}
    if (size>=n2max){break}
  }
  if (size>n2max){size <- n2max}
  return(size)
}

#' @title  Conditional power calculation using Wilson-Hilferty transformation
#' @description
#' This function calculates the conditional power under CMH test using Wilson-Hilferty transformation  as in formula (14) and (15).
#' This function will be called for in the ssr_cp_wh function.
#'
#' @param psi1 the calculated CMH test statistics from the CMH test at interim analysis.
#' @param n2 the initial planning total sample size of two arms at stage 2.
#' @param n the total sample size of two arms at stage 1.
#' @param alphalevel the alpha level of the one-side test.
#'
#' @examples
#' library(SSreEstimate)
#' cp_wh <- conditionalpower_wh(4,100,100,0.025)
#'
#' @export
conditionalpower_wh <- function(psi1,n2,n,alphalevel){
  z_wh1 <- (psi1^(1/3)-(1-2/9))/(2/9)^(1/2)
  gama <- qnorm(1-alphalevel)
  time_fraction <- n/(n2+n)
  value <- (1-time_fraction)^(-1/2)*(z_wh1/time_fraction^(1/2)-gama)
  cp_wh <- pnorm(value)
  return(cp_wh)
}


#' @title Sample size re-estimation based on conditional power using  Wilson-Hilferty transformation
#' @description
#' This function is for sample size re-estimation with CMH test based on conditional power and Wilson-Hilferty transformation as in formula (16).
#' This function is based on the conditionalpower_wh function.
#'
#' @param psi1value the calculated CMH test statistics from the CMH test at interim analysis.
#' @param n2plan the initial planning total sample size of two arms at stage 2.
#' @param n2max the allowed largest sample size considering the medical resource.
#' @param nnow the total sample size of two arms at stage 1.
#' @param alpha the alpha level of the one-side test.
#' @param target_power the target power wish to achieve.
#'
#' @examples
#' library(SSreEstimate)
#' ssr_cp_wh(5,100,300,200,0.025,0.8)
#'
#' @export
ssr_cp_wh <- function(psi1value,
                      n2plan,
                      n2max,
                      nnow,
                      alpha,
                      target_power){
  size <- n2plan
  repeat{
    cp_wh <- conditionalpower_wh(psi1 = psi1value,
                                 n2=size,
                                 n=nnow,
                                 alphalevel = alpha)
    cat("size=",size," conditional power based on WH transformation is ",cp_wh,"\n")
    if (target_power>cp_wh){size <- size+2}
    if (cp_wh>=target_power){break}
    if (size>=n2max){break}
  }
  if (size>n2max){size <- n2max}
  return(size)
}


#' @title Bayesian predictive power calculation
#' @description
#' This function calculates the bayesian predictive power with non-informative prior for binary endpoints as in formula (21) and (22).
#' This function will be called for in the ssr_bpp function.
#'
#' @param a the allocation ratio for patients.The patients are randomized into treatment arm and control arm under the ratio a:1
#' @param n2 the initial planning total sample size of two arms at stage 2.
#' @param n the total sample size of two arms at stage 1.
#' @param pt the calculated response rate in treatment arm at interim analysis.
#' @param pc the calculated response rate in control arm at interim analysis.
#' @param alphalevel the alpha level of the one-side test.
#' @param h0delta the difference of efficacy expected to be detected in the hypothesises, usually 0.
#'
#' @examples
#' library(SSreEstimate)
#' bpp <- bayespredpower(1,100,100,0.55,0.4,0.025,0)
#'
#' @export
bayespredpower <- function(a,n2,n,pt,pc,alphalevel,h0delta){
  time_fraction <- n/(n+n2)
  r <-((a+1)^2/a)^(1/2)
  sn2 <- a/(a+1)*((pt*(1-pt))/a+pc*(1-pc))
  sn <- sn2^(1/2)
  se_theta1hat <- r*sn/n^(1/2)
  z1 <- (pt-pc-h0delta)/se_theta1hat
  gama <- qnorm(1-alphalevel)
  value <- (1-time_fraction)^(-1/2)*(z1/time_fraction^(1/2)-gama)*time_fraction^(1/2)
  bpp <- pnorm(value)
  return(bpp)
}


#' @title Sample size re-estimation based on bayesian predictive power
#' @description
#' This function is for sample size re-estimation based on bayesian predictive power as in formula (23).
#' This function is based on the bayespredpower function.
#'
#' @param allo the allocation rate into the treatment arm.
#' @param n2plan the initial planned sample size of two arms for stage2.
#' @param n2max the allowed largest sample size considering the medical resource.
#' @param nnow the total sample size of two arms at stage 1.
#' @param ptnow the calculated response rate in treatment arm at interim analysis.
#' @param pcnow the calculated response rate in control arm at interim analysis.
#' @param alpha the alpha level of the one-side test.
#' @param delta the difference of efficacy expected to be detected in the hyphothesises, usually 0.
#' @param target_power the target power wish to achieve.
#'
#' @examples
#'  library(SSreEstimate)
#'  ssr_bpp(1,100,400,100,0.58,0.4,0.025,0,0.8)
#'
#' @export
ssr_bpp <- function(allo,#allocation ratio
                    n2plan,
                    n2max,
                    nnow,
                    ptnow,
                    pcnow,
                    alpha,
                    delta,
                    target_power){
  size <- n2plan
  repeat{
    bpp <- bayespredpower(a=allo,
                          n2=size,
                          n=nnow,
                          pt=ptnow,
                          pc=pcnow,
                          alphalevel = alpha,
                          h0delta = delta)
    cat("size=",size," bayesian predictive power=",bpp,"\n")
    if (target_power>bpp){size <- size+2}
    if (bpp>=target_power){break}
    if (size>=n2max){break}
  }
  if (size>n2max){size <- n2max}
  return(size)

}

#' @title Bayesian predictive power calculation using  Wilson-Hilferty transformation
#' @description
#' This function calculates the bayesian predictive power under CMH test using non-informative prior and Wilson-Hilferty transformation as in formula (24) and (25).
#' This function will be called for in the ssr_bpp_wh function.
#'
#' @param psi1 the calculated CMH test statistics from the CMH test at interim analysis.
#' @param n2 the initial planning total sample size of two arms at stage 2.
#' @param n the total sample size of two arms at stage 1.
#' @param alphalevel the alpha level of the one-side test.
#'
#' @examples
#' library(SSreEstimate)
#' bpp_wh <- bayespredpower_wh(4,100,100,0.025)
#'
#' @export
bayespredpower_wh <- function(psi1,n2,n,alphalevel){
  z_wh1 <- (psi1^(1/3)-(1-2/9))/(2/9)^(1/2)
  gama <- qnorm(1-alphalevel)
  time_fraction <- n/(n2+n)
  value <- (1-time_fraction)^(-1/2)*(z_wh1/time_fraction^(1/2)-gama)*time_fraction^(1/2)
  bpp_wh <- pnorm(value)
  return(bpp_wh)
}


#' @title Sample size re-estimation based on bayesian predictive power using  Wilson-Hilferty transformation
#' @description
#' This function is for sample size re-estimation with CMH test based on bayesian predictive power and Wilson-Hilferty transformation as in formula (26).
#' This function is based on the bayespredpower_wh function.
#'
#' @param psi1value the calculated CMH test statistics from the CMH test at interim analysis.
#' @param n2plan the initial planning total sample size of two arms at stage 2.
#' @param n2max the allowed largest sample size considering the medical resource.
#' @param alpha the alpha level of the one-side test.
#' @param nnow the total sample size of two arms at stage 1.
#' @param target_power the target power wish to achieve.
#'
#' @examples
#' library(SSreEstimate)
#' ssr_bpp_wh(5,100,300,200,0.025,0.8)
#'
#' @export

ssr_bpp_wh <- function(psi1value,
                       n2plan,
                       n2max,
                       nnow,
                       alpha,
                       target_power){
  size <- n2plan
  repeat{
    bpp_wh <- bayespredpower_wh(psi1 = psi1value,
                                 n2=size,
                                 n=nnow,
                                 alphalevel = alpha)
    cat("size=",size," bayesian predictive power based on WH transformation is ",bpp_wh,"\n")
    if (target_power>bpp_wh){size <- size+2}
    if (bpp_wh>=target_power){break}
    if (size>=n2max){break}
  }
  if (size>n2max){size <- n2max}
  return(size)

}

