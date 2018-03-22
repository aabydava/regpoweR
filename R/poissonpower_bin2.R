
#' Sample size and power for Poisson regression.
#'
#' A function for calculating sample size and power based on Poisson regression. The function solves for one of the following: alpha, power, incidence rate ratio, N.
#' @param alpha alpha level (usually 0.05). Can range from 0 to 1
#' @param power 1 - Pr(type II error) (usually 0.80). Can range from 0 to 1.
#' @param exp.B0 Baseline rate. The response rate when all covariates have a value of 0
#' @param exp.B1 Incidence rate ratio (IRR). Also can be called the IRR.
#' @param muT Mean exposure time
#' @param phi Measure of over-dispersion
#' @param R2 The square of the multiple correlation coefficient when the covariate of interst is regressed on the other covariates.
#' @param pi.x1 Percent of N with X1=1
#' @param N sample size
#' @return Return one of the following parameters
#' @return \code{alpha}
#' @return \code{power}
#' @return \code{exp.B1}
#' @return \code{N}
#' @note The test is a two-sided test. For one-sided tests, double the significance level. For example, you can set alpha=0.10 to obtain one-sided test at 0.05 significance level.
#' @author David Aaby <david.aaby@@northwestern.edu>
#' @references Signorini, D. 1991. Sample Size for Poisson Regression, Biometrika (1991), 78, 2, pages 446-450..
#' @examples
#' Poissonpower.bin(alpha=NULL, power=.80, exp.B0=.85, exp.B1=1.3, muT=1, phi=1,R2=0, pi.x1=.5, N=406)
#' Poissonpower.bin(alpha=NULL, power=.80, exp.B0=.85, exp.B1=1.3, muT=1, phi=1, R2=0, pi.x1=.5, N=406)
#' Poissonpower.bin(alpha=.05, power=.80, exp.B0=.85, exp.B1=1.3, muT=1, phi=1, R2=0, pi.x1=.5, N=NULL)
#' Poissonpower.bin(alpha=.05, power=.80, exp.B0=.85, exp.B1=NULL, muT=1, phi=1, R2=0, pi.x1=.5, N=406)


Poissonpower.bin <- function(alpha=NULL, power=NULL, exp.B0=NULL, exp.B1=NULL, muT=1,
                             phi=1, R2=0, pi.x1=NULL, N=NULL) {

  ##################
  # error messages #
  ##################
  if(any(alpha < 0 | alpha > 1))    stop('alpha not between 0 and 1')
  if(any(power < 0 | power > 1))    stop('power not between 0 and 1')
  if(any(N < 2))                    stop('N is less than 2')
  if(any(muT < 1))                  stop('muT must be >= 1')
  if(any(phi < 1))                  stop('phi must be >= 1')
  if(any(R2 < 0))                   stop('R2 must be >= 0')
  if(any(pi.x1 < 0 | pi.x1 > 1))    stop('pi.x1 not between 0 and 1')

  if(length(muT) > 1)               stop('muT must be single value')
  if(length(phi) > 1)               stop('phi must be single value')
  if(length(R2) > 1)                stop('R2 must be single value')


  mylist = list(alpha, power, exp.B0, exp.B1, pi.x1, N)
  l.mylist = lengths(mylist)
  if(length(l.mylist[l.mylist>=2]) > 1) stop('Only vary one parameter at a time')


  # solve for alpha level #
  if(is.null(alpha)) {

    B0 = log(exp.B0)
    B1 = log(exp.B1)

    beta = 1 - power

    var.x1 = pi.x1*(1-pi.x1)

    var.beta1.0 = 1 / var.x1
    var.beta1.B1 =  (1 / (1-pi.x1)) + (1 / (pi.x1*exp.B1))

    A1 = sqrt((N * muT * exp.B0 * (B1^2) * (1-R2)) / phi)
    A2 = stats::qnorm(1-beta) * sqrt(var.beta1.B1)

    A3 = sqrt(var.beta1.0)
    A4 = (A1 - A2) / A3

    alpha = 2*(1-stats::pnorm(A4))     # two-sided test
    #alpha = (1-stats::pnorm(A4))       # one-sided test

    # output results #
    results = NULL

    if(length(power) > 1 | length(exp.B0) > 1 | length(exp.B1) > 1 | length(pi.x1) > 1 | length(N) > 1)  {
      rownum = max(length(alpha), length(power), length(exp.B0))
      for(i in 1:rownum) {
        results = cbind(alpha, power, N, exp.B0, exp.B1, pi.x1)
      }
      colnames(results) = c("alpha", "power", "N", "baseline rate (exp.B0)", "IRR (exp.B1)", "X1 proportion (pi.x1)")
    }

    else {
      results = c(alpha, power, N, exp.B0, exp.B1, pi.x1)
      results = data.frame(results)
      results = t(results)
      colnames(results) = c("alpha", "power", "N", "baseline rate (exp.B0)", "IRR (exp.B1)",  "X1 proportion (pi.x1)")
    }
  }



  # solve for power #
  if(is.null(power)) {

    B0 = log(exp.B0)
    B1 = log(exp.B1)

    beta = 1 - power

    var.x1 = pi.x1*(1-pi.x1)

    var.beta1.0 = 1 / var.x1
    var.beta1.B1 =  (1 / (1-pi.x1)) + (1 / (pi.x1*exp.B1))

    A1 = sqrt((N * muT * exp.B0 * (B1^2) * (1-R2)) / phi)
    A2 = stats::qnorm(1-alpha/2) * sqrt(var.beta1.0)    # two-sided test
    #A2 = stats::qnorm(1-alpha) * sqrt(var.beta1.0)      # one-sided test

    A3 = sqrt(var.beta1.B1)
    A4 = (A1 - A2) / A3

    power = stats::pnorm(A4)     # two-sided test


    # output results #
    results = NULL

    if(length(alpha) > 1 | length(exp.B0) > 1 | length(exp.B1) > 1 | length(pi.x1) > 1 | length(N) > 1)  {
      rownum = max(length(alpha), length(power), length(exp.B0))
      for(i in 1:rownum) {
        results = cbind(alpha, power, N, exp.B0, exp.B1, pi.x1)
      }
      colnames(results) = c("alpha", "power", "N", "baseline rate (exp.B0)", "IRR (exp.B1)", "X1 proportion (pi.x1)")
    }

    else {
      results = c(alpha, power, N, exp.B0, exp.B1, pi.x1)
      results = data.frame(results)
      results = t(results)
      colnames(results) = c("alpha", "power", "N", "baseline rate (exp.B0)", "IRR (exp.B1)",  "X1 proportion (pi.x1)")
    }
  }



  # solve for sample size N #
  if(is.null(N)) {

    B0 = log(exp.B0)
    B1 = log(exp.B1)

    beta = 1 - power

    var.x1 = pi.x1*(1-pi.x1)

    var.beta1.0 = 1 / var.x1
    var.beta1.B1 =  (1 / (1-pi.x1)) + (1 / (pi.x1*exp.B1))

    A1 = stats::qnorm(1-alpha/2) * sqrt(var.beta1.0)     # two-sided test
    #A1 = stats::qnorm(1-alpha) * sqrt(var.beta1.0)     # one-sided test
    A2 = stats::qnorm(1-beta) * sqrt(var.beta1.B1)
    A3 =  muT * exp.B0 * (B1^2) * (1-R2)

    N = phi * (A1 + A2)^2 / A3
    N = ceiling(N)



    # output results #
    results = NULL

    if(length(alpha) > 1 | length(power) > 1 | length(exp.B0) > 1 | length(exp.B1) > 1 | length(pi.x1) > 1) {
      rownum = max(length(alpha), length(power), length(exp.B0))
      for(i in 1:rownum) {
        results = cbind(alpha, power, N, exp.B0, exp.B1, pi.x1)
      }
      colnames(results) = c("alpha", "power", "N", "baseline rate (exp.B0)", "IRR (exp.B1)", "X1 proportion (pi.x1)")
    }

    else {
      results = c(alpha, power, N, exp.B0, exp.B1, pi.x1)
      results = data.frame(results)
      results = t(results)
      colnames(results) = c("alpha", "power", "N", "baseline rate (exp.B0)", "IRR (exp.B1)",  "X1 proportion (pi.x1)")
    }
  }

  # solve for incidence rate ratio #
  if(is.null(exp.B1)) {
    N1 = N

    beta = 1 - power
    exp.B1 = seq(1, 5, .00001)
    exp.B1 = exp.B1[-1]

    maxlength = max(c(length(alpha), length(power), length(exp.B0), length(N)))

    vec.exp.B1 = NULL

    if(maxlength==1) {

      for(i in 1:length(exp.B0)) {
        B0 = log(exp.B0[i])
        B1 = log(exp.B1)

        var.x1 = pi.x1*(1-pi.x1)
        var.beta1.0 = 1 / var.x1
        var.beta1.B1 =  (1 / (1-pi.x1)) + (1 / (pi.x1*exp.B1))

        A1 = stats::qnorm(1-alpha/2) * sqrt(var.beta1.0)
        A2 = stats::qnorm(1-beta) * sqrt(var.beta1.B1)
        A3 = muT * exp.B0[i] * (B1^2) * (1-R2)

        n = phi * ((A1 + A2)^2) / A3
        foo = cbind(exp.B1,n)
        N.N1 = abs(n-N1)
        foo = cbind(exp.B1, N, N.N1)

        x = foo[which.min(foo[,3]),]

        exp.B1.new = round(x[[1]],4)

        vec.exp.B1 = c(vec.exp.B1, exp.B1.new)
      }
    } else {
      if (length(exp.B0) > 1) {
        for(i in 1:length(exp.B0)) {
          B0 = log(exp.B0[i])
          B1 = log(exp.B1)

          var.x1 = pi.x1*(1-pi.x1)
          var.beta1.0 = 1 / var.x1
          var.beta1.B1 =  (1 / (1-pi.x1)) + (1 / (pi.x1*exp.B1))

          A1 = stats::qnorm(1-alpha/2) * sqrt(var.beta1.0)
          A2 = stats::qnorm(1-beta) * sqrt(var.beta1.B1)
          A3 = muT * exp.B0[i] * (B1^2) * (1-R2)

          n = phi * ((A1 + A2)^2) / A3
          foo = cbind(exp.B1,n)
          N.N1 = abs(n-N1)
          foo = cbind(exp.B1, N, N.N1)

          x = foo[which.min(foo[,3]),]

          exp.B1.new = round(x[[1]],4)

          vec.exp.B1 = c(vec.exp.B1, exp.B1.new)
        }
      }

      if(length(N) > 1) {
        for(i in 1:length(N)) {
          B0 = log(exp.B0)
          B1 = log(exp.B1)

          var.x1 = pi.x1*(1-pi.x1)
          var.beta1.0 = 1 / var.x1
          var.beta1.B1 =  (1 / (1-pi.x1)) + (1 / (pi.x1*exp.B1))

          A1 = stats::qnorm(1-alpha/2) * sqrt(var.beta1.0)
          A2 = stats::qnorm(1-beta) * sqrt(var.beta1.B1)
          A3 = muT * exp.B0 * (B1^2) * (1-R2)

          n = phi * ((A1 + A2)^2) / A3
          foo = cbind(exp.B1, n)
          N.N1 = abs(n - N1[i])
          foo = cbind(exp.B1, N[i], N.N1)

          x = foo[which.min(foo[,3]),]

          exp.B1.new = round(x[[1]],4)

          vec.exp.B1 = c(vec.exp.B1, exp.B1.new)
        }
      }

      if(length(alpha) > 1) {
        for(i in 1:length(alpha)) {
          B0 = log(exp.B0)
          B1 = log(exp.B1)

          var.x1 = pi.x1*(1-pi.x1)
          var.beta1.0 = 1 / var.x1
          var.beta1.B1 =  (1 / (1-pi.x1)) + (1 / (pi.x1*exp.B1))

          A1 = stats::qnorm(1-alpha[i]/2) * sqrt(var.beta1.0)
          A2 = stats::qnorm(1-beta) * sqrt(var.beta1.B1)
          A3 = muT * exp.B0 * (B1^2) * (1-R2)

          n = phi * ((A1 + A2)^2) / A3
          foo = cbind(exp.B1, n)
          N.N1 = abs(n - N1)
          foo = cbind(exp.B1, N, N.N1)

          x = foo[which.min(foo[,3]),]

          exp.B1.new = round(x[[1]],4)

          vec.exp.B1 = c(vec.exp.B1, exp.B1.new)

        }
      }

      if(length(power) > 1) {
        for(i in 1:length(power)) {
          B0 = log(exp.B0)
          B1 = log(exp.B1)

          var.x1 = pi.x1*(1-pi.x1)
          var.beta1.0 = 1 / var.x1
          var.beta1.B1 =  (1 / (1-pi.x1)) + (1 / (pi.x1*exp.B1))

          A1 = stats::qnorm(1-alpha/2) * sqrt(var.beta1.0)
          A2 = stats::qnorm(1-beta[i]) * sqrt(var.beta1.B1)
          A3 = muT * exp.B0 * (B1^2) * (1-R2)

          n = phi * ((A1 + A2)^2) / A3
          foo = cbind(exp.B1, n)
          N.N1 = abs(n - N1)
          foo = cbind(exp.B1, N, N.N1)

          x = foo[which.min(foo[,3]),]

          exp.B1.new = round(x[[1]],4)

          vec.exp.B1 = c(vec.exp.B1, exp.B1.new)

        }
      }
    }

    # output results #
    results = NULL

    if(length(alpha) > 1 | length(power) > 1 | length(exp.B0) > 1 | length(pi.x1) > 1 | length(N) > 1)  {
      rownum = max(length(alpha), length(power), length(exp.B0))
      for(i in 1:rownum) {
        results = cbind(alpha, power, N, exp.B0, vec.exp.B1, pi.x1)
      }
      colnames(results) = c("alpha", "power", "N", "baseline rate (exp.B0)", "IRR (exp.B1)", "X1 proportion (pi.x1)")
    }

    else {
      results = c(alpha, power, N, exp.B0, vec.exp.B1, pi.x1)
      results = data.frame(results)
      results = t(results)
      colnames(results) = c("alpha", "power", "N", "baseline rate (exp.B0)", "IRR (exp.B1)",  "X1 proportion (pi.x1)")
    }

  }


  return(results)
}



