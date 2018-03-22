
#' Sample size and power for logistic regression with a binary covariate.
#'
#' Calculating power/sample size for simple logistic regression with a binary predictor. The function solves for one of the following: alpha, power, OR, N.
#' @param alpha type I error rate. Can range from 0 to 1 (typically 0.05)
#' @param power 1 - Pr(type II error) Can range from 0 to 1 (typically 0.80)
#' @param P0 baseline probability that Y=1
#' @param OR odds ratio (odds(Y=1|X=1) / odds(Y=1|X=0))
#' @param R proportion of sample (N) with X1=1
#' @param N sample size
#' @return Return one of the following parameters
#' @return \code{alpha}, \code{power}, \code{OR}, \code{N}
#' @note The test is a two-sided test. For one-sided tests, double the significance level. For example, you can set alpha=0.10 to obtain one-sided test at 0.05 significance level.
#' @note \code{alpha}, \code{power}, \code{P0}, \code{OR}, \code{N} can be input as either single values or vectors. Only one parameter can be input as a vector.
#' @author David Aaby <david.aaby@@northwestern.edu>
#' @references Hsieh, F.Y., Block, D.A., and Larsen, M.D. 1998. A Simple Method of Sample Size Calculation for Linear and Logistic Regression, Statistics in Medicine, Volume 17, pages 1623-1634.
#' @export Logpower.bin
#' @import stats
#' @examples
#' Logpower.bin(alpha=0.05, power=.90, P0=0.07, OR=1.5, R=.50)            #outputs N
#' Logpower.bin(alpha=0.05, power=.90, P0=0.07, OR=c(1.5,2.0), R=.50)     #outputs a vector of Ns
#' Logpower.bin(alpha=0.05, power=.90, P0=0.07, R=.50, N=3326)            #outputs OR
#' Logpower.bin(power=.90, P0=0.07, OR=1.5, R=.50, N=3326)                #outputs alpha
#' Logpower.bin(alpha=0.05, P0=0.07, OR=1.5, R=.50, N=3326)               #outputs beta



Logpower.bin <- function(alpha=NULL, power=NULL, P0=NULL, OR=NULL, R=NULL, N=NULL) {


  # error messages #
  if(any(alpha < 0 | alpha > 1)) stop('alpha not between 0 and 1')
  if(any(power < 0 | power > 1)) stop('power not between 0 and 1')
  if(any(P0 < 0 | P0 > 1))       stop('P0 not between 0 and 1')
  if(any(R < 0 | R > 1))         stop('R not between 0 and 1')
  if(any(OR < 0))                stop('OR not a positive value')
  if(any(N < 2))                 stop('N is less than 2')

  mylist = list(alpha, power, P0, OR, N)
  l.mylist = lengths(mylist)
  if(length(l.mylist[l.mylist>=2]) > 1) stop('Only vary one parameter at a time')


  # solve for alpha level #
  if(is.null(alpha)) {
    beta = 1 - power
    P1 = (P0*OR) / (OR*P0 + 1 - P0)
    Pbar = (1-R)*P0 + R*P1

    A1 = sqrt(N*((P0-P1)^2)*(1-R))
    A2 = stats::qnorm(1-beta)*sqrt(P0*(1-P0) + ((P1*(1-P1)*(1-R))/R))
    A3 = sqrt(Pbar*(1-Pbar)/R)
    A4 = (A1 - A2) / A3

    alpha = 2*(1-stats::pnorm(A4))
    alpha = round(alpha,3)

    P1 = round(P1,3)


    results = NULL

    if(length(power) > 1 | length(P0) > 1 | length(OR) > 1 | length(N) > 1) {
      rownum = max(length(alpha), length(power), length(P0), length(OR))
      for(i in 1:rownum) {
        results = cbind(alpha, power, P0, P1, OR, R, N)
      }
    }

    else {
      results = c(alpha, power, P0, P1, OR, R, N)
      results = data.frame(results)
      results = t(results)
      colnames(results) = c("alpha", "power", "P0", "P1", "OR", "R" ,"N")
    }
  }


  # solve for power #
  if(is.null(power)) {

    P1 = (P0*OR) / (OR*P0 + 1 - P0)
    Pbar = (1-R)*P0 + R*P1

    A1 = sqrt(N*((P0-P1)^2)*(1-R))
    A2 = stats::qnorm(1-alpha/2)* sqrt(Pbar*(1-Pbar)/R)
    A3 = sqrt(P0*(1-P0) + ((P1*(1-P1)*(1-R))/R))
    A4 = (A1 - A2) / A3

    power = stats::pnorm(A4)
    power = round(power,3)

    P1 = round(P1,3)


    results = NULL

    if(length(alpha) > 1 | length(P0) > 1 | length(OR) > 1 | length(N) > 1) {
      rownum = max(length(alpha), length(power), length(P0), length(OR))
      for(i in 1:rownum) {
        results = cbind(alpha, power, P0, P1, OR, R, N)
      }
    }

    else {
      results = c(alpha, power, P0, P1, OR, R, N)
      results = data.frame(results)
      results = t(results)
      colnames(results) = c("alpha", "power", "P0", "P1", "OR", "R" ,"N")
    }
  }


  # solve for N #
  if(is.null(N)) {

    beta = 1 - power


    P1 = (P0*OR) / (OR*P0 + 1 - P0)


    Pbar = (1-R)*P0 + R*P1

    A1 = stats::qnorm(1-alpha/2)* sqrt(Pbar*(1-Pbar)/R)
    A2 = stats::qnorm(1-beta)*sqrt(P0*(1-P0) + ((P1*(1-P1)*(1-R))/R))
    A3 = ((P0-P1)^2) * (1-R)

    N = ((A1 + A2)^2) / A3
    N = ceiling(N)

    P1 = round(P1,3)


    results = NULL

    if(length(alpha) > 1 | length(power) > 1 | length(P0) > 1 | length(OR) > 1) {
      rownum = max(length(alpha), length(power), length(P0), length(OR))
      for(i in 1:rownum) {
        results = cbind(alpha, power, P0, P1, OR, R, N)
      }
    }

    else {
      results = c(alpha, power, P0, P1, OR, R, N)
      results = data.frame(results)
      results = t(results)
      colnames(results) = c("alpha", "power", "P0", "P1", "OR", "R" ,"N")
    }
  }

  # solve for odds ratio #
  if(is.null(OR)) {

    beta = 1 - power
    P1 = seq(0, .90, .00001)
    N1 = N
    OR = NULL
    P1.new = NULL

    maxlength = max(c(length(alpha), length(power), length(P0), length(N)))
    if(maxlength==1) {

      for(i in 1:length(P0)) {
        Pbar = (1-R)*P0[i] + R*P1
        n = (((stats::qnorm(1-alpha/2)* sqrt(Pbar*(1-Pbar)/R)) +
              (stats::qnorm(1-beta)*sqrt(P0[i]*(1-P0[i]) + ((P1*(1-P1)*(1-R))/R))))^2) /
        (((P0[i]-P1)^2) * (1-R))

        or = (P1*(1-P0[i])) / (P0[i]*(1-P1))
        N.N1 = abs(n - N1)
        foo = cbind(P1,N,or, N.N1)
        foo = foo[which(foo[,3] >= 1),]
        x = foo[which.min(foo[,4]),]

        or = x[[3]]
        OR = c(OR, or)
        OR = round(OR,3)

        p1 = x[[1]]
        P1.new = c(P1.new, p1)
        P1.new = round(P1.new,3)
      }
    } else {
      if (length(P0) > 1) {
        for(i in 1:length(P0)) {
          Pbar = (1-R)*P0[i] + R*P1
          n = (((stats::qnorm(1-alpha/2)* sqrt(Pbar*(1-Pbar)/R)) +
                (stats::qnorm(1-beta)*sqrt(P0[i]*(1-P0[i]) + ((P1*(1-P1)*(1-R))/R))))^2) /
          (((P0[i]-P1)^2) * (1-R))

          or = (P1*(1-P0[i])) / (P0[i]*(1-P1))
          N.N1 = abs(n - N1)
          foo = cbind(P1,N,or, N.N1)
          foo = foo[which(foo[,3] >= 1),]
          x = foo[which.min(foo[,4]),]

          or = x[[3]]
          OR = c(OR, or)
          OR = round(OR,3)

          p1 = x[[1]]
          P1.new = c(P1.new, p1)
          P1.new = round(P1.new,3)
        }
      }

      if(length(N) > 1) {
        for(i in 1:length(N)) {
          Pbar = (1-R)*P0 + R*P1
          n = (((stats::qnorm(1-alpha/2)* sqrt(Pbar*(1-Pbar)/R)) +
                  (stats::qnorm(1-beta)*sqrt(P0*(1-P0) + ((P1*(1-P1)*(1-R))/R))))^2) /
            (((P0-P1)^2) * (1-R))

          or = (P1*(1-P0)) / (P0*(1-P1))
          N.N1 = abs(n - N1[i])
          foo = cbind(P1, N[i], or, N.N1)
          foo = foo[which(foo[,3] >= 1),]
          x = foo[which.min(foo[,4]),]

          or = x[[3]]
          OR = c(OR, or)
          OR = round(OR,3)

          p1 = x[[1]]
          P1.new = c(P1.new, p1)
          P1.new = round(P1.new,3)
        }
      }

      if(length(alpha) > 1) {
        for(i in 1:length(alpha)) {
          Pbar = (1-R)*P0 + R*P1
          n = (((stats::qnorm(1-alpha[i]/2)* sqrt(Pbar*(1-Pbar)/R)) +
                  (stats::qnorm(1-beta)*sqrt(P0*(1-P0) + ((P1*(1-P1)*(1-R))/R))))^2) /
            (((P0-P1)^2) * (1-R))

          or = (P1*(1-P0)) / (P0*(1-P1))
          N.N1 = abs(n - N1)
          foo = cbind(P1, N, or, N.N1)
          foo = foo[which(foo[,3] >= 1),]
          x = foo[which.min(foo[,4]),]

          or = x[[3]]
          OR = c(OR, or)
          OR = round(OR,3)

          p1 = x[[1]]
          P1.new = c(P1.new, p1)
          P1.new = round(P1.new,3)
        }
      }

      if(length(power) > 1) {
        for(i in 1:length(power)) {
          Pbar = (1-R)*P0 + R*P1
          n = (((stats::qnorm(1-alpha/2)* sqrt(Pbar*(1-Pbar)/R)) +
                  (stats::qnorm(1-beta[i])*sqrt(P0*(1-P0) + ((P1*(1-P1)*(1-R))/R))))^2) /
            (((P0-P1)^2) * (1-R))

          or = (P1*(1-P0)) / (P0*(1-P1))
          N.N1 = abs(n - N1)
          foo = cbind(P1, N, or, N.N1)
          foo = foo[which(foo[,3] >= 1),]
          x = foo[which.min(foo[,4]),]

          or = x[[3]]
          OR = c(OR, or)
          OR = round(OR,3)

          p1 = x[[1]]
          P1.new = c(P1.new, p1)
          P1.new = round(P1.new,3)
        }
      }

    }

    results = NULL

    if(length(alpha) > 1 | length(power) > 1 | length(P0) > 1 | length(N) > 1) {
      rownum = max(length(alpha), length(power), length(P0), length(OR))
      for(i in 1:rownum) {
        results = cbind(alpha, power, P0, P1.new, OR, R, N)
      }
      colnames(results) = c("alpha", "power", "P0", "P1", "OR", "R" ,"N")
    }

    else {
      results = c(alpha, power, P0, P1.new, OR, R, N)
      results = data.frame(results)
      results = t(results)
      colnames(results) = c("alpha", "power", "P0", "P1", "OR", "R" ,"N")
    }
  }

  return(results)
}



