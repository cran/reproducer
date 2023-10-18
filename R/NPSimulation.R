#' @title checkIfValidDummyVariable
#' @description This helper function checks whether a vector variable comprises
#' only zeros and 1's.
#' @author Barbara Kitchenham and Lech Madeyski
#' @param vector Vector variable
#' @return The logical value:
#' TRUE - if a vector variable passed as the function's parameter represents a
#' valid dummy variable, i.e., comprises only zeros and 1's.
#' FALSE - if a vector variable passed as the function's parameter does not
#' represent a valid dummy variable.
#' @examples
#' print(reproducer:::checkIfValidDummyVariable(c(0, 1, 0, 0)))
#' # [1] TRUE
#' print(reproducer:::checkIfValidDummyVariable(c(0, 1, 2, 0)))
#' # [1] FALSE
checkIfValidDummyVariable <- function(vector) {
  result <- all(as.logical(levels(factor(vector)) %in% c(0, 1)))
  return(result)
}

#' @title CatchError
#' @description This is a helper function to stop simulations failing if the
#' metafor function rma fails for example cannot converge properly for
#' a specific dataset.
#' @author Barbara Kitchenham and Lech Madeyski
#' @param expr The expression that is being monitored
#' @return A message confirming whether the expression has performed
#' successfully
#' @examples
#' ES=c(0.2,0.3)
#' ESvar=c(0.04,0.03)
#' outcome=reproducer:::CatchError(metafor::rma(ES,ESvar,method='Meth'))
#' outcome
#' # [1] 'Failure'
CatchError <- function(expr) {
  tryCatch(
    expr,
    error = function(e) {
      "Failure"
    },
    warning = function(w) {
      "Failure"
    },
    finally = "Success"
  )
}

#' @title testfunctionParameterChecks
#' @description This is a helper function that ensures parameter values used
#' for performing special statistical tests are valid.
#' @author Barbara Kitchenham and Lech Madeyski
#' @param alternative The type of statistical test. Valid values are one of
#' c('two.sided', 'greater', 'less')
#' @param alpha The test level. Valid values are between 0.0001 and 0.2
#' @param stderr The standard error of a parameter whose confidence intervals
#' is to be calculated
#' @return 'Success' or an error message.
#' @examples
#' #reproducer:::testfunctionParameterChecks(alternative='larger',alpha=0.1,stderr=0.002)
#' #Error in testfunctionParameterChecks(alternative = 'larger', alpha = 0.1) :
#' #  Invalid alternative parameter, choose one of two.sided, greater or less
#' reproducer:::testfunctionParameterChecks(alternative='greater',alpha=0.1,stderr=0.002)
#' #[1] 'Success'
#' #reproducer:::testfunctionParameterChecks(alternative='greater',alpha=0.1,stderr=0.000)
#' #Error in testfunctionParameterChecks(alternative = 'greater', alpha = 0.1,  :
#' #  Improbably small variance, data are essentially constant
testfunctionParameterChecks <-
  function(alternative, alpha, stderr) {
    alternative.allowed <- c("two.sided", "greater", "less")
    Alt.test <-
      CatchError(match.arg(arg = alternative, choices = alternative.allowed))
    if (Alt.test == "Failure") {
      stop("Invalid alternative parameter, choose one of two.sided, greater or less")
    }
    if (alpha < 1e-04 | alpha > 0.2) {
      stop("Invalid alpha parameter, select alpha in range (0.0001,0.2)")
    }
    if (stderr < 10 * .Machine$double.eps) {
      stop("Improbably small variance, data are essentially constant")
    }
    return("Success")
  }

#' @title calcEffectSizeConfidenceIntervals
#' @description This function provides single-sided and two-sided confidence interval of an effect size (assuming that the null hypothesis value is zero).
#' @author Barbara Kitchenham and Lech Madeyski
#' @param effectsize The effect size
#' @param effectsize.variance The effect size variance
#' @param effectsize.df The degrees of freedom for confidence intervals based on the t- distribution. If df=0 (default), the confidence interval is based on the normal distribution
#' @param alpha The significance level of the confidence interval
#' (default 0.05).
#' @param alternative This defines whether a one-sided test or a two-sided
#' (default) test is required. For a one-sided test use parameter values
#' greater' or 'less' to define whether the d-value should be greater or less
#' than zero.
#' @param UpperValue The maximum legal value of the effect size (default Inf).
#' Used to ensure that confidence intervals of effect sizes such as correlation
#' coefficients are restricted to sensible values.
#' @param LowerValue The minimum legal value of the effect size (default -Inf).
#' Used to ensure that confidence intervals of effect sizes
#' such as correlation coefficients are restricted to sensible values
#' @return The value of the test statistic, the p.value of test statistic, the upper and lower confidence interval of the effect size, a logical value specifying whether the effect size is significantly different from zero based on the confidence interval and the lower and upper confidence interval bounds.
#' @examples
#' reproducer:::calcEffectSizeConfidenceIntervals(
#'   effectsize=0.37,effectsize.variance=0.00847,effectsize.df=11.1,
#'   alpha=0.05,alternative='two.sided',UpperValue=0.5,LowerValue=-0.5)
#' # A tibble: 1 x 5
#' #  ES.test ES.pvalue ES.sig ES.ci.lower ES.ci.upper
#' #    <dbl>     <dbl> <lgl>        <dbl>       <dbl>
#' #1    4.02   0.00198 TRUE         0.168         0.5
calcEffectSizeConfidenceIntervals <-
  function(effectsize,
           effectsize.variance,
           effectsize.df = 0,
           alpha = 0.05,
           alternative = "two.sided",
           UpperValue = Inf,
           LowerValue = -Inf) {
    testfunctionParameterChecks(
      alternative = alternative,
      alpha = alpha,
      stderr = sqrt(effectsize.variance)
    )

    ES.se <- sqrt(effectsize.variance)
    ES.test <- effectsize / ES.se
    useTTest <- (effectsize.df > 0)


    if (alternative == "two.sided") {
      if (useTTest) {
        vv <- stats::qt(alpha / 2, effectsize.df)
        ES.pvalue <- 2 * (1 - stats::pt(abs(ES.test), effectsize.df))
      } else {
        vv <- stats::qnorm(alpha / 2)
        ES.pvalue <- 2 * (1 - stats::pnorm(abs(ES.test)))
      }
      ES.ci.lower <- effectsize + vv * ES.se
      ES.ci.upper <- effectsize - vv * ES.se
      ES.sig <- (ES.ci.upper < 0 | ES.ci.lower > 0)
    } else {
      if (useTTest) {
        vv <- stats::qt(alpha, effectsize.df)
      } else {
        vv <- stats::qnorm(alpha)
      }

      if (alternative == "greater") {
        ES.ci.lower <- effectsize + vv * ES.se
        ES.ci.upper <- UpperValue
        ES.sig <- ES.ci.lower > 0
        if (useTTest) {
          ES.pvalue <- (1 - stats::pt(ES.test, effectsize.df))
        } else {
          ES.pvalue <- (1 - stats::pnorm(ES.test))
        }
      } else {
        ES.ci.upper <- effectsize - vv * ES.se
        ES.ci.lower <- LowerValue
        ES.sig <- ES.ci.upper < 0

        if (useTTest) {
          ES.pvalue <- (stats::pt(ES.test, effectsize.df))
        } else {
          ES.pvalue <- (stats::pnorm(ES.test))
        }
      }
    }

    if (ES.ci.upper > UpperValue) {
      ES.ci.upper <- UpperValue
    }
    if (ES.ci.lower < LowerValue) {
      ES.ci.lower <- LowerValue
    }
    out <-
      tibble::tibble(
        ES.test = ES.test,
        ES.pvalue = ES.pvalue,
        ES.sig = ES.sig,
        ES.ci.lower = ES.ci.lower,
        ES.ci.upper = ES.ci.upper
      )

    return(out)
  }


#' @title calculateCliffd
#' @description This function implements finds Cliff's d and its confidence intervals. The null hypothesis is that for two independent group, P(X<Y)=P(X>Y). The function reports a 1-alpha confidence interval for P(X>Y)-P(X<Y). The algorithm computes a confidence interval for Cliff's d using the method in Cliff, 1996, p. 140, eq 5.12. The function is based on code produce by Rand Wilcox but has been amended. The plotting function has been removed and the dependency on Wilcox's binomci function has been removed. Construction of confidence intervals if values in one group are all larger than values in the other group has been amended to use the smallest non-zero variance method. Upper and lower confidence interval bounds cannot assume invalid values, i.e. values <-1 or >1.
#' @author Rand Wilcox, amendments Barbara Kitchenham and Lech Madeyski
#' @export calculateCliffd
#' @param x is a vector of values from group 1
#' @param y is a vector of values from group 2
#' @param alpha is the Type 1 error level for statistical tests
#' @param sigfig is the number of significant digit. If sigfig>0 the data in x and y is truncated to the specified value.
#' @return list including the value of Cliffs d its consistent variance and confidence intervals and the equivalent probability of superiority value and its confidence intervals.
#' @examples
#' x=c(1.2,3,2.2,4,2.5,3)
#' y=c(3,4.2,4,6,7,5.9)
#' calculateCliffd(x,y)
#' #  $n1
#' # [1] 6
#' # $n2
#' # [1] 6
#' # $d
#' # [1] -0.8611111
#' # $sqse.d
#' # [1] 0.02017931
#' # $phat
#' # [1] 0.06944444

#' z=c(1,2,3,4)
#' y=c(5,6,7,8)
#' calculateCliffd(z,y)
#' # $n1
#' # [1] 4
#' # $n2
#' # [1] 4
#' # $d
#' # [1] -1
#' # $sqse.d
#' # [1] 0.009765625
#' # $phat
#' # [1] 0


calculateCliffd <- function(x,
                            y,
                            alpha = 0.05,
                            sigfig = -1) {
  # Check that the data is valid
  if (length(x) <= 1) {
    stop("Too few data points")
  }
  if (length(y) <= 1) {
    stop("Too few data points")
  }
  if (length(x) != length(x[!is.na(x)])) {
    stop("Missing values not permitted")
  }
  if (length(y) != length(y[!is.na(y)])) {
    stop("Missing values not permitted")
  }

  # Truncate the data if necessary
  if (sigfig > 0) {
    x <- signif(x, sigfig)
    y <- signif(y, sigfig)
  }

  m <- base::outer(x, y, FUN = "-")
  msave <- m
  m <- sign(m)

  d <- base::mean(m)

  phat <- (1 + d) / 2

  flag <- TRUE
  if (d == -1 || d == 1) {
    flag <- FALSE
  }
  q0 <- sum(msave == 0) / length(msave)
  qxly <- sum(msave < 0) / length(msave)
  qxgy <- sum(msave > 0) / length(msave)
  c.sum <- base::matrix(c(qxly, q0, qxgy), nrow = 1, ncol = 3)
  dimnames(c.sum) <- list(NULL, c("P(X<Y)", "P(X=Y)", "P(X>Y)"))
  if (flag) {
    # This is appropriate for the consistent variance of d
    sigdih <- sum((m - d) ^ 2) / (length(x) * length(y) - 1)

    di <- NA
    for (i in 1:length(x))
      di[i] <- sum(x[i] > y) / length(y) - sum(x[i] < y) / length(y)

    dh <- NA
    for (i in 1:length(y))
      dh[i] <- sum(y[i] < x) / length(x) - sum(y[i] > x) / length(x)


    sdi <- stats::var(di)
    sdh <- stats::var(dh)
    # sh is the consistent variance of d
    sh <-
      ((length(y) - 1) * sdi + (length(x) - 1) * sdh + sigdih) / (length(x) * length(y))
  }
  if (!flag) {
    # phat=0 and d=-1 or phat=1 and d=1 Cannot calculate standard error of d, use alternative minimum change approach to calculate a 'close'
    # slightly too large estimate. Calculates the smallest non-zero consistent variance
    sh <- 0
    di <- NA
    dh <- NA
    nx <- length(x)
    ny <- length(y)

    tempn <- nx + ny - 1

    if (d == -1) {
      tempn <- nx + ny - 1
      disturbedx <- c(nx:tempn)
      disturbedy <- c(1:nx)
      disturbedanalysis <- calculateCliffd(disturbedx, disturbedy)
    }

    if (d == 1) {
      disturbedy <- c(nx:tempn)
      disturbedx <- c(1:nx)
      disturbedanalysis <- calculateCliffd(disturbedx, disturbedy)
    }

    sh <- disturbedanalysis$sqse
  }

  cid.results <-
    list(
      n1 = length(x),
      n2 = length(y),
      d = d,
      sqse.d = sh,
      phat = phat
    )
  return(cid.results)
}


#' @title calculatePhat
#' @description This function calculates the probability of superiority (i.e., Phat) and its confidence interval based on Brunner and Munzel (2000) heteroscedastic analog of WMW test. It is based on Wilcox'x bmp function with some amendments. It does not include a plotit facility. It uses the smallest non-zero variance to identify confidence intervals and statistical significance for values of Phat=0 and Phat=1. It ensure that confidence intervals do not take on invalid values such as values <0 or >1.
#' @author Rand Wilcox amendments by Barbara Kitchenham and Lech Madeyski
#' @export calculatePhat
#' @param x is a vector of values from group 1
#' @param y is a vector of values from group 2
#' @param alpha is the Type 1 error level for statistical tests
#' @param sigfig is the number of significant digits. If sigfig>0 the data in x
#' and y is truncated to the specified number of significant digits.
#' @return list including the value of the t-test for PHat, the estimate of PHat and Cliff's d, and the confidence intervals for PHat.
#' @examples
#' x <- c(1.2, 3.0, 2.2, 4.0, 2.5, 3.0)
#' y <- c(3, 4.2, 4, 6, 7, 5.9)
#' reproducer:::calculatePhat(x, y)
#' # $test.stat
#' # [1] 6.381249
#' # $phat
#' # [1] 0.9305556
#' # $dhat
#' # [1] 0.8611111
#' # $sig.level
#' # [1] 0.0001191725
#' # $s.e.
#' # [1] 0.06747199
#' # $ci.p
#' # [1] 0.7783001 1.0000000
#' # $df
#' # [1] 9.148489
#' # Another example:
#' z <- c(1, 2, 3, 4)
#' y <- c(5, 6, 7, 8)
#' reproducer:::calculatePhat(z, y)
#' # $test.stat
#' # [1] 10.6066
#' # $phat
#' # [1] 1
#' # $dhat
#' # [1] 1
#' # $sig.level
#' # [1] 4.135921e-05
#' # $s.e.
#' # [1] 0.04419417
#' # $ci.p
#' # [1] 0.8918608 1.0000000
#' # $df
#' # [1] 6
calculatePhat <- function(x,
                          y,
                          alpha = 0.05,
                          sigfig = -1) {
  if (sigfig > 0) {
    x <- signif(x, sigfig)
    y <- signif(y, sigfig)
  }

  x <- x[!is.na(x)] # Remove any missing values
  y <- y[!is.na(y)]
  n1 <- length(x)
  n2 <- length(y)
  N <- n1 + n2
  n1p1 <- n1 + 1
  flag1 <- c(1:n1)
  flag2 <- c(n1p1:N)
  R <- rank(c(x, y))
  R1 <- mean(R[flag1])
  R2 <- mean(R[flag2])
  phat <- (R2 - (n2 + 1) / 2) / n1
  dhat <- 2 * phat - 1

  flag <- TRUE # TRUE if phat ne 0 and phat ne 1
  if (phat == 0 | phat == 1) {
    flag <- FALSE
  }

  if (flag) {
    Rg1 <- base::rank(x)
    Rg2 <- base::rank(y)
    S1sq <- sum((R[flag1] - Rg1 - R1 + (n1 + 1) / 2) ^ 2) / (n1 - 1)
    S2sq <- sum((R[flag2] - Rg2 - R2 + (n2 + 1) / 2) ^ 2) / (n2 - 1)
    sig1 <- S1sq / n2 ^ 2
    sig2 <- S2sq / n1 ^ 2
    se <- sqrt((sig1 / n1 + sig2 / n2))

    df <-
      (S1sq / n2 + S2sq / n1) ^ 2 / ((S1sq / n2) ^ 2 / (n1 - 1) + (S2sq / n1) ^
                                       2 / (n2 - 1))
  } else {
    # Calculate the smallest non-negative variance and use results to approximate variance and confidence intervals of phat. Gives a more
    # realistic confidence interval than other methods
    Nl1 <- N - 1
    newx <- c(1:n1)
    newy <- c(n1:Nl1)
    newres <- calculatePhat(newx, newy)
    se <- newres$s.e.
    df <- newres$df
  }


  list(
    phat = phat,
    dhat = dhat,
    s.e. = se,
    df = df
  )
}


#' @title calcCliffdConfidenceIntervals
#' @description This functions is a helper function. It assesses the significance one-sided and two-sided statistical of Cliff's d based on its confidence interval. The type of test is determined by the parameter One.Sided.Tests, the direction of one-sided tests is determined by the parameter Positive.MD.
#' @author Barbara Kitchenham and Lech Madeyski
#' @param d.value This is the value of Cliff's d.
#' @param d.variance This is the estimated variance of Cliff's d
#' @param d.df The degrees of freedom.
#' @param alpha This is the alpha level required for the statistical tests (default 0.05)
#' @param alternative This defines whether a one-sided test or a two-sided
#' (default) test is required. For a one-sided test use parameter values
#' greater' or 'less' to define whether the d-value should be greater or less
#' than zero.
#' @return The function returns a Boolean variable identifying whether the effect size is significant and the confidence interval bounds.
#' @examples
#' reproducer:::calcCliffdConfidenceIntervals(d.value=0.5, d.variance=0.04,d.df=18)
#' # A tibble: 1 x 5
#' #  d.tvalue d.pvalue d.ci.lower d.ci.upper d.sig
#' #     <dbl>    <dbl>      <dbl>      <dbl> <lgl>
#' #1      2.5   0.0223     0.0479      0.782 TRUE
#'
#' reproducer:::calcCliffdConfidenceIntervals(
#'   d.value=0.5,d.variance=0.04,d.df=18,alternative='greater')
#' # A tibble: 1 x 5
#' #  d.tvalue d.pvalue d.ci.lower d.ci.upper d.sig
#' #   <dbl>    <dbl>      <dbl>      <dbl> <lgl>
#' #1    2.5   0.0112      0.123          1 TRUE
#'
#' reproducer:::calcCliffdConfidenceIntervals(
#'   d.value=0.2,d.variance=0.04,d.df=18,alternative='greater')
#' # A tibble: 1 x 3
#' #  d.tvalue d.pvalue d.ci.lower d.ci.upper d.sig
#' #     <dbl>    <dbl>      <dbl>      <dbl> <lgl>
#' #1        1    0.165     -0.133          1 FALSE
#'
#' reproducer:::calcCliffdConfidenceIntervals(
#'   d.value=-0.5,d.variance=0.04,d.df=18,alternative='less')
#' # A tibble: 1 x 5
#' #  d.tvalue d.pvalue d.ci.lower d.ci.upper d.sig
#' #   <dbl>    <dbl>      <dbl>      <dbl> <lgl>
#' #1     -2.5   0.0112         -1     -0.123 TRUE
calcCliffdConfidenceIntervals <-
  function(d.value,
           d.variance,
           d.df,
           alpha = 0.05,
           alternative = "two.sided") {
    d <- d.value
    vard <- d.variance
    d.tvalue <- d / sqrt(vard)

    if (alternative == "two.sided") {
      zv <- stats::qnorm(alpha / 2)
      d.cu <-
        (d - d ^ 3 - zv * sqrt(vard) * sqrt((1 - d ^ 2) ^ 2 + zv ^ 2 * vard)) / (1 - d ^
                                                                                   2 + zv ^ 2 * vard)
      d.cl <-
        (d - d ^ 3 + zv * sqrt(vard) * sqrt((1 - d ^ 2) ^ 2 + zv ^ 2 * vard)) / (1 - d ^
                                                                                   2 + zv ^ 2 * vard)
      d.sig <- d.cu < 0 | d.cl > 0

      d.pvalue <- 2 * (1 - stats::pt(abs(d.tvalue), d.df))
    } else {
      zv <- stats::qnorm(alpha)
      if (alternative == "greater") {
        d.cu <- 1
        d.cl <-
          (d - d ^ 3 + zv * sqrt(vard) * sqrt((1 - d ^ 2) ^ 2 + zv ^ 2 * vard)) / (1 - d ^
                                                                                     2 + zv ^ 2 * vard)
        d.sig <- d.cl > 0
        d.pvalue <- 1 - stats::pt(d.tvalue, d.df)
      } else {
        d.cl <- -1
        d.cu <-
          (d - d ^ 3 - zv * sqrt(vard) * sqrt((1 - d ^ 2) ^ 2 + zv ^ 2 * vard)) / (1 - d ^
                                                                                     2 + zv ^ 2 * vard)

        d.sig <- d.cu < 0
        d.pvalue <- stats::pt(d.tvalue, d.df)
      }
    }

    if (d.cl < -1) {
      d.cl <- -1
    }
    if (d.cu > 1) {
      d.cu <- 1
    }

    out <-
      tibble::tibble(
        d.tvalue = d.tvalue,
        d.pvalue = d.pvalue,
        d.ci.lower = d.cl,
        d.ci.upper = d.cu,
        d.sig = d.sig
      )

    return(out)
  }

#' @title Cliffd.test
#' @description This function provides single-sided and two-sided tests of Cliff's d
#' @author Barbara Kitchenham and Lech Madeyski
#' @export Cliffd.test
#' @param x The data from one group
#' @param y The data from the alternative group
#' @param alpha The significance level of tests which also controls
#' the values of the confidence interval (default 0.05)
#' @param alternative This defines whether a one-sided test or a two-sided test
#' is required (default "two.sided"). For a one-sided test use parameter values
#' 'greater' or 'less' to define whether the d-value should be greater or
#' less than zero.
#' @param sigfig is the number of significant digits. If sigfig>0 the data in x
#' and y is truncated to the specified number of significant digits.
#' @return The values of Cliff's d and its standard error, the t-value,
#' its pvalue and the upper and lower confidence interval.
#' @examples
#' a=c(1.2,3,2.2,4,2.5,3)
#' b=c(3,4.2,4,6,7,5.9)
#' Cliffd.test(a,b,alpha = .05,alternative='two.sided',sigfig = -1)
#' # A tibble: 1 x 7
#' #       d sqse.d d.tvalue d.pvalue d.ci.lower d.ci.upper d.sig
#' #   <dbl>  <dbl>    <dbl>    <dbl>      <dbl>      <dbl> <lgl>
#' #1 -0.861 0.0202    -42.7 1.20e-12     -0.896     -0.816 TRUE
#'
#' Cliffd.test(b,a,alpha = .05,alternative='greater',sigfig = -1)
#' # A tibble: 1 x 7
#' #      d sqse.d d.tvalue d.pvalue d.ci.lower d.ci.upper d.sig
#' #  <dbl>  <dbl>    <dbl>    <dbl>      <dbl>      <dbl> <lgl>
#' #1 0.861 0.0202     42.7 5.99e-13      0.824          1 TRUE
#'
Cliffd.test <-
  function(x,
           y,
           alpha = 0.05,
           alternative = "two.sided",
           sigfig = -1) {
    Cliffd.stats <-
      calculateCliffd(x, y, alpha = alpha, sigfig = sigfig)


    d <- Cliffd.stats$d
    sqse.d <- Cliffd.stats$sqse.d # SQuared Standard Error
    df <- length(x) + length(y) - 2

    testfunctionParameterChecks(alternative = alternative,
                                alpha = alpha,
                                stderr = sqrt(sqse.d))

    CliffsCI <-
      calcCliffdConfidenceIntervals(
        d.value = d,
        d.variance = sqse.d,
        d.df = df,
        alpha = alpha,
        alternative = alternative
      )
    out <- tibble::tibble(
      d = d,
      sqse.d = sqse.d,
      d.tvalue = CliffsCI$d.tvalue,
      d.pvalue = CliffsCI$d.pvalue,
      d.ci.lower = CliffsCI$d.ci.lower,
      d.ci.upper = CliffsCI$d.ci.upper,
      d.sig = CliffsCI$d.sig
    )
    return(out)
  }



#' @title calcPHatConfidenceIntervals
#' @description This functions is a helper function. It assesses the significance one-sided and two-sided statistical of the probability of superiority based on its confidence interval. The type of test and the direction of the test is determined by the parameter alternative which takes one of the values 'two.sided', 'greater' or 'less'.
#' @author Barbara Kitchenham and Lech Madeyski
#' @param phat This is the value of the probability of superiority.
#' @param phat.variance This is the estimated variance of the probability of superiority.
#' @param phat.df The degrees of freedom associated with phat value.
#' @param alpha This is the alpha level required for the statistical tests
#' (default 0.05)
#' @param alternative This defines whether a one-sided test or a two-sided
#' (i.e. default) test is required. For a one-sided test use parameter values
#' greater' or 'less' to define whether the phat-value should be greater or
#' less than 0.5.
#' @return The function returns a Boolean variable identifying whether the effect size is significant and the confidence interval bounds.
#' @examples
#' reproducer:::calcPHatConfidenceIntervals(.65,0.005,8)
#' # A tibble: 1 x 5
#' #  phat.test pvalue phat.sig phat.ci.lower phat.ci.upper
#' #     <dbl>  <dbl> <lgl>            <dbl>         <dbl>
#' # 1      2.12 0.0667 FALSE            0.487         0.813
#' reproducer:::calcPHatConfidenceIntervals(.65,0.005,8,alternative='greater')
#' # A tibble: 1 x 5
#' #  phat.test pvalue phat.sig phat.ci.lower phat.ci.upper
#' #      <dbl>  <dbl> <lgl>            <dbl>         <dbl>
#' # 1      2.12 0.0333 TRUE             0.519             1

calcPHatConfidenceIntervals <-
  function(phat,
           phat.variance,
           phat.df,
           alpha = 0.05,
           alternative = "two.sided") {
    phat.se <- sqrt(phat.variance)
    phat.test <- (phat - 0.5) / phat.se


    if (alternative == "two.sided") {
      vv <- stats::qt(alpha / 2, phat.df)
      phat.ci.lower <- phat + vv * phat.se
      phat.ci.upper <- phat - vv * phat.se
      phat.sig <- phat.ci.upper < 0.5 | phat.ci.lower > 0.5
      phat.pvalue <- 2 * (1 - stats::pt(abs(phat.test), phat.df))
    } else {
      vv <- stats::qt(alpha, phat.df)
      if (alternative == "greater") {
        phat.ci.lower <- phat + vv * phat.se
        phat.ci.upper <- 1
        phat.sig <- phat.ci.lower > 0.5
        phat.pvalue <- 1 - stats::pt(phat.test, phat.df)
      } else {
        phat.ci.upper <- phat - vv * phat.se
        phat.ci.lower <- 0
        phat.sig <- phat.ci.upper < 0.5
        phat.pvalue <- stats::pt(phat.test, phat.df)
      }
    }

    if (phat.ci.upper > 1) {
      phat.ci.upper <- 1
    }
    if (phat.ci.lower < 0) {
      phat.ci.lower <- 0
    }
    out <-
      tibble::tibble(
        phat.test = phat.test,
        pvalue = phat.pvalue,
        phat.sig = phat.sig,
        phat.ci.lower = phat.ci.lower,
        phat.ci.upper = phat.ci.upper
      )

    return(out)
  }

#' @title PHat.test
#' @description This function provides single-sided and two-sided tests of the probability of superiority (phat).
#' @author Barbara Kitchenham and Lech Madeyski
#' @export PHat.test
#' @param x The data from one group
#' @param y The data from the alternative group
#' @param alpha The significance level of tests which also controls the values
#' of the confidence interval (default 0.05)
#' @param alternative This defines whether a one-sided test or a two-sided
#' (default) test is required. For a one-sided test use parameter values
#' greater' or 'less' to define whether the d-value should be greater or less
#' than zero.
#' @param sigfig is the number of significant digits in the data.
#' @return The values of phat and its standard error,the t-value, its pvalue and the upper and lower confidence interval.
#' @examples
#' set.seed(456)
#' x <- rnorm(10, 0, 1)
#' y <- rnorm(10, 0.8, 1)
#' PHat.test(x, y, alpha = .05, alternative = "greater", sigfig = -1)
#' # A tibble: 1 x 8
#' #    phat sqse.phat phat.df phat.tvalue phat.pvalue phat.ci.lower phat.ci.upper phat.sig
#' #   <dbl>     <dbl>   <dbl>       <dbl>       <dbl>         <dbl>         <dbl> <lgl>
#' # 1  0.79    0.0118    13.6        2.67     0.00924         0.599             1 TRUE
#' PHat.test(x, y, alpha = .05, alternative = "two.sided", sigfig = -1)
#' # A tibble: 1 x 8
#' # phat sqse.phat phat.df phat.tvalue phat.pvalue phat.ci.lower phat.ci.upper phat.sig
#' #  <dbl>     <dbl>   <dbl>       <dbl>       <dbl>         <dbl>         <dbl> <lgl>
#' # 1  0.79    0.0118    13.6        2.67      0.0185         0.557             1 TRUE
#'
PHat.test <-
  function(x,
           y,
           alpha = 0.05,
           alternative = "two.sided",
           sigfig = -1) {
    PHat.stats <- calculatePhat(x, y, alpha = alpha, sigfig = sigfig)

    phat <- PHat.stats$phat
    sqse.phat <- PHat.stats$s.e ^ 2
    df <- PHat.stats$df

    testfunctionParameterChecks(
      alternative = alternative,
      alpha = alpha,
      stderr = sqrt(sqse.phat)
    )


    PhatCI <-
      calcPHatConfidenceIntervals(
        phat = phat,
        phat.variance = sqse.phat,
        phat.df = df,
        alpha = alpha,
        alternative = alternative
      )
    out <- tibble::tibble(
      phat = phat,
      sqse.phat = sqse.phat,
      phat.df = df,
      phat.tvalue = PhatCI$phat.test,
      phat.pvalue = PhatCI$pvalue,
      phat.ci.lower = PhatCI$phat.ci.lower,
      phat.ci.upper = PhatCI$phat.ci.upper,
      phat.sig = PhatCI$phat.sig
    )
    return(out)
  }

#' @title calc.a
#' @description This function is a helper function that calculates one element of the standardized mean difference effect size variance based on Hedges and Olkin p128-131.
#' @author Barbara Kitchenham and Lech Madeyski
#' @param f a vector defining the degrees of freedom for the effect sizes
#' @param A a vector defining the constants that relate each StdMD to its related t-variable where t^2=Ad^2
#' @return The value of a.
#' @examples
#' reproducer:::calc.a(10,2/10)
#' # [1] 0.2128649
calc.a <- function(f, A) {
  c <- reproducer::calculateSmallSampleSizeAdjustment(f)
  a <- c ^ 2 * (f) * A / (f - 2)
  return(a)
}

#' @title calc.b
#' @description This function is a helper function that calculates one element of the standardized mean difference effect size variance based on Hedges and Olkin p128-131.
#' @author Barbara Kitchenham and Lech Madeyski
#' @param f a vector defining the degrees of freedom for the effect sizes
#' @return The value of b.
#' @examples
#' reproducer:::calc.b(8)
#' #0.08649774
#'
calc.b <- function(f) {
  c <- reproducer::calculateSmallSampleSizeAdjustment(f)
  b <- c ^ 2 * (f) / (f - 2) - 1
  return(b)
}

#' @title metaanalyseSmallSampleSizeExperiments
#' @description Implements analysis of small sample size experiments based on Hedges and Olkin p128-131.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export metaanalyseSmallSampleSizeExperiments
#' @param d a vector of standardized mean differences for different experiments (not adjusted for small sample size)
#' @param f a vector defining the degrees  for each experiment
#' @param A a vector defining the relationship between d and its related t value for each experiment
#' @return UnweightedMean The unweighted mean of the small size adjusted standardized mean differences
#' @return WeightedMean The weighted mean of the small size adjusted standardized mean differences
#' @return VarWeightedMean The variance of the weighted mean
#' @examples
#' d=c(.461,.782,.513,.612,-0.131,-0.018,0.774,0.138,-0.482,0.333,.701,-0.222,
#'   .399,.538,-0.19,0.833,0.512,0.601,-0.366,.510)
#' A=c(2/5,2/5,2/7,2/10,2/8,2/10,2/6,2/6,2/4,2/9,2/9,2/7,2/7,2/5,2/5,2/4,2/10,2/8,2/4,2/8)
#' f=c(8,8,12,18,14,18,10,10,6,16,16,12,12,8,8,6,18,14,6,14)
#' metaanalyseSmallSampleSizeExperiments(d,f,A)
# A tibble: 1 x 3 UnweightedMean WeightedMean VarWeightedMean <dbl> <dbl> <dbl> 1 0.294 0.320 0.0156
metaanalyseSmallSampleSizeExperiments <- function(d, f, A) {
  c <- reproducer::calculateSmallSampleSizeAdjustment(f)
  g <- d * c
  meang.unweighted <- mean(g)

  a <- calc.a(f, A)
  b <- calc.b(f)

  # Hedges recommended analysis
  temp.gvar <- meang.unweighted ^ 2 * b + a
  gweight <- 1 / temp.gvar / sum(1 / temp.gvar)

  meang.weighted <- sum(g * gweight)
  g.exact.var <- 1 / sum(1 / temp.gvar)


  output <-
    tibble::as_tibble(
      dplyr::bind_cols(
        UnweightedMean = meang.unweighted,
        WeightedMean = meang.weighted,
        VarWeightedMean = g.exact.var
      )
    )

  return(output)
}


#' @title AnalyseResiduals
#' @description The function calculates sample statistics based on the residuals from a specified experiment
#' @author Barbara Kitchenham and Lech Madeyski
#' @importFrom grDevices boxplot.stats
#' @import stats
#' @export AnalyseResiduals
#' @param Residuals a vector of residuals
#' @param ExperimentName a character string identifying the data set
#' @return  A dataframe identifying the ExperimentName and its associated sample parameter: Length, Mean, Median, Variance, Standard deviation, skewness, kurtosis, the outcome of the Shapiro and Anderson-Darling normality test and the number of outliers.
#' @examples
#' ExpData=rnorm(30,0,1)
#' set.seed(123)
#' AnalyseResiduals(Residuals=ExpData,ExperimentName='ExpName')
#' #  ExperimentName       Mean      Median  Variance   Skewness Kurtosis ShapiroTest AndersonDarling
#' #1        ExpName -0.1396192 -0.01943395 0.8424521 -0.1964175 4.559587   0.1608315       0.1316835
#' #  NumOut
#' #  1
AnalyseResiduals <-
  function(Residuals, ExperimentName = "ExpName") {
    n <- length(Residuals)
    Res.Mean <- mean(Residuals)
    Res.Median <- median(Residuals)
    Res.Var <- var(Residuals)
    Res.Std <- sqrt(Res.Var)
    Res.Skew <- sum((Res.Mean - Residuals) ^ 3) / n
    Res.Skew <- Res.Skew / Res.Std ^ 3
    Res.Kurt <- sum((Res.Mean - Residuals) ^ 4) / n
    Res.Kurt <- Res.Kurt / Res.Std ^ 4
    # Res.ShapTes=rstatix::shapiro_test(Residuals)$p.value
    Res.ShapTes <- stats::shapiro.test(Residuals)$p.value
    Res.AndersonDarling <- nortest::ad.test(Residuals)$p.value
    # Calculate number of outliers

    NumOut <- length(boxplot.stats(Residuals)$out)
    ResStats <- data.frame(
      Mean = Res.Mean,
      Median = Res.Median,
      Variance = Res.Var,
      Skewness = Res.Skew,
      Kurtosis = Res.Kurt,
      ShapiroTest = Res.ShapTes,
      AndersonDarling = Res.AndersonDarling,
      NumOut = NumOut
    )

    ResStats <- cbind(ExperimentName, ResStats)

    return(ResStats)
  }


#' @title LaplaceDist
#' @description Returns a sample of N observations from a Laplace distribution with specified mean and spread.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export LaplaceDist
#' @param N is the required sample size
#' @param mean is the required mean
#' @param spread is the spread of the function
#' @param max the upper limit of the distribution. Must be finite.
#' @param min the lower limit of the distribution. Must be finite.
#' @return N values from a Laplace distribution
#' @examples
#' set.seed(123)
#' LaplaceDist(10, 0, 1)
#' #  [1] -0.55311564  0.85946218 -0.20094937  1.45258293  2.12808209 -2.39565480  0.05785263
#' #   [8]  1.53636446  0.10855453 -0.09076809
#'
LaplaceDist <- function(N,
                        mean,
                        spread,
                        max = 0.5,
                        min = -0.5) {
  y <-
    stats::runif(N, min, max) # Get data from a uniform distribution
  x <- mean - spread * sign(y) * log(1 - 2 * abs(y))

  return(x)
}


#' @title simulate2GExperimentData
#' @description The function returns a two group data set based on one of four different distributions.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export simulate2GExperimentData
#' @param mean The mean (or rate for gamma data) of the baseline distribution
#' @param sd The standard deviation (or shape for gamma data) of the baseline distribution
#' @param diff The adjustment to the baseline mean for the alternative distribution.
#' @param GroupSize An integer defining the number of data items in each group.
#' @param type A string identifying the distribution used to simulate the data: 'n' for normal, 'ln' for log-normal, 'g' for gamma, 'lap' for Laplace.
#' @param ExpAdj An additional adjustment factor that is added to both the mean value. Defaults to zero.
#' @param StdAdj An additional adjustment factor that is added to both group variance (or rate for gamma data). Defaults to zero.
#' @param BlockEffect An additional factor that is added to the mean of the both groups (shape for the gamma distribution). Defaults to zero.
#' @param BlockStdAdj An additional factor that is added to the variance of both groups (shape for the gamma distribution). Defaults to zero.
#' @return A table with two columns (BaselineData and AlternativeData) holding the data for each group. For lognormal data an additional two columns are added which return the log transformed data.
#' @examples
#' set.seed(236)
#' simulate2GExperimentData(mean = 0, sd = 1, diff = 0.5, GroupSize = 10,
#'   type = "n", ExpAdj = 0, StdAdj = 0, BlockEffect = 0, BlockStdAdj = 0)
#' # A tibble: 10 x 2
#' #    BaselineData AlternativeData
#' #           <dbl>           <dbl>
#' #          <dbl>           <dbl>
#' # 1      -0.285           -0.255
#' # 2      -0.972            0.112
#' # 3      -0.549            1.36
#' # 4       1.05             1.47
#' # 5      -0.267            0.107
#' # 6      -0.137            0.395
#' # 7       1.30             1.27
#' # 8      -0.722            1.70
#' # 9      -0.525            0.264
#' # 10      -0.0222           0.787
#' set.seed(345)
#' simulate2GExperimentData(mean = 0, sd = 1, diff = 0.5, GroupSize = 10,
#'   type = "l", ExpAdj = 0, StdAdj = 0, BlockEffect = 0, BlockStdAdj = 0)
#' # A tibble: 10 x 4
#' #    BaselineData AlternativeData transBaselineData transAlternativeData
#' #          <dbl>           <dbl>             <dbl>                <dbl>
#' # 1        0.456          10.7             -0.785                 2.37
#' # 2        0.756           0.407           -0.280                -0.900
#' # 3        0.851           0.705           -0.161                -0.350
#' # 4        0.748           2.27            -0.291                 0.818
#' # 5        0.935           4.07            -0.0675                1.40
#' # 6        0.531           0.405           -0.634                -0.903
#' # 7        0.395           2.91            -0.928                 1.07
#' # 8        5.53            4.69             1.71                  1.55
#' # 9        5.23            0.602            1.65                 -0.508
#' # 10        6.11            2.23             1.81                  0.802
#'
simulate2GExperimentData <-
  function(mean,
           sd,
           diff,
           GroupSize,
           type = "n",
           ExpAdj = 0,
           StdAdj = 0,
           BlockEffect = 0,
           BlockStdAdj = 0) {
    if (type == "n") {
      BaselineData <-
        stats::rnorm(GroupSize, mean + ExpAdj + BlockEffect, sd + BlockStdAdj)
      AlternativeData <-
        stats::rnorm(GroupSize,
                     mean + ExpAdj + BlockEffect + diff,
                     sd + StdAdj + BlockStdAdj)
    }
    if (type == "g") {
      rate1 <- mean + ExpAdj
      rate2 <- rate1 + diff + StdAdj
      shape1 <- sd + BlockStdAdj + BlockEffect
      BaselineData <- stats::rgamma(GroupSize, shape1, rate1)
      AlternativeData <- stats::rgamma(GroupSize, shape1, rate2)
    }
    if (type == "l") {
      BaselineData <-
        stats::rlnorm(GroupSize, mean + ExpAdj + BlockEffect, sd + BlockStdAdj)
      AlternativeData <-
        stats::rlnorm(GroupSize,
                      mean + ExpAdj + BlockEffect + diff,
                      sd + StdAdj + BlockStdAdj)
      transBaselineData <- log(BaselineData)
      transAlternativeData <- log(AlternativeData)
    }

    if (type == "lap") {
      # Changed to use built-in defaults for max and min
      BaselineData <-
        LaplaceDist(GroupSize, mean + ExpAdj + BlockEffect, sd + BlockStdAdj)
      AlternativeData <-
        LaplaceDist(GroupSize,
                    mean + ExpAdj + diff + BlockEffect,
                    sd + StdAdj + BlockStdAdj)
    }

    OutputData <-
      dplyr::bind_cols(BaselineData = BaselineData, AlternativeData = AlternativeData)
    if (type == "l") {
      OutputData <-
        dplyr::bind_cols(OutputData,
                         transBaselineData = transBaselineData,
                         transAlternativeData = transAlternativeData)
    }

    return(OutputData)
  }


#' @title calculateLargeSampleRandomizedDesignEffectSizes
#' @description The function simulates a large experiment  to estimate the asymptotic values of the probability of superiority, Cliff's d and the standardized mean difference data for a two group randomized experiment for four different distributions: Normal (i.e. type="n"), log-normal (i.e. type="l"), gama (i.e. tyep="g") and Laplace (i.e., type="lap").
#' @author Barbara Kitchenham and Lech Madeyski
#' @export calculateLargeSampleRandomizedDesignEffectSizes
#' @param meanC to act as the mean of the distribution used to generate the control group data (default 0) (note for the gamma distribution this is the rate parameter and must not be zero)
#' @param sdC the variance/spread of the distribution used to generate the control group data (default 1).
#' @param diff a value added to meanC to generate the treatment group data (default 0).
#' @param N the size of each group (default 5000000)
#' @param type the distribution of the data to be generated (default "n").
#' @param StdAdj a value that can be added to sdC to introduce heterogeneity into the treatment group (default 0).
#' @param reporttrans If set to "Yes" AND type="l" the algorithm returns the values obtained by analysing applying the logarithmic transformation to the simulated data (default "No").
#' @return A tibble identifying the sample statistics and the values of the probability of superiority, Cliff's d and StdMD (labelled StdES)
#' @examples
#'set.seed=400
#'calculateLargeSampleRandomizedDesignEffectSizes(meanC=0, sdC=1, diff=.5,
#' N=10000, type="n",StdAdj = 0) #N=100000, type="n",StdAdj = 0)
#'# A tibble: 1 x 9
#'#     MeanC   SdC MeanT   SdT  Phat Cliffd   UES   Var StdES
#'#     <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl>
#'#1  0.00642  1.00 0.519 0.995 0.642  0.284 0.513 0.996 0.514
#'#1     0     1   0.5     1 0.637  0.275 0.499  1.01 0.497
#'as.data.frame(calculateLargeSampleRandomizedDesignEffectSizes(meanC=0, sdC=1, diff=0.707104,
#' N=100000, type="l",StdAdj = 0,reporttrans="Yes"))
#' #N=1000000, type="l",StdAdj = 0,reporttrans="Yes"))
#'#     MeanC     StdC    MeanT     StdT      Phat    Cliffd      UES      Var
#'#1 1.647446 2.219114  3.33124 4.404537 0.6926779 0.3853558 1.683795 12.1622
#'#       StdES   MeanCTrans MeanTTrans StdCTrans StdTTrans PhatTrans CliffdTrans
#'#1  0.4828175 -0.004298487  0.7066049  1.001199 0.9963736 0.6926779   0.3853558
#'#   UESTrans VarTrans StdESTrans
#'#1 0.7109034  0.99758  0.7117651

calculateLargeSampleRandomizedDesignEffectSizes <-
  function (meanC = 0,
            sdC = 1,
            diff = 0,
            N = 5e+06,
            type = "n",
            StdAdj = 0,
            reporttrans = "No")
  {
    DataSet <-
      reproducer::simulate2GExperimentData(
        mean = meanC,
        sd = sdC,
        diff = diff,
        GroupSize = N,
        type = type,
        ExpAdj = 0,
        StdAdj = StdAdj,
        BlockEffect = 0,
        BlockStdAdj = 0
      )
    UES <-
      base::mean(DataSet$AlternativeData) - base::mean(DataSet$BaselineData)
    Var <-
      (stats::var(DataSet$AlternativeData) + stats::var(DataSet$BaselineData)) /
      2
    StdES <- UES / sqrt(Var)
    B1res <-
      calculatePhat(DataSet$BaselineData, DataSet$AlternativeData)
    Phat <- B1res$phat
    Cliffd <- B1res$dhat
    MeanC <- base::mean(DataSet$BaselineData)
    MeanT <- base::mean(DataSet$AlternativeData)
    StdT = sqrt(stats::var(DataSet$AlternativeData))
    StdC = sqrt(stats::var(DataSet$BaselineData))

    res <-
      as.data.frame(
        cbind(
          MeanC = MeanC,
          StdC = StdC,
          MeanT = MeanT,
          StdT = StdT,
          Phat = Phat,
          Cliffd = Cliffd,
          UES = UES,
          Var = Var,
          StdES = StdES
        )
      )

    if (type == "l" & reporttrans == "Yes") {
      UESTrans = base::mean(DataSet$transAlternativeData) - base::mean(DataSet$transBaselineData)
      VarTrans <-
        (
          stats::var(DataSet$transAlternativeData) + stats::var(DataSet$transBaselineData)
        ) / 2
      StdESTrans <- UESTrans / sqrt(VarTrans)
      B1resTrans <-
        calculatePhat(DataSet$transBaselineData, DataSet$transAlternativeData)
      PhatTrans <- B1resTrans$phat
      CliffdTrans <- B1resTrans$dhat
      MeancTrans <- base::mean(DataSet$transBaselineData)
      MeantTrans <- base::mean(DataSet$transAlternativeData)
      StdTTrans <- sqrt((stats::var(DataSet$transAlternativeData)))
      StdCTrans <- sqrt((stats::var(DataSet$transBaselineData)))
      res <-
        as.data.frame(
          cbind(
            res,
            MeanCTrans = MeancTrans,
            MeanTTrans = MeantTrans,
            StdCTrans = StdCTrans,
            StdTTrans = StdTTrans,
            PhatTrans = PhatTrans,
            CliffdTrans = CliffdTrans,
            UESTrans = UESTrans,
            VarTrans = VarTrans,
            StdESTrans = StdESTrans
          )
        )

    }

    res <- tibble::as_tibble(res)
    return(res)
  }



#'  @title calculateKendalltaupb
#'  @description  Computes point bi-serial version of  Kendall's tau plus a 1-alpha confidence interval using the method recommended by Long and Cliff (1997).  The algorithm is based on Wilcox's code but was extended to return the consistent variance and the confidence intervals based on the t-distribution. Also added a Diagnostic parameter to output internal calculations.
#' @author Rand Wilcox, Barbara Kitchenham and Lech Madeyski
#' @export calculateKendalltaupb
#' @param x either a matrix with two columns containing two correlated variables or a vector of variables
#' @param y if y=NULL, assume x is a matrix with two columns, otherwise y is a vector of variables with x[i] and y[i] being from the same experimental unit
#' @param alpha the Type 1 error level used for statistical tests (default 0.05)
#' @param alternative The type of statistical test. Valid values are one of
#' c('two.sided', 'greater', 'less')
#' @return list containing the estimate of Kendall's tau, the consistent variance of tau and its confidence intervals based on the t-test (recommended by Long and Cliff)
#' @examples
#' set.seed(123)
#' a <- c(1.2, 3, 1.8, 2, 2, 0.5, 0.5, 1, 3, 1)
#' b <- c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0)

#' calculateKendalltaupb(a,b,alpha=.05)
#' #$cor
#' #[1] 0.3555556

#' #$cit
#' #[1] -0.1240567  0.5555556

#' #$n
#' #[1] 10

#' #$df
#' #[1] 7

#' #$consistentvar
#' #[1] 0.04113925

#' #$sig
#' #[1] FALSE

#' set.seed(234)
#' a2=c(rnorm(10,0,1),rnorm(10,0.5,1))
#' b2=c(rep(0,10),rep(1,10))
#' calculateKendalltaupb(a2,b2,alpha=.05,alternative='greater')
#' #$cor
#' #[1] 0.2842105

#' #$cit
#' #[1] 0.06517342 0.52631579
#' #$n
#' #[1] 20

#' #$df
#' #[1] 17

#' #consistentvar
#' #1] 0.01585379

#' #$sig
#' #[1] TRUE

#' calculateKendalltaupb(a2,b2,alpha=.05,alternative='less')
#' #$cor
#' #[1] 0.2842105

#' #$cit
#' #[1] -0.5263158  0.5032476

#' #$n
#' #[1] 20

#' #$df
#' #[1] 17

#' #$consistentvar
#' #[1] 0.01585379

#' #$sig
#' #[1] FALSE

calculateKendalltaupb <-
  function(x,
           y = NULL,
           alpha = 0.05,
           alternative = "two.sided") {
    if (!checkIfValidDummyVariable(x) & !checkIfValidDummyVariable(y)) {
      stop("No valid dummy variable - use only 1 and 0")
    }
    if (length(x) <= 3) {
      stop("Too few data points")
    }
    if (length(x) != length(y)) {
      stop("Invalid input vectors")
    }
    if (length(x) != length(x[!is.na(x)])) {
      stop("Missing values not permitted")
    }
    if (length(y) != length(y[!is.na(y)])) {
      stop("Missing values not permitted")
    }

    # Needs a test to ensure one of the variables contains only 1 or 0 values

    m <- cbind(x, y)

    x <- m[, 1]
    y <- m[, 2]
    xdif <- base::outer(x, x, FUN = "-")
    ydif <- base::outer(y, y, FUN = "-")
    tv <- sign(xdif) * sign(ydif)

    # Corrects error in Wilcox's algorithm
    n <- length(x)
    dbar <- base::apply(tv, 1, sum) / (n - 1)

    tau <- sum(tv) / (n * (n - 1))

    A <- sum((dbar - tau) ^ 2) / (n - 1)
    B <- (n * (n - 1) * (-1) * tau ^ 2 + sum(tv ^ 2)) / (n ^ 2 - n - 1)
    C <-
      (4 * (n - 2) * A + 2 * B) / (n * (n - 1)) # C is the consistent variance

    # Confidence interval based on t distribution - recommended by Long and Cliff
    vv <- stats::qt(alpha / 2, 2 * n - 3)

    cilowt <- tau + vv * sqrt(C)
    if (cilowt < (-n / (2 * (n - 1)))) {
      cilowt <- -n / (2 * (n - 1))
    }
    cihit <- tau - vv * sqrt(C)
    if (cihit > n / (2 * (n - 1))) {
      cihit <- n / (2 * (n - 1))
    }

    testfunctionParameterChecks(alternative = alternative,
                                alpha = alpha,
                                stderr = sqrt(C))

    tau.Stats <- calcEffectSizeConfidenceIntervals(
      effectsize = tau,
      effectsize.variance = C,
      effectsize.df = n - 3,
      alpha = alpha,
      alternative = alternative,
      UpperValue = n / (2 * (n - 1)),
      LowerValue = -n / (2 * (n - 1))
    )

    list(
      cor = tau,
      cit = c(tau.Stats$ES.ci.lower, tau.Stats$ES.ci.upper),
      n = n,
      df = n - 3,
      consistentvar = C,
      sig = tau.Stats$ES.sig
    )
  }





#' @title simulateRandomizedDesignEffectSizes
#' @description This simulates one of four data distributions (normal, log-normal, gamma and Laplace), and finds the values of phat and Cliffs d and their variances. It assumes equal group sizes. It returns values of the effect sizes and their variance for a simulated randomized experiment with two treatments.  It returns whether or not each non-parametric effect size was significant. It also returns the parametric (standardized and unstandardized) Effect Size and the whether the t-test was significant.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export simulateRandomizedDesignEffectSizes
#' @param mean The mean used for one of the treatment groups (this is the rate for the gamma data)
#' @param sd The spread used for both treatment groups. It mus be a real value greater than 0 (this is the shape for the gamma data).
#' @param diff This is added to the parameter mean, to define the mean of the other treatment group. It can be a real value avd can take the value zero.
#' @param N this is the number of observations in each group. It must be an integer greater than 3.
#' @param type this specifies the underlying distribution used to generate the data. it takes the values 'n' for a normal distribution, 'l' for lognormal distribution,'g' for a gamma distribution, 'lap' for a Laplace distribution.
#' @param StdAdj this specifies the extent of variance instability to be introduced.
#' @param alpha the level for all statistical tests (default 0.05)
#' @param AlwaysTwoSidedTests if set to FALSE (i.e. default) the algorithms uses one-sided tests if diff!=0 and two-sided tests otherwise. If set to TRUE the algorithm always uses two-sided tests.
#' @param Return.Data if set to true the algorithm returns the data not the effect sizes (default FALSE).
#' @return data frame incl. the non-parametric and parametric effect sizes and whether the effect sizes are significant at the specified alpha level. For log-normal data the function returns the effect sizes for the transformed data.
#' @examples
#' set.seed(123)
#' as.data.frame(
#'   simulateRandomizedDesignEffectSizes(
#'     mean = 0, sd = 1, diff = 0.8, N = 10, type = "n", StdAdj = 0))
#' #   phat    varphat   dfphat sigphat   d       vard sigd       cor     varcor sigCVt  t.value
#' # 1 0.75 0.01522222 17.46405    TRUE 0.5 0.06237576 TRUE 0.2631579 0.01754995   TRUE 2.095142
#' #      t.se     t.df      t.lb t.ub t.sig        ES  Variance     StdES  MedDiff
#' # 1 0.4457915 17.87244 0.1606665  Inf  TRUE 0.9339963 0.9936502 0.9369759 1.260127
#' set.seed(123)
#' as.data.frame(
#'   simulateRandomizedDesignEffectSizes(
#'     mean = 0, sd = 1, diff = 0.8, N = 10, type = "n", StdAdj = 0,
#'     AlwaysTwoSidedTests = TRUE))
#' #  phat    varphat   dfphat sigphat   d       vard  sigd       cor
#' # 1 0.75 0.01522222 17.46405   FALSE 0.5 0.06237576 FALSE 0.2631579
#' #      varcor sigCVt  t.value      t.se     t.df         t.lb     t.ub t.sig
#' # 1 0.01754995  FALSE 2.095142 0.4457915 17.87244 -0.003056196 1.871049 FALSE
#' #         ES  Variance     StdES  MedDiff
#' # 1 0.9339963 0.9936502 0.9369759 1.260127
#' set.seed(456)
#' as.data.frame(
#'   simulateRandomizedDesignEffectSizes(
#'     mean = 0, sd = 1, diff = 0.8, N = 10, type = "l", StdAdj = 0))
#' # phat     varphat  dfphat sigphat    d      vard sigd       cor     varcor
#' # 1 0.87 0.008466667 11.1111    TRUE 0.74 0.0350497 TRUE 0.3894737 0.01039674
#' #  sigCVt  t.value     t.se     t.df     t.lb t.ub t.sig       ES Variance
#' # 1   TRUE 3.599375 2.148297 9.312472 3.809448  Inf  TRUE 7.732529 23.07591
#' #    StdES MedDiff transttest  EStrans StdEStrans VarTrans
#' # 1 1.60969 7.77893   0.998772 1.731323   1.598065 1.173728
#'
#' set.seed(123)
#' as.data.frame(
#'   simulateRandomizedDesignEffectSizes(
#'     mean = 0, sd = 1, diff = 0.8, N = 10, type = "n", StdAdj = 0,
#'     Return.Data = TRUE))
#' #   BaselineData AlternativeData
#' # 1   -0.69470698       1.0533185
#' # 2   -0.20791728       0.7714532
#' # 3   -1.26539635       0.7571295
#' # 4    2.16895597       2.1686023
#' # 5    1.20796200       0.5742290
#' # 6   -1.12310858       2.3164706
#' # 7   -0.40288484      -0.7487528
#' # 8   -0.46665535       1.3846137
#' # 9    0.77996512       0.9238542
#' # 10  -0.08336907       1.0159416
simulateRandomizedDesignEffectSizes <-
  function(mean,
           sd,
           diff,
           N,
           type = "n",
           StdAdj = 0,
           alpha = 0.05,
           AlwaysTwoSidedTests = FALSE,
           Return.Data = FALSE) {
    # 12-03-2022 Corrected spelling of function name

    conf.level <- 1 - alpha
    # Allowed for one-sided tests
    if (AlwaysTwoSidedTests | (diff == 0)) {
      alternative <- "two.sided"
    } else {
      # Find direction of MD

      if (type == "g") {
        # A negative value for diff decreases the rate & the MD is positive
        Positive.MD <- (diff < 0)
      } else {
        Positive.MD <- (diff > 0)
      }
      if (Positive.MD) {
        alternative <- "greater"
      } else {
        alternative <- "less"
      }
    }

    DataSet <- simulate2GExperimentData(
      mean = mean,
      sd = sd,
      diff = diff,
      GroupSize = N,
      type = type,
      ExpAdj = 0,
      StdAdj = StdAdj,
      BlockEffect = 0,
      BlockStdAdj = 0
    )


    if (Return.Data) {
      output <- tibble::as_tibble(DataSet)
    } else {
      # Remove_ktau Calculate the Kendall's tau value for the data dummy = c(rep(0, N), rep(1, N)) xy = c(DataSet$BaselineData,
      # DataSet$AlternativeData) Calculate ktau for the data set ktest = calculateKendalltaupb( dummy, xy, alpha = alpha, alternative =
      # alternative )


      # Calculate Cliff's d for the data set and determine whether the value is significant
      d.test <-
        Cliffd.test(
          DataSet$AlternativeData,
          DataSet$BaselineData,
          alternative = alternative,
          alpha = alpha
        )
      # Calculate phat statistics for the data set


      ptest <-
        PHat.test(
          DataSet$BaselineData,
          DataSet$AlternativeData,
          alpha = alpha,
          alternative = alternative
        )

      # Add the results of a t-test as a baseline


      res <-
        t.test(
          DataSet$AlternativeData,
          DataSet$BaselineData,
          alternative = alternative,
          conf.level = conf.level
        )

      UES <-
        base::mean(DataSet$AlternativeData) - base::mean(DataSet$BaselineData)
      Var <-
        (stats::var(DataSet$BaselineData) + stats::var(DataSet$AlternativeData)) / 2
      MedDiff <-
        stats::median(DataSet$AlternativeData) - stats::median(DataSet$BaselineData)
      EffectSize <- UES / sqrt(Var)
      pval <- res$p.value
      t.sig <- (pval < 0.05)

      se <- as.numeric(res$stderr)
      df <- as.numeric(res$parameter)
      t.val <- UES / se

      StandardMetrics <- tibble::tibble(
        phat = ptest$phat,
        varphat = ptest$sqse.phat,
        dfphat = ptest$phat.df,
        sigphat = ptest$phat.sig,
        d = d.test$d,
        vard = d.test$sqse.d,
        sigd = d.test$d.sig,
        t.value = t.val,
        t.se = as.numeric(res$stderr),
        t.df = as.numeric(res$parameter),
        t.lb = as.numeric(res$conf.int[1]),
        t.ub = as.numeric(res$conf.int[2]),
        t.sig = t.sig,
        ES = UES,
        Variance = Var,
        StdES = EffectSize,
        MedDiff = MedDiff
      )

      if (type == "l") {
        # Check that log-normal data gives appropriate values after transformation
        trans.t <-
          t.test(
            DataSet$transBaselineData,
            DataSet$transAlternativeData,
            alt = alternative,
            conf.level = conf.level
          )
        pval.trans <- trans.t$p.value
        ES.trans <-
          base::mean(DataSet$transAlternativeData) - base::mean(DataSet$transBaselineData)
        VarTrans <-
          (
            stats::var(DataSet$transBaselineData) + stats::var(DataSet$transAlternativeData)
          ) / 2
        StdES.trans <- ES.trans / sqrt(VarTrans)

        AdditionalMetrics <-
          cbind(
            transttest = pval.trans,
            EStrans = ES.trans,
            StdEStrans = StdES.trans,
            VarTrans = VarTrans
          )

        AdditionalMetrics <- tibble::as_tibble(AdditionalMetrics)
        StandardMetrics <-
          dplyr::bind_cols(StandardMetrics, AdditionalMetrics)
      }
      output = StandardMetrics
    }
    return(output)
  }


#' @title RandomExperimentSimulations
#' @description This function performs multiple simulations of two-group balanced experiments for one of four distributions and a specific group size. It identifies the average value of phat, Cliff' d and their variances. It either returns the effect sizes for each non-parametric effect size or it reports the number of times the each non-parametric effect size is assessed to be significantly different from zero. We also present the values for the t-test as a comparison. For log-normal data the results of analysing the transformed data are also reported.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export RandomExperimentSimulations
#' @param mean The default mean used for both groups (one treatment group and one control group). It can be changed for the treatment group using  the parameter diff
#' @param sd This is the default spread for both groups. It must be a real value greater than 0. It can be adjusted for the treatment group using the parameter StdAdj
#' @param diff This is added to the treatment group mean. It can be a real value avd can take the value zero.
#' @param N this is the number of observations in each group. It must be an integer greater than 3.
#' @param reps this identifies the number of times each experiment simulation is replicated.
#' @param type this specifies the underlying distribution used to generate the data. It takes the values 'n' for a normal distribution, 'l' for lognormal distribution,'g' for a gamma distribution, 'lap' for a Laplace distribution.
#' @param seed This specifies the initial seed for the set of replications
#' (default 123).
#' @param StdAdj this specifies the extent of variance instability introduced by the treatment and it must be non-negative but can be 0.
#' @param alpha This specifies the level of significance used for statistical
#' tests (default 0.05).
#' @param returnData If TRUE, the function returns the individual effect sizes and their variances, otherwise it returns summary statistics (default FALSE).
#' @param AlwaysTwoSidedTests If set to FALSE (default) the algorithms uses one-sided tests if diff!=0 and two-sided tests if diff=0. If set to TRUE the algorithm always uses two-sided tests.
#' @examples
#' as.data.frame(
#'   RandomExperimentSimulations(
#'     mean = 0, sd = 1, diff = 0.5, N = 20, reps = 50, type = "n",
#'     seed = 123, StdAdj = 0, alpha = 0.05))
#' #        phat     phatvar sigphat emp.phat.var       d       dvar sigd
#' # 1  0.636675 0.007980072    0.38  0.006413391 0.27335 0.03257962 0.36
#' #    emp.d.var   tpower        ES Variance     StdES   MedDiff
#' #1  0.02565356     0.41 0.4849609 0.988889 0.4982554 0.4666802
#' #as.data.frame(
#'  # RandomExperimentSimulations(
#'  #   mean = 0, sd = 1, diff = 0.5, N = 20, reps = 500, type = "n",
#'  #   seed = 123, StdAdj = 0, alpha = 0.05))
#' #     phat     phatvar sigphat emp.phat.var      d       dvar  sigd  emp.d.var
#' # 1 0.63915 0.007925803   0.444  0.007904962 0.2783 0.03235111 0.414 0.03161985
#' #     tpower        ES Variance
#' # 1     0.444 0.4999034 1.002012
#' # 1      StdES   MedDiff
#' # 1 0.5099792 0.4901394
#'
#' #as.data.frame(
#' #   RandomExperimentSimulations(
#' #     mean = 0, sd = 1, diff = 0.2, N = 20, reps = 500, type = "n",
#' #     seed = 123, StdAdj = 0, alpha = 0.05, AlwaysTwoSidedTests = TRUE))
#' #     phat     phatvar sigphat emp.phat.var       d       dvar  sigd emp.d.var
#' # 1 0.55762 0.008596555   0.092  0.008457325 0.11524 0.03505528 0.076 0.0338293
#' #     tpower        ES Variance     StdES   MedDiff
#' # 1       0.1 0.1999034 1.002012 0.2043908 0.1901394
#'
#' #as.data.frame(
#' #   RandomExperimentSimulations(
#' #     mean = 0, sd = 1, diff = 0.2, N = 20, reps = 500, type = "n",
#' #     seed = 123, StdAdj = 0, alpha = 0.05, AlwaysTwoSidedTests = FALSE))
#' #     phat     phatvar sigphat emp.phat.var       d       dvar  sigd emp.d.var
#' # 1 0.55762 0.008596555   0.154  0.008457325 0.11524 0.03505528 0.146 0.0338293
#' #        tpower        ES Variance
#' # 1         0.16 0.1999034 1.002012
#' #      StdES   MedDiff
#' # 1 0.2043908 0.1901394
#'
#' RandomExperimentSimulations(
#'   mean = 0, sd = 1, diff = 0.5, N = 20, reps = 10, type = "l", seed = 456,
#'   StdAdj = 0, alpha = 0.05, returnData = TRUE, AlwaysTwoSidedTests = FALSE)
#' # A tibble: 10 x 6
#' #   Cliffd CliffdSig  PHat PHatSig  StdES ESSig
#' #    <dbl>     <dbl> <dbl>   <dbl>  <dbl> <dbl>
#' # 1 -0.185         0 0.407       0 -0.246     0
#' # 2 -0.08          0 0.46        0  0.185     0
#' # 3  0.1           0 0.55        0  0.149     0
#' # 4  0.42          1 0.71        1  0.885     1
#' # 5  0.51          1 0.755       1  0.827     1
#' # 6  0.185         0 0.592       0  0.628     1
#' # 7  0.465         1 0.732       1  0.818     1
#' # 8  0.42          1 0.71        1  0.341     0
#' # 9  0.37          1 0.685       1  0.419     0
#' # 10  0.115         0 0.557       0  0.273     0
#'
RandomExperimentSimulations <-
  function(mean,
           sd,
           diff,
           N,
           reps,
           type = "n",
           seed = 123,
           StdAdj = 0,
           alpha = 0.05,
           returnData = FALSE,
           AlwaysTwoSidedTests = FALSE) {
    phatsum <-
      0 # This is used to sum the value of phat across the replications
    phatvarsum <-
      0 # This is used to sum the value of the variance of phat across the replications
    sig.phat <-
      0 # This is used to sum the number of times pat is significant across the replications
    phatsquare <-
      0 # This  is used to sum phat^2 and construct an empirical variance of phat
    dsum <-
      0 # This is used to sum the value of Cliff's d across the replications
    dvarsum <-
      0 # This is used to sum the value of the variance of Cliff's d across the replications
    sig.d <-
      0 # This is used to sum the number of times d is significant across the replications
    dsquare <-
      0 # This  is used to sum d^2 and construct an empirical variance of d.

    # Remove_ktau ksum = 0 # This is used to sum the value of the point biserial tau across the replications kvarsum = 0 # This is used to sum the
    # value of the variance of the point biserial tau across the replications ksquare = 0 # This is used to sum the square of tau_pb across the
    # replications and construct an empirical variance ksigCVt = 0

    tsig <-
      0 # This is used to count the number of significant t values across the replications
    ES <-
      0 # This is used to sum the value of the parametric effect size (unstandardized) across the replications
    StdES <-
      0 # This is used to sum the value of the parametric effect size (standardized) across the replications
    ESSig <- 0 # Used to assess whether the t-test was significant
    ESCISig <- 0 # Check the CI results
    Var <-
      0 # This is used to sum the value of the variance across replications
    ES.l.trans <-
      0 # This is used to sum the transformed unstandardized effect size for lognormal data sets
    StdES.l.trans <-
      0 # This is used to sum the transformed standardized effect size for lognormal data sets
    Var.l.trans <-
      0 # This is used to sum the transformed variance for lognormal data sets


    MedDiff <- 0 # Used to hold the median difference

    DataTable <- NULL

    base::set.seed(seed)
    for (i in 1:reps) {
      # Call the program that generates the random data sets and calculates the sample statistics.  12-03-2022 Corrected name of function res =
      # simulateRandomzedDesignEffectSizes(mean, sd, diff, N, type, StdAdj)
      res <-
        simulateRandomizedDesignEffectSizes(
          mean = mean,
          sd = sd,
          diff = diff,
          N = N,
          type = type,
          StdAdj = StdAdj,
          alpha = alpha,
          AlwaysTwoSidedTests = AlwaysTwoSidedTests
        )
      if (returnData == FALSE) {
        # Aggregate data to provide counts of significance and overall effect size averages Cliff's d
        dsum <- dsum + res$d
        dvarsum <- dvarsum + res$vard
        if (res$sigd) {
          sig.d <- sig.d + 1
        }
        dsquare <- dsquare + res$d ^ 2

        # Probability of superiority
        phatsum <- phatsum + res$phat
        phatvarsum <- phatvarsum + res$varphat
        if (res$sigphat) {
          sig.phat <- sig.phat + 1
        }
        phatsquare <- phatsquare + res$phat ^ 2

        # Parametric statistics This needs a change to simulateRandomizedDesignEffectSizes to return the significance not the pvalue

        ES <- ES + res$ES
        StdES <- StdES + res$StdES
        MedDiff <- MedDiff + res$MedDiff
        ESSig <- if (res$t.sig) {
          1
        } else {
          0
        }

        tsig <- tsig + ESSig

        Var <- Var + res$Variance

        if (type == "l") {
          ES.l.trans <- ES.l.trans + res$EStrans
          StdES.l.trans <- StdES.l.trans + res$StdEStrans
          Var.l.trans <- Var.l.trans + res$VarTrans
        }
      } else {
        # Store the outcome from each replication
        ESSig <- if (res$t.sig) {
          1
        } else {
          0
        }

        DataTable <-
          tibble::tibble(dplyr::bind_rows(
            DataTable,
            dplyr::bind_cols(
              Cliffd = res$d,
              CliffdSig = as.numeric(res$sigd),
              PHat = res$phat,
              PHatSig = as.numeric(res$sigphat),
              StdES = res$StdES,
              ESSig = ESSig
            )
          ))
      }
    }

    if (returnData == FALSE) {
      # Calculate averages of the statistics across the replications
      d <- dsum / reps
      dvar <- dvarsum / (reps)
      sigd <- sig.d / (reps)
      emp.d.var <- (dsquare - reps * d ^ 2) / (reps - 1)

      phat <- phatsum / reps
      phatvar <- phatvarsum / (reps)
      sigphat <- sig.phat / (reps)
      emp.phat.var <- (phatsquare - reps * phat ^ 2) / (reps - 1)

      ES <- ES / reps
      StdES <- StdES / reps
      MedDiff <- MedDiff / reps
      Variance <- Var / reps
      tpower <- tsig / reps

      outcome <-
        tibble::tibble(
          phat,
          phatvar,
          sigphat,
          emp.phat.var,
          d,
          dvar,
          sigd,
          emp.d.var,
          tpower,
          ES,
          Variance,
          StdES,
          MedDiff
        )


      if (type == "l") {
        # This is used for validation that the algorithms are consistent. The statistics from the transformed lognormal data can be compared
        # with the statistics from the normal data.
        ESLog <- ES.l.trans / reps
        StdESLog <- StdES.l.trans / reps
        VarLog <- Var.l.trans / reps
        AdditionalMetrics <-
          cbind(ESLog = ESLog,
                StdESLog = StdESLog,
                VarLog = VarLog)

        AdditionalMetrics <- tibble::as_tibble(AdditionalMetrics)
        outcome <- dplyr::bind_cols(outcome, AdditionalMetrics)
      }
    } else {
      outcome <- DataTable
    }
    return(outcome)
  }




#' @title simulate4GExperimentData
#' @description The function returns a four group data set based on one of four different distributions.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export simulate4GExperimentData
#' @param mean The mean (or rate for gamma data) of the baseline distribution
#' @param sd The standard deviation (or shape for gamma data) of the baseline distribution
#' @param diff The adjustment to the baseline mean for the alternative distribution.
#' @param GroupSize An integer defining the number of data items in each group.
#' @param type A string identifying the distrubtion used to simulate the data: 'n' for normal, 'l' for log-normal, 'g' for gamma, 'lap' for Laplace.
#' @param ExpAdj An additional adjument factor that is added to both the mean values. Defaults to zero.
#' @param StdAdj An aditional adjustment factor that is added to the second group variance (or rate for gamma data). Defaults to zero.
#' @param BlockEffect An additional factor that is added to the mean of the second group groups (shape for the gamma distribution). Defaults to zero.
#' @param BlockStdAdj An additional factor that is added to the variance of the second group (shape for the gamma distribution). Defaults to zero.
#' @return A table with four columns (BaselineData.B1, AlternativeData.B1,BaselineData.B2, AlternativeData.B2,) holding the data for each group and block. For lognormal data an additional four columns are added which return the log transformed data for each group.
#' @examples
#' set.seed(246)
#' simulate4GExperimentData(mean = 0, sd = 1, diff = 0.5, GroupSize = 5,
#'   type = "n", ExpAdj = 0, StdAdj = 0, BlockEffect = 0.5, BlockStdAdj = 0)
#' # A tibble: 5 x 4
#' #  BaselineData.B1 AlternativeData.B1 BaselineData.B2 AlternativeData.B2
#' #            <dbl>              <dbl>           <dbl>              <dbl>
#' # 1          0.533               1.84            0.749              3.98
#' # 2          0.251               2.03            1.56               1.09
#' # 3         -0.290               0.929           0.213              3.94
#' # 4         -1.48                1.17            1.13               0.106
#' # 5          0.0340              0.895           0.399              0.879

#' as.data.frame(
#'   simulate4GExperimentData(
#'     mean=0, sd=1, diff=0.5, GroupSize=5, type='l', ExpAdj=0, StdAdj=0,
#'     BlockEffect = 0.5, BlockStdAdj = 0))
#' #  BaselineData.B1 AlternativeData.B1 transBaselineData.B1 transAlternativeData.B1
#' #1       1.4019869           1.049158            0.3378905               0.0479875
#' #2       3.8514120           0.769227            1.3484398              -0.2623692
#' #3       6.5162726           1.574126            1.8743025               0.4537002
#' #4       1.3309218           1.082774            0.2858718               0.0795259
#' #5       0.2772234           1.630194           -1.2829316               0.4886992
#' #  BaselineData.B2 AlternativeData.B2 transBaselineData.B2 transAlternativeData.B2
#' #1       5.4656049          4.6095688            1.6984748               1.5281343
#' #2       1.6149559          2.0244244            0.4793077               0.7052854
#' #3       1.7718620          0.5504016            0.5720310              -0.5971070
#' #4       0.6774067          1.5434812           -0.3894834               0.4340404
#' #5       0.4507284          5.4987830           -0.7968903               1.7045268


simulate4GExperimentData <-
  function(mean,
           sd,
           diff,
           GroupSize,
           type = "n",
           ExpAdj = 0,
           StdAdj = 0,
           BlockEffect = 0,
           BlockStdAdj = 0) {
    Output1 <- simulate2GExperimentData(
      mean = mean,
      sd = sd,
      diff = diff,
      GroupSize = GroupSize,
      type = type,
      ExpAdj = ExpAdj,
      StdAdj = StdAdj,
      BlockEffect = 0,
      BlockStdAdj = 0
    )

    Output1 <-
      reshape::rename(
        Output1,
        c(BaselineData = "BaselineData.B1", AlternativeData = "AlternativeData.B1")
      )

    if (type == "l") {
      Output1 <-
        reshape::rename(
          Output1,
          c(
            transBaselineData = "transBaselineData.B1",
            transAlternativeData = "transAlternativeData.B1"
          )
        )
    }

    Output2 <- simulate2GExperimentData(
      mean = mean,
      sd = sd,
      diff = diff,
      GroupSize = GroupSize,
      type = type,
      ExpAdj = ExpAdj,
      StdAdj = StdAdj,
      BlockEffect = BlockEffect,
      BlockStdAdj = BlockStdAdj
    )

    Output2 <-
      reshape::rename(
        Output2,
        c(BaselineData = "BaselineData.B2", AlternativeData = "AlternativeData.B2")
      )

    if (type == "l") {
      Output2 <-
        reshape::rename(
          Output2,
          c(
            transBaselineData = "transBaselineData.B2",
            transAlternativeData = "transAlternativeData.B2"
          )
        )
    }


    OutputData <- dplyr::bind_cols(Output1, Output2)

    return(OutputData)
  }



#' @title calculateLargeSampleRandomizedBlockDesignEffectSizes
#' @description The function uses a simulates a large experiment  to estimate the asymptotic values of the probability of superiority, Cliff's d and the standardized mean difference data for a four group randomized blocks experiment for four different distributions: Normal (i.e. type='n'), log-normal (i.e. type='l'), gama (i.e. type='g') and Laplace (i.e., type='lap').
#' @author Barbara Kitchenham and Lech Madeyski
#' @export calculateLargeSampleRandomizedBlockDesignEffectSizes
#' @param meanC to act as the mean of the distribution (default 0) used to
#' generate the control group data (note for the gamma distribution this is the
#' rate parameter and must not be zero)
#' @param sdC the variance/spread of the distribution (default 1) used to
#' generate the control group data.
#' @param diff a value added to meanC to generate the treatment group data
#' (default 0).
#' @param N the size of each group (default 5000000)
#' @param type the distribution of the data to be generated. One of: 'n' for
#' normal (default), 'l' for log-normal, 'g' for gamma, and 'lap' for Laplace.
#' @param Blockmean, a value that can be added one of the blocks to represent
#' a fixed block effect (default 0).
#' @param StdAdj a value that can be added to sdC to introduce heterogeneity
#' into the treatment group (default 0).
#' @return A tibble identifying the sample statistics and the values of the
#' probability of superiority, Cliff's d and StdMD (labelled StdES)
#' @examples
#' set.seed=400
#' calculateLargeSampleRandomizedBlockDesignEffectSizes(
#'   meanC=0, sdC=1, diff=.5, N=100000, type='n',Blockmean=0.5,StdAdj = 0)
#' #  MeanC   SdC MeanT   SdT    BE  Phat Cliffd   UES   Var StdES
#' #  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl>
#' #1     0     1   0.5     1   0.5 0.638  0.277 0.501 0.998 0.502
#'
calculateLargeSampleRandomizedBlockDesignEffectSizes <-
  function(meanC = 0,
           sdC = 1,
           diff,
           N = 5e+06,
           type = "n",
           Blockmean = 0,
           StdAdj = 0) {
    MeanT <- meanC + diff
    sdT <- sdC + StdAdj

    if (type == "g") {
      MeanT <- meanC + diff + StdAdj
      sdT <- sdC + Blockmean
    }
    DataSet <- simulate4GExperimentData(
      mean = meanC,
      sd = sdC,
      diff = diff,
      GroupSize = N,
      type = type,
      ExpAdj = 0,
      StdAdj = StdAdj,
      BlockEffect = Blockmean,
      BlockStdAdj = 0
    )


    UES <-
      (
        base::mean(DataSet$AlternativeData.B1) + base::mean(DataSet$AlternativeData.B2) - base::mean(DataSet$BaselineData.B1) - base::mean(DataSet$BaselineData.B2)
      ) / 2

    Var <-
      (
        stats::var(DataSet$AlternativeData.B1) + stats::var(DataSet$AlternativeData.B2) + stats::var(DataSet$BaselineData.B1) + stats::var(DataSet$BaselineData.B2)
      ) / 4
    StdES <- UES / sqrt(Var)

    B1res <-
      calculatePhat(DataSet$BaselineData.B1, DataSet$AlternativeData.B1)
    B2res <-
      calculatePhat(DataSet$BaselineData.B2, DataSet$AlternativeData.B2)

    Phat <- (B1res$phat + B2res$phat) / 2

    Cliffd <- (B1res$dhat + B2res$dhat) / 2

    res <- as.data.frame(
      cbind(
        MeanC = meanC,
        SdC = sdC,
        MeanT = MeanT,
        SdT = sdT,
        BE = Blockmean,
        Phat = Phat,
        Cliffd = Cliffd,
        UES = UES,
        Var = Var,
        StdES = StdES
      )
    )
    res <- tibble::as_tibble(res)
    return(res)
  }

#' @title  Calc4GroupNPStats
#' @description This function does a non-parametric analysis of a randomized
#' blocks experiment assuming 2 blocks and 2 treatment conditions.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export Calc4GroupNPStats
#' @param x1 is the data associated with  treatment A in one block 1
#' @param x2 is the data associated with treatment B in block 1
#' @param x3 is the data associated with treatment A in block 2
#' @param x4  is the data associated with treatment B in block 2
#' @param sigfig is the number of significant digits in the data. If >0
#' the datav will be appropriately truncated.
#' @param alpha is the significance level for all statistical tests
#' @param alternative The type of statistical test. Valid values are one of
#' c('two.sided', 'greater', 'less')
#' @return The function returns the Cliff's d and its variance, the probability
#' of superiority, phat, and its variance for the 4 group experiment experiment.
#' @examples
#' set.seed(123)
#' x <- list()
#' x[[1]] <- rnorm(10, 0, 1)
#' x[[2]] <- rnorm(10, 0.8, 1)
#' x[[3]] <- rnorm(10, 0.5, 1)
#' x[[4]] <- rnorm(10, 1.3, 1)
#' as.data.frame(
#' Calc4GroupNPStats(x[[1]], x[[2]], x[[3]], x[[4]], sigfig = -1, alpha = 0.05)
#' )
#' #   N phat    phat.var  phat.df phat.test  phat.pvalue phat.sig phat.ci.upper
#' # 1 40 0.17 0.004966667 31.00131 -4.682539 5.324252e-05     TRUE     0.3137336
#' #  phat.ci.lower     d       vard d.sig d.ci.lower d.ci.upper        cor       sqse
#' # 1    0.02626639 -0.66 0.02060121  TRUE -0.8545073 -0.3031667 -0.3473684 0.01315789
#' #        ctvar n1 n2 sigCVt
#' # 1 0.005990797 20 20   TRUE
#' as.data.frame(
#' Calc4GroupNPStats(x[[1]], x[[2]], x[[3]], x[[4]], sigfig = -1, alpha = 0.05,
#' alternative = "less")
#' )
#' #   N phat    phat.var  phat.df phat.test  phat.pvalue phat.sig phat.ci.upper
#' # 1 40 0.17 0.004966667 31.00131 -4.682539 2.662126e-05     TRUE     0.2894908
#' #  phat.ci.lower     d       vard d.sig d.ci.lower d.ci.upper        cor       sqse
#' # 1             0 -0.66 0.02060121  TRUE         -1 -0.3677704 -0.3473684 0.01315789
#' #        ctvar n1 n2 sigCVt
#' # 0.005990797 20 20   TRUE
#' as.data.frame(
#' Calc4GroupNPStats(x[[2]], x[[1]], x[[4]], x[[3]], sigfig = -1, alpha = 0.05,
#' alternative = "greater")
#' )
#' #   N phat    phat.var  phat.df phat.test  phat.pvalue phat.sig phat.ci.upper
#' # 1 40 0.83 0.004966667 31.00131  4.682539 2.662126e-05     TRUE             1
#' #  phat.ci.lower    d       vard d.sig d.ci.lower d.ci.upper       cor       sqse
#' # 1     0.7105092 0.66 0.02060121  TRUE  0.3677704          1 0.3473684 0.01315789
#' #        ctvar n1 n2 sigCVt
#' # 1 0.005990797 20 20   TRUE
#'
#' #as.data.frame(
#' #Calc4GroupNPStats(x[[1]],x[[2]],x[[3]],x[[4]],sigfig=-1,alpha=0.00))
#' #Error in testfunctionParameterChecks(alternative = alternative, alpha = alpha,  :
#' #  Invalid alpha parameter, select alpha in range (0.0001,0.2)
#'
Calc4GroupNPStats <-
  function(x1,
           x2,
           x3,
           x4,
           sigfig = -1,
           alpha = 0.05,
           alternative = "two.sided") {
    # Check the significant digits to ensure that equal values are properly detected
    if (sigfig > 0) {
      x1 <- signif(x1, sigfig)
      x2 <- signif(x2, sigfig)
      x3 <- signif(x3, sigfig)
      x4 <- signif(x4, sigfig)
    }
    # Set up a dummy variable such that the observations using treatment A in block 1 are associated with the value 1 and observations using
    # treatment B in block 1 are associated with the value 0
    dummy1 <- c(rep(1, length(x1)), rep(0, length(x2)))
    # Concatenate the observations in block 1
    xCO1 <- c(x1, x2)
    # Use Wilcox's function to find tau and its two variances for block 1
    tau1 <- calculateKendalltaupb(xCO1, dummy1)
    n1 <- length(x1) + length(x2)
    # Set up a dummy variable such that the observations using treatment A in block 2 are associated with the value 1 and observations using
    # treatment B in block 2 are associated with the value 0

    dummy2 <- c(rep(1, length(x3)), rep(0, length(x4)))
    # Concatenate the observations in block 2
    xCO2 <- c(x3, x4)

    # Use Wilcox's function to find tau and its two variances for block 2

    tau2 <- calculateKendalltaupb(xCO2, dummy2)
    n2 <- length(x3) + length(x4)
    N <- n1 + n2

    average.tau <- (tau1$cor + tau2$cor) / 2

    ctvar <- (tau1$consistentvar + tau2$consistentvar) / 4

    # Find the confidence limits on the combined tau using t-distribution
    testfunctionParameterChecks(alternative = alternative,
                                alpha = alpha,
                                stderr = sqrt(ctvar))
    taupb.Analysis <- calcEffectSizeConfidenceIntervals(
      effectsize = average.tau,
      effectsize.variance = ctvar,
      effectsize.df = N - 6,
      alpha = 0.05,
      alternative = alternative,
      UpperValue = n1 / (2 * (n1 - 1)),
      LowerValue = -n1 / (2 * (n1 - 1))
    )

    # Find the average d and the combined variances for the full experiment

    d <- (calculateCliffd(x1, x2)$d + calculateCliffd(x3, x4)$d) / 2
    vard <-
      (calculateCliffd(x1, x2)$sqse.d + calculateCliffd(x3, x4)$sqse.d) / 4

    testfunctionParameterChecks(alternative = alternative,
                                alpha = alpha,
                                stderr = sqrt(vard))

    Cliffd.Analysis <-
      calcCliffdConfidenceIntervals(
        d.value = d,
        d.variance = vard,
        d.df = N - 4,
        alpha = alpha,
        alternative = alternative
      )


    # Find the average phat and the combined variances for the full experiment

    B1.BMP <- calculatePhat(x2, x1)
    B2.BMP <- calculatePhat(x4, x3)
    phat1 <- B1.BMP$phat
    phat2 <- B2.BMP$phat
    phat <- (phat1 + phat2) / 2
    se1 <- B1.BMP$s.e.
    se2 <- B2.BMP$s.e.

    phat.se <- sqrt((B1.BMP$s.e ^ 2 + B2.BMP$s.e. ^ 2) / 4)
    phat.var <- phat.se ^ 2
    phat.df <- B1.BMP$df + B2.BMP$df

    testfunctionParameterChecks(
      alternative = alternative,
      alpha = alpha,
      stderr = sqrt(phat.var)
    )
    phat.Analysis <-
      calcPHatConfidenceIntervals(
        phat = phat,
        phat.variance = phat.var,
        phat.df = phat.df,
        alpha = alpha,
        alternative = alternative
      )


    output <- tibble::tibble(
      N = N,
      phat = phat,
      phat.var = phat.var,
      phat.df = phat.df,
      phat.test = phat.Analysis$phat.test,
      phat.pvalue = phat.Analysis$pvalue,
      phat.sig = phat.Analysis$phat.sig,
      phat.ci.upper = phat.Analysis$phat.ci.upper,
      phat.ci.lower = phat.Analysis$phat.ci.lower,
      d = d,
      vard = vard,
      d.sig = Cliffd.Analysis$d.sig,
      d.ci.lower = Cliffd.Analysis$d.ci.lower,
      d.ci.upper = Cliffd.Analysis$d.ci.upper,
      cor = average.tau,
      ctvar = ctvar,
      n1 = n1,
      n2 = n2,
      sigCVt = taupb.Analysis$ES.sig
    )

    return(output)
  }


#' @title RandomizedBlocksAnalysis
#' @description The function performs a heteroscedastic test of a two treatment by J blocks randomized blocks effect size. The data are assumed to be stored in $x$ in list mode. All groups are assumed to be independent. Missing values are not permitted.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export RandomizedBlocksAnalysis
#' @param x the structure holding the data. In list format, for a 2 treatment by J block randomized blocks experiments, there are  2J list elements each one specifying the outcome for a specific block and a specific treatment.
#' @param con is a 2J list containing the contrast coefficients that are used to calculate the mean effect size.
#' @param alpha is the Type 1 error level used for the test of significance
#' (default 0.05)
#' @param alternative The type of statistical test. Valid values are one of
#' c('two.sided', 'greater', 'less')
#' @return The t-test and its associated metrics (i.e., critical value standard error and degrees of freedom) and the estimate of the contrast with its upper and lower confidence interval bounds and p-value.
#' @examples
#' set.seed(123)
#' x <- list()
#' x[[1]] <- rnorm(10, 0, 1)
#' x[[2]] <- rnorm(10, 0.8, 1)
#' x[[3]] <- rnorm(10, 0.5, 1)
#' x[[4]] <- rnorm(10, 1.3, 1)
#' vec <- c(-1, 1, -1, 1) / 2
#' RandomizedBlocksAnalysis(x, con = vec, alpha = 0.05)
#' # $n
#' # [1] 10 10 10 10
#' # $test
#' #      test     crit        se       df
#' # [1,] 4.432644 2.038622 0.2798104 31.33793
#' # $psihat
#' #      psihat  ci.lower ci.upper      p.value
#' # [1,] 1.2403 0.6698721 1.810728 0.0001062952
#' # $sig
#' # [1] TRUE

#' RandomizedBlocksAnalysis(x,con=vec,alpha=0.05,alternative='greater')
#' # n
#' # [1] 10 10 10 10

#' # $test
#' #          test     crit        se       df
#' # [1,] 4.432644 1.694956 0.2798104 31.33793

#' # $psihat
#' # psihat  ci.lower ci.upper      p.value
#' #[1,] 1.2403 0.7660336      Inf 5.314762e-05

#' # $sig
#' # [1] TRUE

#' RandomizedBlocksAnalysis(x,con=-vec,alpha=0.05,alternative='greater')
#' #$n
#' #[1] 10 10 10 10

#' #$test
#' #          test     crit        se       df
#' #[1,] -4.432644 1.694956 0.2798104 31.33793

#' #$psihat
#' #      psihat  ci.lower ci.upper   p.value
#' #[1,] -1.2403 -1.714566      Inf 0.9999469

#' #$sig
#' #[1] FALSE


#' x[[5]]=rnorm(10,-0.2,1)
#' x[[6]]=rnorm(10,0.6,1)
#' vec=c(1,-1,1,-1,1,-1)/3
#' RandomizedBlocksAnalysis(x,con=vec,alpha=0.05,alternative='less')
#' #$n
#' #[1] 10 10 10 10 10 10

#' #$test
#' #          test     crit       se       df
#' #[1,] -4.946987  1.677021 0.236575 48.29776

#' #$psihat
#' #        psihat ci.lower   ci.upper     p.value
#' #[1,] -1.170334     -Inf -0.7735925 4.76961e-06

#' #$sig
#' #[1] TRUE

RandomizedBlocksAnalysis <-
  function(x,
           con = c(-0.5, 0.5, -0.5, 0.5),
           alpha = 0.05,
           alternative = "two.sided") {
    if (!base::is.list(x)) {
      stop("Data must be stored list mode.")
    }
    con <- base::as.matrix(con)

    if (ncol(con) > 1) {
      stop("Only one linear contrast permitted")
    }
    J <- length(x)
    sam <- NA
    h <- base::vector("numeric", J)
    w <- base::vector("numeric", J)
    xbar <- base::vector("numeric", J)
    for (j in 1:J) {
      if (sum(as.numeric(is.na(x[[j]]))) > 0) {
        stop("Missing values are not permitted")
      }
      sam[j] <- length(x[[j]])
      h[j] <- length(x[[j]])
      # h is the number of observations in the jth group.
      w[j] <-
        ((length(x[[j]]) - 1) * stats::var(x[[j]])) / (h[j] * (h[j] - 1)) # The variance of the jth group
      xbar[j] <- base::mean(x[[j]]) # The mean of the jth group
    }

    if (nrow(con) != length(x)) {
      stop("The number of groups does not match the number of contrast coefficients.")
    }
    psihat <- base::matrix(0, 1, 4)
    dimnames(psihat) <-
      list(NULL, c("psihat", "ci.lower", "ci.upper", "p.value"))
    test <- matrix(0, 1, 4)
    dimnames(test) <- list(NULL, c("test", "crit", "se", "df"))
    df <- 0


    psihat[1, 1] <- sum(con[, 1] * xbar)
    sejk <-
      sqrt(sum(con[, 1] ^ 2 * w)) # The pooled standard error of the contrast

    test[1, 1] <-
      sum(con[, 1] * xbar) / sejk # The value of the t-test
    df <-
      (sum(con[, 1] ^ 2 * w)) ^ 2 / sum(con[, 1] ^ 4 * w ^ 2 / (h - 1)) # Degrees of freedom allowing for heterogeneity
    test[1, 3] <- sejk
    test[1, 4] <- df

    if (alternative == "two.sided") {
      vv <- (1 - alpha / 2)
      crit <-
        stats::qt(vv, df) # The critical value of the t-test for the degrees of freedom
      test[1, 2] <- crit

      psihat[1, 2] <- psihat[1, 1] - crit * sejk
      psihat[1, 3] <- psihat[1, 1] + crit * sejk
      psihat[1, 4] <- 2 * (1 - stats::pt(abs(test[1, 1]), df))

      sig <- psihat[1, 2] > 0 | psihat[1, 3] < 0
    } else {
      vv <- (1 - alpha)
      crit <- stats::qt(vv, df)
      test[1, 2] <- crit
      stats::qt(vv, df)
      lowerbound <- psihat[1, 1] - crit * sejk
      upperbound <- psihat[1, 1] + crit * sejk

      if (alternative == "greater") {
        # The upper bound is irrelevant
        upperbound <- Inf
        sig <- lowerbound > 0
        psihat[1, 4] <- (1 - stats::pt((test[1, 1]), df))
      } else {
        # The lower bound is irrelevant
        lowerbound <- -Inf
        sig <- upperbound < 0
        psihat[1, 4] <- (stats::pt((test[1, 1]), df))
      }
      psihat[1, 2] <- lowerbound
      psihat[1, 3] <- upperbound
    }


    list(
      n = sam,
      sig = sig,
      test = test,
      psihat = psihat
    )
  }



#' @title simulateRandomizedBlockDesignEffectSizes
#' @description This simulates a two-block and two-treatment design based on one of four distributions, and finds the values of ktau and Cliffs d and their variances. It simulates a randomised blocks experiment with two treatment groups and two control groups each of which being divided into two blocks. It assumes equal group sizes but  group spread (standard deviation can be changed, see StAdj). It returns values of both parametric and non-parametric effect sizes and their variance and significance. For the logarithmic distribution it calculates effect sizes based on the log transformed data as well as the raw data.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export simulateRandomizedBlockDesignEffectSizes
#' @param mean The default value for all groups which can be changed for the two treatment groups using the parameter diff and for the two block 2 groups using the parameter Blockmean
#' @param sd The default spread used for all four groups unless adjusted by the StdAdj. It must be a real value greater than 0.
#' @param diff This is added to the parameter mean to obtain the required mean for treatment groups. It can be a real value and can take the value zero.
#' @param N this is the number of observations in each group. It must be an integer greater than 3.
#' @param type this specifies the underlying distribution used to generate the data. it takes the values 'n' for a normal distribution, 'l' for lognormal distribution,'g' for a gamma distribution, 'lap' for a Laplace distribution.
#' @param alpha  The level used for statistical tests (default 0.05).
#' @param Blockmean if >0 an adjustment made to both group means in Block 2
#' @param BlockStdAdj if >0, an adjustment that can be made to the sd of each group in block 2
#' @param StdAdj this specifies the extent of variance instability introduced by the treatment and if >0 will be used to amend the sd parameter for both treatment groups. This value must be positive and less than 0.5
#' @param AlwaysTwoSidedTests Logical varable (default FALSE) if TRUE the
#' function always performs two-sided tests. Otherwise if the parameter diff is
#' not equal to zero, the function performs one-sided tests.
#' @param ReturnData Logical variable, If TRUE, the function simply returns the
#' generated data. If false (which is default value) the function returns
#' various effect sizes and whether the effect sizes are statistically
#' significant.
#' @return data frame incl. either the non-parametric and parametric effect
#' sizes and whether the effect sizes are significant at the 0.05 level or the
#' generated data depending on the value of the ReturnData parameter.
#' @examples
#' set.seed(123)
#' as.data.frame(
#'   simulateRandomizedBlockDesignEffectSizes(
#'     mean = 0, sd = 1, diff = .5, N = 10, type = "n", alpha = 0.05,
#'     Blockmean = 0.5, BlockStdAdj = 0, StdAdj = 0))
#' #  N phat    phat.var  phat.df phat.test  phat.pvalue phat.sig phat.ci.upper phat.ci.lower    d
#' # 1 40 0.79 0.005866667 30.15715  3.786189 0.0003403047     TRUE             1     0.6600213 0.58
#' #     vard d.sig d.ci.lower d.ci.upper       cor       sqse       ctvar n1 n2 sigCVt sigCVn
#' # 1 0.02430788  TRUE  0.2775601          1 0.3052632 0.01315789 0.006953352 20 20   TRUE   TRUE
#' #    ttest.sig     ES  Variance   StdES BlockEffect MedianDiff
#' # 1      TRUE 0.9402999 0.7829385 1.06268    0.307119   1.313642
#' set.seed(123)
#' as.data.frame(
#'   simulateRandomizedBlockDesignEffectSizes(
#'     mean = 0, sd = 1, diff = 0.5, N = 10, type = "n", alpha = 0.05,
#'     Blockmean = 0.5, BlockStdAdj = 0, StdAdj = 0, AlwaysTwoSidedTests = TRUE)
#'     )
#' #   N phat    phat.var  phat.df phat.test  phat.pvalue phat.sig phat.ci.upper phat.ci.lower
#' # 1 40 0.79 0.005866667 30.15715  3.786189 0.0006806094     TRUE      0.946392      0.633608
#' #     d       vard d.sig d.ci.lower d.ci.upper       cor       sqse       ctvar n1 n2 sigCVt
#' # 1 0.58 0.02430788  TRUE  0.2135334  0.8033737 0.3052632 0.01315789 0.006953352 20 20   TRUE
#' #  ttest.sig        ES  Variance   StdES BlockEffect MedianDiff
#' # 1      TRUE 0.9402999 0.7829385 1.06268    0.307119   1.313642
#' set.seed(123)
#' as.data.frame(
#'   simulateRandomizedBlockDesignEffectSizes(
#'     mean = 0, sd = 1, diff = .5, N = 10, type = "l", alpha = 0.05,
#'     Blockmean = 0.5, BlockStdAdj = 0, StdAdj = 0, ReturnData = TRUE))
#' #   BaselineData.B1 AlternativeData.B1 transBaselineData.B1 transAlternativeData.B1
#' # 1        0.5709374          5.6073700          -0.56047565              1.72408180
#' # 2        0.7943926          2.3627208          -0.23017749              0.85981383
#' # 3        4.7526783          2.4615013           1.55870831              0.90077145
#' # 4        1.0730536          1.8416883           0.07050839              0.61068272
#' # 5        1.1380175          0.9456894           0.12928774             -0.05584113
#' # 6        5.5570366          9.8445021           1.71506499              2.28691314
#' # 7        1.5855260          2.7124451           0.46091621              0.99785048
#' # 8        0.2822220          0.2307046          -1.26506123             -1.46661716
#' # 9        0.5031571          3.3246217          -0.68685285              1.20135590
#' # 10       0.6404002          1.0275821          -0.44566197              0.02720859
#' #   BaselineData.B2 AlternativeData.B2 transBaselineData.B2 transAlternativeData.B2
#' # 1        0.5667575           4.163950           -0.5678237               1.4264642
#' # 2        1.3258120           2.023702            0.2820251               0.7049285
#' # 3        0.5909615           6.653384           -0.5260044               1.8951257
#' # 4        0.7954150           6.541284           -0.2288912               1.8781335
#' # 5        0.8824622           6.181624           -0.1250393               1.8215811
#' # 6        0.3052289           5.412117           -1.1866933               1.6886403
#' # 7        3.8106015           4.729964            1.3377870               1.5539177
#' # 8        1.9220131           2.555092            0.6533731               0.9380883
#' # 9        0.5282757           2.001781           -0.6381369               0.6940373
#' # 10       5.7765980           1.858053            1.7538149               0.6195290
simulateRandomizedBlockDesignEffectSizes <-
  function(mean,
           sd,
           diff,
           N,
           type = "n",
           alpha = 0.05,
           Blockmean = 0,
           BlockStdAdj = 0,
           StdAdj = 0,
           AlwaysTwoSidedTests = FALSE,
           ReturnData = FALSE) {
    # 29-06-2022 Changed to support one-sided tests 14-03-2022 BlockStdadj parameter changed to BlockStdAdj for consistency with other functions.
    # Changes to calls to generate data chnaged to cater for new label.  Generate data. x and x2 hold control data, y and y2 to hold treatment
    # data.  x and y are ib block 1 and x2 and y2 are in block 2

    alpha.set <- alpha

    # Allowed for one-sided tests
    if (AlwaysTwoSidedTests | (diff == 0)) {
      alternative <- "two.sided"
    } else {
      # Find direction of MD

      if (type == "g") {
        # A negative value for diff decreases the rate & the MD is positive
        Positive.MD <- (diff < 0)
      } else {
        Positive.MD <- (diff > 0)
      }

      if (Positive.MD) {
        alternative <- "greater"
      } else {
        alternative <- "less"
      }
    }

    DataSet <- simulate4GExperimentData(
      mean = mean,
      sd = sd,
      diff = diff,
      GroupSize = N,
      type = type,
      ExpAdj = 0,
      StdAdj = StdAdj,
      BlockEffect = Blockmean,
      BlockStdAdj = 0
    )

    if (ReturnData) {
      res <- DataSet
    } else {
      res <-
        Calc4GroupNPStats(
          DataSet$AlternativeData.B1,
          DataSet$BaselineData.B1,
          DataSet$AlternativeData.B2,
          DataSet$BaselineData.B2,
          alpha = alpha.set,
          alternative = alternative
        )

      StandardMetrics <- as.data.frame(res)

      # Add the results of a t-test as a baseline
      UES <-
        (
          base::mean(DataSet$AlternativeData.B1) + base::mean(DataSet$AlternativeData.B2) - base::mean(DataSet$BaselineData.B1) - base::mean(DataSet$BaselineData.B2)
        ) / 2

      BlockEffect <-
        (
          base::mean(DataSet$BaselineData.B2) + base::mean(DataSet$AlternativeData.B2) - base::mean(DataSet$BaselineData.B1) - base::mean(DataSet$AlternativeData.B1)
        ) / 2

      MedDiff <-
        (
          stats::median(DataSet$AlternativeData.B1) + stats::median(DataSet$AlternativeData.B2) - stats::median(DataSet$BaselineData.B1) -
            stats::median(DataSet$BaselineData.B2)
        ) / 2



      Var <-
        (
          stats::var(DataSet$AlternativeData.B1) + stats::var(DataSet$AlternativeData.B2) + stats::var(DataSet$BaselineData.B1) + stats::var(DataSet$BaselineData.B2)
        ) / 4
      # Estimate of Cohen's d
      StdES <- UES / sqrt(Var)


      # Need to use the linear contrast method for the significance test of randomized blocks data that allows for variance heterogeneity.
      newlist <- list()
      newlist[[1]] <- DataSet$BaselineData.B1
      newlist[[2]] <- DataSet$AlternativeData.B1
      newlist[[3]] <- DataSet$BaselineData.B2
      newlist[[4]] <- DataSet$AlternativeData.B2
      vec <- c(-1, 1, -1, 1) / 2


      res.t <-
        RandomizedBlocksAnalysis(newlist,
                                 con = vec,
                                 alpha = alpha.set,
                                 alternative = alternative)

      ttest.sig <- as.logical(res.t$sig)


      StandardMetrics <- dplyr::bind_cols(
        StandardMetrics,
        ttest.sig = ttest.sig,
        ES = UES,
        Variance = Var,
        StdES = StdES,
        BlockEffect = BlockEffect,
        MedianDiff = MedDiff
      )


      if (type == "l") {
        loglist <- list()

        loglist[[1]] <- DataSet$transBaselineData.B1
        loglist[[2]] <- DataSet$transAlternativeData.B1
        loglist[[3]] <- DataSet$transBaselineData.B2
        loglist[[4]] <- DataSet$transAlternativeData.B2

        vec <- c(-1, 1, -1, 1) / 2
        logres.t <-
          RandomizedBlocksAnalysis(loglist,
                                   con = vec,
                                   alpha = alpha.set,
                                   alternative = alternative)
        Log.sig <- as.logical(logres.t$sig)

        ES.trans <-
          (
            base::mean(DataSet$transAlternativeData.B1) + base::mean(DataSet$transAlternativeData.B2) - base::mean(DataSet$transBaselineData.B1) -
              base::mean(DataSet$transBaselineData.B2)
          ) / 2

        VarTrans <-
          (
            stats::var(DataSet$transBaselineData.B1) + stats::var(DataSet$transBaselineData.B2) + stats::var(DataSet$transAlternativeData.B1) +
              stats::var(DataSet$transAlternativeData.B2)
          ) / 4

        StdES.trans <- ES.trans / sqrt(VarTrans)

        AdditionalMetrics <-
          list(
            Log.sig = Log.sig,
            ES.Trans = ES.trans,
            StdES.Trans = StdES.trans,
            VarTrans = VarTrans
          )

        AdditionalMetrics <- tibble::as_tibble(AdditionalMetrics)
        StandardMetrics <-
          dplyr::bind_cols(StandardMetrics, AdditionalMetrics)
      }
      res <- StandardMetrics
    }
    return(res)
  }



#' title RandomizedBlocksExperimentSimulations
#' description This function performs multiple simulations of 4 group balanced randomised Block experiments with two control groups and two treatment groups where one control group and one treatment group are assigned to block 1 and the other control group and treatment group are assigned to block 2.  The simulations are based on one of four distributions and a specific group size. The function identifies the average value of the non-parametric effect sizes P-hat, Cliff' d and their variances and whether ot not the statistics were significant at the 0.05 level. We also present the values of the t-test as a comparison.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export RandomizedBlocksExperimentSimulations
#' @param mean The default mean for all 4 groups. The default for the two treatment groups can be altered using the parameter diff and the block mean for block 2 can be altered using the parameter Blockmean.
#' @param sd The default spread for all 4 groups. It must be a real value greater than 0. If can be altered for treatment groups using the parameter StdAdj and for Block 2 groups using BlockStdAdj
#' @param diff The is is added to the parameter mean, to define the mean of the other treatment group. It can be a real value ad can take the value zero.
#' @param N this is the number of observations in each group. It must be an integer greater than 3.
#' @param reps this identifies the number of times the simulation is replicated.
#' @param type this specifies the underlying distribution used to generate the data. it takes the values 'n' for a normal distribution, 'l' for lognormal distribution,'g' for a gamma distribution, 'lap' for a Laplace distribution.
#' @param alpha is the Type 1 error level used for constructing confidence
#' intervals and statistical tests (default 0.05)
#' @param Blockmean is the effect of having two different blocks
#' @param BlockStdAdj is the variance associated with the Block mean. If Blockvar is zero it means we are treat the block effect as a fixed effect. If BlockStdAdj>0, we treat the block effect as a random effect.
#' @param StdAdj The value used to introduce heterogeneity into the treatment groups variance if required.
#' @param seed this specifies the seed value for the simulations and allows the experiment to be repeated.
#' @param returnData if TRUE the function returns the generated data otherwise it returns summary statistics.
#' @param AlwaysTwoSidedTests A boolean variable. If TRUE the simulations always used two-sided tests otherwise the simulations use one-sided tests.
#' return depending on the parameter returnData it returns the generated nonparametric and parametric values and their statistical significance (1 for significant, 0 for not significant) or the summary statistics (averages of effect sizes and their variances and the proportion significant effect sizes)
#' @examples
#' as.data.frame(
#'   RandomizedBlocksExperimentSimulations(
#'     mean = 0, sd = 1, diff = 0.5, N = 10, reps = 50, type = "n",
#'     alpha = 0.05, Blockmean = 0.5, BlockStdAdj = 0, StdAdj = 0, seed = 123,
#'     AlwaysTwoSidedTests = FALSE))
#' #     phat     varphat sigphat emp.phat.var      d      vard sigd  emp.d.var
#' #1 0.64415 0.008271389    0.45  0.005888917 0.2883 0.0340919 0.41 0.02355567
#' #        StdES        ES       Var emp.StdESvar   MedDiff tpower
#' #1   0.5413961 0.5264245 0.9904726   0.08811262 0.5538213   0.46
#' #as.data.frame(
#'  # RandomizedBlocksExperimentSimulations(
#'  #   mean = 0, sd = 1, diff = 0.5, N = 10, reps = 500, type = "n",
#'  #   alpha = 0.05, Blockmean = 0.5, BlockStdAdj = 0, StdAdj = 0, seed = 123,
#'  #   AlwaysTwoSidedTests = FALSE))
#' #  phat    varphat       sigphat emp.phat.var  d       vard        sigd  emp.d.var
#' # 1  0.63967  0.008322856  0.436   0.007728698   0.27934 0.03430328  0.416 0.03091479
#' #       StdES        ES      Var emp.StdESvar   MedDiff
#' # 1 0.5130732 0.5029075 1.001602    0.1116687 0.5110203
#' #  tpower
#' # 1   0.45
#'
#' #as.data.frame(
#'  # RandomizedBlocksExperimentSimulations(
#'  #   mean = 0, sd = 1, diff = 0.5, N = 10, reps = 500, type = "n",
#'  #   alpha = 0.05, Blockmean = 0.5, BlockStdAdj = 0, StdAdj = 0, seed = 123,
#'  #   AlwaysTwoSidedTests = TRUE))
#' #       phat     varphat sigphat emp.phat.var        d       vard   sigd
#' # 1  0.63967 0.008322856   0.326  0.007728698  0.27934 0.03430328  0.282
#' #     emp.d.var        StdES        ES      Var
#' # 1  0.03091479    0.5130732 0.5029075 1.001602
#' # emp.StdESvar   MedDiff tpower
#' # 1    0.1116687 0.5110203  0.334
#'
#' #RandomizedBlocksExperimentSimulations(
#'  # mean = 0, sd = 1, diff = 0.5, N = 10, reps = 10, type = "n", alpha = 0.05,
#'  #Blockmean = 0.5, BlockStdAdj = 0, StdAdj = 0, seed = 123, returnData = TRUE)
#' # A tibble: 10 x 6
#' #   Cliffd  PHat StdES CliffdSig PHatSig ESSig
#' #    <dbl> <dbl> <dbl>     <dbl>   <dbl> <dbl>
#' # 1   0.58 0.79  1.06          1       1     1
#' # 2   0.21 0.605 0.383         0       0     0
#' # 3   0.37 0.685 0.761         1       1     1
#' # 4   0.44 0.72  0.821         1       1     1
#' # 5   0.13 0.565 0.240         0       0     0
#' # 6   0.16 0.58  0.222         0       0     0
#' # 7   0.38 0.69  0.580         1       1     1
#' # 8   0.48 0.74  0.882         1       1     1
#' # 9   0.11 0.555 0.181         0       0     0
#' # 10  -0.03 0.485 0.124        0       0     0
#'
RandomizedBlocksExperimentSimulations <-
  function(mean,
           sd,
           diff,
           N,
           reps,
           type = "n",
           alpha = 0.05,
           Blockmean = 0,
           BlockStdAdj = 0,
           StdAdj = 0,
           seed = 123,
           returnData = FALSE,
           AlwaysTwoSidedTests = FALSE) {
    # 20-06-2022 Changed to cater for one-sided tests 14-03-2022 Changed the parameter name BlockStdadj to BlockStdAdj for consistency

    dsum <-
      0 # This is used to sum the value of Cliff's d across the replications
    dvarsum <-
      0 # This is used to aggregate the variance of d across the replications
    d.sig <-
      0 # This is used to sum the number of significant d values  across the replications
    d.ss <-
      0 # This is used to sum the squared value of Cliff's d across the replications. It can be used to calculate the empirical variance of the average d.

    phatsum <-
      0 # This is used to sum the value of phat across the replications
    phatvarsum <-
      0 # This is used to aggregate the variance of phat across the replications
    phat.sig <-
      0 # This is used to sum the number of significant phat values  across the replications
    phat.ss <-
      0 # This is used to sum the squared value of phat across the replications. It can be used to calculate the empirical variance of the average phat.

    tsig <-
      0 # This is used to count the number of significant t values across the replications
    ES <-
      0 # This is used to sum the value of the parametric effect size (unstandardized) across the replications
    StdES <-
      0 # This is used to sum the value of the parametric effect size (standardized) across the replications
    StdES.ss <-
      0 # This is used to sum the squared value of the parametric effect size (standardized) across the replications. It can be used to calculate the empirical variance of the overall mean value

    Var <-
      0 # This is used to sum the estimate of the pooled variance across replications
    MedDiff <- 0 # This is used to sum the median across replications

    ES.l.trans <-
      0 # This is used to sum the unstandardized effect size of the log-normal data after being transformed
    StdES.l.trans <-
      0 # This is used to sum the standardized effect size of the log-normal data after being transformed
    Var.l.trans <-
      0 # This is used to sum the variance of the log-normal data after being transformed

    trans.sig <-
      0 # This is used to sum the number of signficant t-tests of the transformed lognormal data

    base::set.seed(seed)

    DataTable <- NULL

    for (i in 1:reps) {
      # Call the program than generates the randomized block experiment data sets and calculates the sample statistics
      res <- simulateRandomizedBlockDesignEffectSizes(
        mean = mean,
        sd = sd,
        diff = diff,
        N = N,
        type = type,
        alpha = alpha,
        Blockmean = Blockmean,
        BlockStdAdj = BlockStdAdj,
        StdAdj = StdAdj,
        AlwaysTwoSidedTests = AlwaysTwoSidedTests
      )

      if (!returnData) {
        # Calculate the averages of the effect sizes across the replications Cliff's d
        dsum <- dsum + res$d
        dvarsum <- dvarsum + res$vard
        if (res$d.sig) {
          d.sig <- d.sig + 1
        }
        d.ss <- d.ss + res$d ^ 2

        # Probability of Speriority (phat)
        phatsum <- phatsum + res$phat
        if (res$phat.sig) {
          phat.sig <- phat.sig + 1
        }
        phat.ss <- phat.ss + res$phat ^ 2
        phatvarsum <- phatvarsum + res$phat.var

        # Standard parametric effect sizes
        ES <- ES + res$ES
        StdES <- StdES + res$StdES
        StdES.ss <- StdES.ss + res$StdES ^ 2
        Var <- Var + res$Variance
        MedDiff <- MedDiff + res$MedianDiff
        tsig <- tsig + if (res$ttest.sig) {
          1
        } else {
          0
        }

        if (type == "l") {
          trans.sig <- trans.sig + if (res$Log.sig) {
            1
          } else {
            0
          }

          ES.l.trans <- ES.l.trans + res$ES.Trans
          StdES.l.trans <- StdES.l.trans + res$StdES.Trans
          Var.l.trans <- Var.l.trans + res$VarTrans
        }
      } else {
        # Store the outcome from each replication 18-5-2022 Included significance in output
        DataTable <-
          tibble::tibble(dplyr::bind_rows(
            DataTable,
            dplyr::bind_cols(
              Cliffd = res$d,
              PHat = res$phat,
              StdES = res$StdES,
              CliffdSig = as.numeric(res$d.sig),
              PHatSig = as.numeric(res$phat.sig),
              ESSig = as.numeric(res$ttest.sig)
            )
          ))
      }
    }
    if (!returnData) {
      # Calculate averages.
      d <- dsum / reps
      vard <- dvarsum / reps
      sigd <- d.sig / reps
      emp.d.var <- (d.ss - reps * d ^ 2) / (reps - 1)

      phat <- phatsum / reps
      varphat <- phatvarsum / reps
      sigphat <- phat.sig / reps
      emp.phat.var <- (phat.ss - reps * phat ^ 2) / (reps - 1)

      ES <- ES / reps
      StdES <- StdES / reps
      emp.StdESvar <- (StdES.ss - reps * StdES ^ 2) / (reps - 1)
      Var <- Var / reps
      tpower <- tsig / reps
      MedDiff <- MedDiff / reps

      if (type == "l") {
        ESLog <- ES.l.trans / reps
        StdESLog <- StdES.l.trans / reps
        VarLog <- Var.l.trans / reps
        Log.sig <- trans.sig / reps
      }
    }
    if (!returnData) {
      if (type == "l") {
        outcome <- tibble::tibble(
          phat,
          varphat,
          sigphat,
          emp.phat.var,
          d,
          vard,
          sigd,
          emp.d.var,
          StdES,
          ES,
          Var,
          emp.StdESvar,
          MedDiff,
          tpower,
          ESLog,
          StdESLog,
          VarLog,
          Log.sig
        )
      } else {
        outcome <-
          tibble::tibble(
            phat,
            varphat,
            sigphat,
            emp.phat.var,
            d,
            vard,
            sigd,
            emp.d.var,
            StdES,
            ES,
            Var,
            emp.StdESvar,
            MedDiff,
            tpower
          )
      }
    } else {
      outcome <- DataTable
    }

    return(outcome)
  }


#' @title calculateNullESAccuracy
#' @description The function uses simulation to assess the accuracy when the mean difference is zero, and the type 1 error rates of parametric and non-parametric effect sizes for both two group randomized designs and four group randomized block designs, for each of four different distributions.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export calculateNullESAccuracy
#' @param mean The mean of the baseline distribution.
#' @param sd The standard deviation or shape of the baseline distribution
#' @param N The number of observations per group for two group experiments and N/2 the sample sizes for four group experiments. N must be even to ensure equal N/2 defines appropriate sample sizes per group for 4 group experiments
#' @param reps The number of replications (i.e. two-group and four group experiments) to be simulated
#' @param type A string parameter defining the distribution being simulated i.e. 'n' for normal data, 'l' for log-normal data, 'g' for gamma data and 'lap' for LaPlace data.
#' @param seed A starting value for the simulations
#' @param StdAdj A numerical parameter that can be used to add additional variance for normal, lognormal and Laplce data and to change the shape parameter for gamma data.
#' @param Blockmean A numerical parameter used to introduce a fixed Block effect for four group experiments
#' @return A tibble identifying the median absolute error for the effect sizes Cliff's d, phat and StdMD and the Type 1 error rate, estimated from the proportion of significant effect sizes in the simulated experiments.
#' @examples
#' as.data.frame(
#'   calculateNullESAccuracy(
#'     mean=0,sd=1,N=10,reps=30,type='n',seed=123,StdAdj = 0,Blockmean = 0.5))
#'#   Design Obs CliffdAbsError PHatAbsError StdESdAbsError  varCliffd    varPHat
#'# 1   2G_n  20           0.20         0.10      0.2624447 0.05530851 0.01382713
#'# 2   4G_n  20           0.16         0.08      0.1848894 0.05447540 0.01361885
#'#    varStdES    ObsCliffd   ObsPHat     ObsStdES CliffdType1ER PHatType1ER
#'# 1 0.1425374  0.021333333 0.5106667 0.0001190251             0           0
#'# 2 0.1484728 -0.009333333 0.4953333 0.0295002335             0           0
#'#   StdESType1ER
#'# 1   0.03333333
#'# 2   0.03333333
#' #as.data.frame(
#'  # calculateNullESAccuracy(
#'  #   mean=0,sd=1,N=10,reps=100,type='n',seed=123,StdAdj = 0,Blockmean = 0.5))
#' #  Design Obs CliffdAbsError PHatAbsError StdESdAbsError  varCliffd    varPHat  varStdES ObsCliffd
#' #1   2G_n  20           0.21        0.105      0.3303331 0.08064949 0.02016237 0.2488365   -0.0010
#' #2   4G_n  20           0.16        0.080      0.2565372 0.05933430 0.01483358 0.1769521    0.0052
#' #  ObsPHat    ObsStdES CliffdType1ER PHatType1ER StdESType1ER
#' #1  0.4995 -0.02395895          0.07        0.08         0.08
#' #2  0.5026  0.03769940          0.01        0.01         0.02

calculateNullESAccuracy <-
  function(mean = 0,
           sd = 1,
           N = 10,
           reps = 10,
           type = "n",
           seed = 123,
           StdAdj = 0,
           Blockmean = 0.5) {
    NullESAccuracyTable <- NULL

    Out1 <-
      RandomExperimentSimulations(
        mean = mean,
        sd = sd,
        diff = 0,
        N = N,
        reps = reps,
        type = type,
        seed = seed,
        StdAdj = StdAdj,
        returnData = TRUE
      )

    CliffdMdMRE <- median(abs(Out1$Cliffd))

    PHatMdMRE <- median(abs(Out1$PHat - 0.5))


    StdESMdMRE <- median(abs(Out1$StdES))

    Design <- "2G"
    Design <- paste(Design, type, sep = "_")


    NullESAccuracyTable <-
      tibble::tibble(dplyr::bind_rows(
        NullESAccuracyTable,
        dplyr::bind_cols(
          Design = Design,
          Obs = as.character(2 * N),
          CliffdAbsError = CliffdMdMRE,
          PHatAbsError = PHatMdMRE,
          StdESdAbsError = StdESMdMRE,
          varCliffd = var(Out1$Cliffd),
          varPHat = var(Out1$PHat),
          varStdES = var(Out1$StdES),
          ObsCliffd = mean(Out1$Cliffd),
          ObsPHat = mean(Out1$PHat),
          ObsStdES = mean(Out1$StdES),
          CliffdType1ER = mean(Out1$CliffdSig),
          PHatType1ER = mean(Out1$PHatSig),
          StdESType1ER = mean(Out1$ESSig)
        )
      ))


    Design <- "4G"
    Design <- paste(Design, type, sep = "_")

    Out1 <- RandomizedBlocksExperimentSimulations(
      mean = mean,
      sd = sd,
      diff = 0,
      N = N / 2,
      reps = reps,
      type = type,
      alpha = 0.05,
      Blockmean = Blockmean,
      BlockStdAdj = 0,
      StdAdj = StdAdj,
      seed = seed + 50,
      returnData = TRUE
    )

    CliffdMdMRE <- median(abs(Out1$Cliffd))


    PHatMdMRE <- median(abs(Out1$PHat - 0.5))

    StdESMdMRE <- median(abs(Out1$StdES))

    NullESAccuracyTable <-
      tibble::tibble(dplyr::bind_rows(
        NullESAccuracyTable,
        dplyr::bind_cols(
          Design = Design,
          Obs = as.character(2 * N),
          CliffdAbsError = CliffdMdMRE,
          PHatAbsError = PHatMdMRE,
          StdESdAbsError = StdESMdMRE,
          varCliffd = var(Out1$Cliffd),
          varPHat = var(Out1$PHat),
          varStdES = var(Out1$StdES),
          ObsCliffd = mean(Out1$Cliffd),
          ObsPHat = mean(Out1$PHat),
          ObsStdES = mean(Out1$StdES),
          CliffdType1ER = mean(Out1$CliffdSig),
          PHatType1ER = mean(Out1$PHatSig),
          StdESType1ER = mean(Out1$ESSig)
        )
      ))

    return(NullESAccuracyTable)
  }







#' @title calculateSmallSampleSizeAdjustment
#' @description Function calculates the small sample size adjustment for standardized mean effect sizes
#' @author Barbara Kitchenham and Lech Madeyski
#' @export calculateSmallSampleSizeAdjustment
#' @param df A vector of degrees of freedom
#' @param exact Default value=TRUE, if exact==TRUE the function returns the exact value of the adjustment(s) which is suitable for small values of df, if exact==FALSE the function returns the approximate version of the adjustment(s). See Hedges and Olkin 'Statistical methods for Meta-Analysis' Academic Press 1985.
#' @return small sample size adjustment value
#' @examples
#' df <- 2
#' a <- calculateSmallSampleSizeAdjustment(df)
#' # > a
#' # [1] 0.5641896
#'
#' df <- c(5, 10, 17)
#' adjexact <- calculateSmallSampleSizeAdjustment(df)
#' # > adjexact
#' # [1] 0.8407487 0.9227456 0.9551115
#' # Hedges and Olkin values 0.8408, 0.9228,0.9551
#' adjapprox <- calculateSmallSampleSizeAdjustment(df, FALSE)
#' # > adjapprox
#' # [1] 0.8421053 0.9230769 0.9552239
#' # Another example:
#' df <- c(10, 25, 50)
#' calculateSmallSampleSizeAdjustment(df, exact = TRUE)
#' # [1] 0.9227456 0.9696456 0.9849119
#' calculateSmallSampleSizeAdjustment(df, exact = FALSE)
#' # [1] 0.9230769 0.9696970 0.9849246
calculateSmallSampleSizeAdjustment <- function(df, exact = TRUE) {
  exactvec <- c(rep(exact, length(df)))
  # If exact is TRUE but the df is too large gamma cannot be calculated and the approximate value is used
  c <-
    ifelse(exactvec &
             df < 340,
           sqrt(2 / df) * gamma(df / 2) / gamma((df - 1) / 2),
           (1 - 3 / (4 * df - 1)))
  return(c)
}


#' @title varStandardizedEffectSize
#' @description Function calculates the exact variance of a standardized effect size based on the relationship between t and the standardized effect size, see Morris and DeShon, Combining Effect Size Estimates in Meta-Analysis With Repeated Measures and Independent-Groups Designs, Psychological Methods, 7 (1), pp 105-125.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export varStandardizedEffectSize
#' @param d An unadjusted standardized effect size
#' @param A The squared constant linking t and d i.e. t*sqrt(A)=d
#' @param f The degrees of freedom of the t value
#' @param returnVarg if set to TRUE return the variance of the small sample size adjusted standardized effect size (g), otherwise returns var(d) where d is the input parameter
#' @return if returnVarg if set to TRUE, return var(g) otherwise var(d)
#' @examples
#' d <- 0.5
#' varStandardizedEffectSize(d, 2 / 20, 38, returnVarg = FALSE)
#' # [1]  0.1047567
#' varStandardizedEffectSize(d, 2 / 20, 38, returnVarg = TRUE)
#' # [1] 0.1090516
varStandardizedEffectSize <- function(d, A, f, returnVarg = TRUE) {
  c <- reproducer::calculateSmallSampleSizeAdjustment(f)
  g <-
    d * c # g is a better estimate of the population standardized effect size delta than d
  var <- (f / (f - 2)) * (A + g ^ 2) - d ^ 2
  if (returnVarg) {
    var <- c ^ 2 * var
  } # best estimate of the variance of g
  return(var)
}




#' @title ExtractMAStatistics
#' @description This function extracts summary statistics from meta-analysis results obtained from the rma function of the metafor R package. If required the function transform back to standardized mean difference (effect size type 'd' i.e. Hg) or point biserial correlations (effect size type 'r').
#' Warning: the `ExtractMAStatistics` function works with `metafor` version 2.0-0, but changes to metafor's method of providing access to its individual results may introduce errors into the function.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export ExtractMAStatistics
#' @param maresults is the output from the rma function.
#' @param Nc is the number of participants in the control condition group.
#' @param Nt is the number of participants in the treatment condition group.
#' @param Transform is a boolean value indicating whether the outcome values need to be transformed back to standardized mean difference ('d' i.e. Hg or d) or point biserial correlations ('r'). It is defaulted to TRUE. If this parameter is set to FALSE, no transformation will be applied.
#' @param type this indicates the type of transformation required - it defaults to 'd' which requests transformation from Zr to Hg, using 'r' requests transformation from Zr to r.
#' @param sig indicates the number of significant digits requested in the output, the default is 4; it rounds the values of mean, pvalue, upper and lower bound to the specified number of significant digits.
#' @param returnse if set to TRUE returns the standard error of the effect size (default: returnse=FALSE)
#' @return data frame incl. summary statistics from meta-analysis results: overall mean value for the effect sizes, the p-value of the mean, the upper and lower confidence interval bounds (UB and LB), QE which is the heterogeneity test statistic and QEp which the the p-value of the heterogeneity statistic
#' @examples
#' ExpData <- reproducer::KitchenhamMadeyskiBrereton.ExpData
#' # Extract the experiment basic statics
#' S1data <- subset(ExpData, ExpData == "S1")
#' # Use the descriptive data to construct effect size
#' S1EffectSizes <- reproducer::PrepareForMetaAnalysisGtoR(
#'   S1data$Mc, S1data$Mt, S1data$SDc, S1data$SDt, S1data$Nc, S1data$Nt
#' )
#' # Do a random effect meta-analysis of the transformed r_pbs effect size
#' S1MA <- metafor::rma(S1EffectSizes$zr, S1EffectSizes$vi)
#' # Extract summary statistics from meta-analysis results and transform back to Hg scale
#' ExtractMAStatistics(S1MA, sum(S1data$Nc), sum(S1data$Nt), TRUE, "d", 4)
#' #     mean   pvalue    UB     LB QE  QEp
#' # 1 0.6658 0.002069 1.122 0.2384  4 0.41
#' ExtractMAStatistics(S1MA, sum(S1data$Nc), sum(S1data$Nt), FALSE, "d", 4)
#' # A tibble: 1 x 6
#' #   mean  pvalue    UB    LB    QE   QEp
#' #  <dbl>   <dbl> <dbl> <dbl> <dbl> <dbl>
#' # 1 0.327 0.00207 0.535 0.119     4  0.41
#'
ExtractMAStatistics <-
  function(maresults,
           Nc,
           Nt,
           Transform = TRUE,
           type = "d",
           sig = 4,
           returnse = FALSE) {
    pvalue <- as.numeric(maresults$pval)

    se <- as.numeric(maresults$se)

    QE <- as.numeric(maresults$QE)
    QEp <- as.numeric(maresults$QEp)
    mean <- as.numeric(maresults$beta)
    UB <- as.numeric(maresults$ci.ub)
    LB <- as.numeric(maresults$ci.lb)

    if (Transform & type == "d") {
      mean <- reproducer::transformZrtoHg(mean, Nc, Nt)
      se <- reproducer::transformZrtoHg(se, Nc, Nt)
      UB <- reproducer::transformZrtoHg(UB, Nc, Nt)
      LB <- reproducer::transformZrtoHg(LB, Nc, Nt)
    }
    if (Transform & type == "r") {
      mean <- reproducer::transformZrtoR(mean)
      se <- reproducer::transformZrtoR(se)

      UB <- reproducer::transformZrtoR(UB)
      LB <- reproducer::transformZrtoR(LB)
    }
    mean <- signif(mean, sig)
    pvalue <- signif(pvalue, sig)
    se <- signif(se, sig)

    UB <- signif(UB, sig)
    LB <- signif(LB, sig)
    QE <- signif(QE, 2)
    QEp <- signif(QEp, 2)
    if (returnse) {
      metaanalysisresults <-
        tibble::tibble(mean, pvalue, se, UB, LB, QE, QEp)
    } else {
      metaanalysisresults <- tibble::tibble(mean, pvalue, UB, LB, QE, QEp)
    }
    return(metaanalysisresults)
  }

#' @title NP2GMetaAnalysisSimulation
#' @description This function simulates data from a family of experiments. The parameter Exp determines the number of experiments in the family. The function simulates data from one of four distributions and uses the data to construct two of groups of equal size (GroupSize). The distribution for one  of the groups corresponds to the control and is based on the given mean and spread, the distribution for the other group corresponds to the treatment group and  is based on the mean+diff and the spread plus any variance adjustment requested (determined by the parameter StdAdj). The data from each experiment is analysed separately to estimate three non-parametric effect sizes: the Cliff's d and the probability of superiority referred to as phat and their variances. Parametric effect sizes Cohen's d (also known as the standarized means difference, SMD) and the small sample size adjusted standardized mean difference g are also calculated together with their variances. The effect sizes are then meta-analysed using various methods: the simple average of the effect size and the variance weighted averages (using the exact and approximate normal variance and the weighted and unweighted standardized mean difference). The function uses the metafor package for formal meta-analysis, and the specific method of formal meta-analysis used is determined by the MAMethod. All tests of significance are done at the 0.05 level. If the parameter returnES is TRUE, the function returns the effect sizes for each experiment in the family, otherwise it returns the meta-analysis results.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export NP2GMetaAnalysisSimulation
#' @param mean the value used for the mean of control group in the simulated data. It can be any real number including zero.
#' @param sd the value used for the spread of the control group and the spread of the treatment group in the simulated data. The value must be a real value greater than 0.
#' @param diff mean+diff is the value used for the mean of the treatment group. It can be zero.
#' @param GroupSize is the size of each of the 2 groups comprising one experiment. Groupsize should be an integer of 4 or more
#' @param Exp is the number of experiments being simulated. Exp should be an integer of 2 or more. It defaults to 5.
#' @param type specifies the distribution being simulated. The permitted values are "n" for the normal distribution,  "l" for the lognormal distribution, "g" for the gamma distribution and "lap" for the Laplace distribution. The parameter defaults to "n".
#' @param StdAdj specifies a level used to adjust the treatment variance. It allows heterogeneity to be modelled. It defaults to zero meaning no variance heterogeneity is introduced.
#' @param seed specifies the seed to be used to initiate the simulation, so the simulation is repeatable. It defauls to 123.
#' @param StdExp defines whether any additional heterogeneity is introduced between families. The value (set to 0 or 0.5 for our simulations) is used when we generate a deviation to be added to the control mean (control rate for gamma data) for each family. The deviation is generated from a Normal distribution with mean 0 and standard deviation=0.5. If StdExp=0 we do not add any deviations to the mean.
#' @param alpha the Type 1 error rate level use for statistical tests.
#' @param MAMethod the meta-analysis method needed for the call to the metafor package rma algorithm
#' @param returnES Determines the format of the output. It defaults to FALSE which causes the function to output the meta-analysis results for the family of experiments. If set to TRUE it returns the effect sizes for each experiment.
#' @param AlwaysTwoSidedTests If FALSE the function performs one-sided tests if diff!=0, and two-sided tests if diff=0. If set to TRUE the function alsways does two-sided tests.
#' @return Depending on the value of the returnES parameter, the function either returns the effect sizes for each experiment or the aggregated results for the family
#' @examples
#' as.data.frame(NP2GMetaAnalysisSimulation(mean=0,sd=1,diff=0.5,GroupSize=10,
#'   Exp=5,type="n",StdAdj=0,alpha=0.05,seed=457,StdExp=1,MAMethod="PM",
#'   returnES=FALSE))
#' #  NumExp GroupSize AveCliffd AveCliffdvar AveCliffdsig Avephat  Avephatvar Avephatsig AveMDStd..
#' #      5        10     0.252   0.01499003         TRUE   0.626 0.003645333       TRUE 0.4883188..
#' #  AveMDStdsig MAphat   MAphatvar MAphatsig MACliffd MACliffdvar MACliffdsig StdMDAdjUnweighted..
#' #1        TRUE 0.6288 0.003620188      TRUE   0.2575  0.01490134        TRUE          0.4748065..
#' #  StdMDAdjUnweightedvar StdMDAdjUnweightedsig StdMDUnweighted StdMDUnweightedvar StdMDUnweight..
#' #1            0.04065614                  TRUE       0.4980148         0.04157691            TRUE
#' #  HedgesMA.Weighted HedgesMA.Weightedvar HedgesMA.Weightedsig StdMDAdjMAexact StdMDAdjMAexactvar
#' #1         0.4755316           0.04307274                 TRUE       0.4725834         0.04315211
#' #  StdMDAdjMAexactsig StdMDAdjMAapprox StdMDAdjMAapproxvar StdMDAdjMAapproxsig StdMDMAapprox St..
#' #1               TRUE           0.4716          0.03762363                TRUE     0.4955783 ..
#' #  StdMDMAapproxsig StdMDMAexact StdMDMAexactvar StdMDMAexactsig
#' #1             TRUE    0.4966121      0.04756193            TRUE

#' as.data.frame(NP2GMetaAnalysisSimulation(mean=0,sd=1,diff=0.5,GroupSize=10,Exp=5,type="n",
#'   StdAdj=0,alpha=0.05,seed=457,StdExp=1,MAMethod="PM",returnES=TRUE))
#' #    MeanExp   VarExp     StdMD       df      tval t.sig Cliffd  Cliffdvar Cliffd.sig PHat PHat..
#' #1 0.5641594 1.437447 0.4705502 17.77980 1.0521822 FALSE   0.26 0.08149818      FALSE 0.63 0.02..
#' #2 0.6400936 1.081352 0.6155452 17.23411 1.3764009 FALSE   0.36 0.06527192      FALSE 0.68 0.01..
#' #3 0.8199650 1.698610 0.6291418 15.42141 1.4068038 FALSE   0.28 0.07362909      FALSE 0.64 0.01..
#' #4 0.2970819 1.709441 0.2272214 13.87833 0.5080824 FALSE   0.04 0.07936485      FALSE 0.52 0.01..
#' #5 0.5688567 1.079082 0.5476154 16.79899 1.2245053 FALSE   0.32 0.07498667      FALSE 0.66 0.01..
#' #  Phat.sig  StdMDAdj StdMDAdjvar.exact StdMDAdjvar.approx StdMDvar.exact StdMDvar.approx
#' #1    FALSE 0.4503698         0.2129598          0.1884384      0.2324722       0.2057040
#' #2    FALSE 0.5882961         0.2182075          0.1918563      0.2388898       0.2100409
#' #3    FALSE 0.5979539         0.2211428          0.1911344      0.2448130       0.2115926
#' #4    FALSE 0.2146782         0.2105671          0.1800107      0.2358918       0.2016604
#' #5    FALSE 0.5227345         0.2162500          0.1896495      0.2373259       0.2081330

#' as.data.frame(NP2GMetaAnalysisSimulation(mean=0,sd=1,diff=0.724,GroupSize=10,Exp=5,type="l",
#'   StdAdj=0,alpha=0.05,seed=123,StdExp=1,MAMethod="PM",returnES=FALSE))
#' #  NumExp GroupSize AveCliffd AveCliffdvar AveCliffdsig Avephat  Avephatvar Avephatsig  AveMDSt..
#' #1      5        10     0.344   0.01288023         TRUE   0.672 0.003118222       TRUE 0.483665..
#' #  AveMDStdsig MAphat   MAphatvar MAphatsig MACliffd MACliffdvar MACliffdsig StdMDAdjUnweighted
#' #1        TRUE 0.7014 0.004229764      TRUE    0.403  0.01690867        TRUE          0.5722448
#' #  StdMDAdjUnweightedvar StdMDAdjUnweightedsig StdMDUnweighted StdMDUnweightedvar StdMDUnweight..
#' #1            0.04146189                  TRUE       0.6046947         0.04260837            TRUE
#' #  HedgesMA.Weighted HedgesMA.Weightedvar HedgesMA.Weightedsig StdMDAdjMAexact StdMDAdjMAexactvar
#' #1         0.5742311           0.04453436                 TRUE       0.5405307          0.0450343
#' #  StdMDAdjMAexactsig StdMDAdjMAapprox StdMDAdjMAapproxvar StdMDAdjMAapproxsig StdMDMAapprox S..
#' #1               TRUE           0.5411          0.03819079                TRUE     0.5737401 0...
#' #  StdMDMAapproxsig StdMDMAexact StdMDMAexactvar StdMDMAexactsig
#' #1             TRUE    0.5727409      0.05042801            TRUE

NP2GMetaAnalysisSimulation = function(mean,
                                      sd,
                                      diff,
                                      GroupSize,
                                      Exp = 5,
                                      type = "n",
                                      StdAdj = 0,
                                      alpha = 0.05,
                                      seed = 123,
                                      StdExp = 0,
                                      MAMethod,
                                      returnES = FALSE,
                                      AlwaysTwoSidedTests = FALSE) {
  # 4-07-2022 Allow for one-sided tests
  base::set.seed(seed)
  N = GroupSize

  if (AlwaysTwoSidedTests | (diff == 0))
  {
    alternative = "two.sided"
    PositiveMD = FALSE
  }
  else
  {
    # Find direction of MD

    if (type == "g")
    {
      # A negative value for diff decreases the rate & the MD is positive
      Positive.MD = (diff < 0)
    }
    else {
      Positive.MD = (diff > 0)
    }
    if (Positive.MD) {
      alternative = "greater"
    }
    else {
      alternative  = "less"
    }
  }

  Family.Ktau = rep(NA, Exp) # This will hold the Kendall's tau for each experiment
  Family.Ktauvar = rep(NA, Exp) # This will hold the estimated consistent variance of tau for each experiment
  Family.PBScor = rep(NA, Exp) # This holds the Pearson point bi-serial correlation
  Family.PBScorvar = rep(NA, Exp) # This holds the variance of normal transformation of the Pearson correlation

  Family.Var = rep(NA, Exp) # This holds the pooled variance for each experiment
  Family.MD = rep(NA, Exp) # This holds the mean difference for each experiment

  Family.StdMD = rep(NA, Exp) # This holds the standardized mean difference for each experiment
  Cliffd = rep(NA, Exp) #This holds Cliff's d
  Cliffdvar = rep(NA, Exp) # This holds the variance of Clff's d
  Cliffd.sig = rep(NA, Exp) #This holds significance of Cliff's d based on the Cliff d confidence interval

  PHat = rep(NA, Exp) #This holds the probability of superiority phat
  PHatvar = rep(NA, Exp) #This holds the variance of phat
  PHatdf = rep(NA, Exp) # This holds the degrees of freedom of the t-test for phat
  PHat.sig = rep(NA, Exp) # This holds the significance of the t-test for phat

  Family.StdMDAdj = rep(NA, Exp) # This holds the small sample size adjusted standardized mean difference
  Family.StdMDAdjvar.exact = rep(NA, Exp) # This holds the exact variance of the small sample size adjusted standardized mean difference
  Family.StdMDAdjvar.approx = rep(NA, Exp) # This holds the approximate variance of the small sample size adjusted standardized mean difference.
  df = rep(NA, Exp) # This holds the degrees of freedom of the t-test of the mean difference
  c = rep(NA, Exp) # This holds the small sample size adjustment factor for the standardized mean difference
  tval = rep(NA, Exp) # This holds the t-test value for each experiment
  t.sig = rep(NA, Exp)  # This holds the significance of the t-test for each experiment
  Family.StdMDvar.exact = rep(NA, Exp) # This holds the exact variance of the standardized mean difference for each experiment
  Family.StdMDvar.approx = rep(NA, Exp) # This holds the large sample size approximate variance of the standardized mean difference for each experiment


  # For  the point-biserial tau we need a dummy variable that takes the value zero for control group data points and 1 for treatment group data points
  dummy = c(rep(0, N), rep(1, N))


  for (i in 1:Exp)
  {
    DataError = TRUE
    Numerrs = 0

    while (DataError) {
      if (StdExp == 0)
        ExpAdj = 0
      else {
        ExpAdj = rnorm(1, 0, StdExp)

        if ((type == "g") & ((mean + ExpAdj + diff) <= 0))
        {
          while ((mean + ExpAdj + diff) <= 0)
          {
            ExpAdj = rnorm(1, 0, StdExp)
          }
        }
      }

      DataSet = simulate2GExperimentData(
        mean = mean,
        sd = sd,
        diff = diff,
        GroupSize = GroupSize,
        type = type,
        ExpAdj = ExpAdj,
        StdAdj = StdAdj,
        BlockEffect = 0,
        BlockStdAdj = 0
      )

      NumNAsBD = sum(as.numeric(is.na(DataSet$BaselineData)))
      NumNAsAD = sum(as.numeric(is.na(DataSet$AlternativeData)))

      if (NumNAsBD == 0 & NumNAsAD == 0)  {
        DataError = FALSE
      }

      else {
        Numerrs = Numerrs + 1
        if (Numerrs > 50) {
          stop()
        }
        seed = seed + 1
        set.seed(seed)
      }

    }


    # Dataset holds the simulated data set with each observation assigned to a group in an experiment

    # The data in each experiment is analysed to obtain the value of Kendall's tau, its consistent variance and its hypothesis testing variance

    xy = c(DataSet$BaselineData, DataSet$AlternativeData)


    expdata = base::data.frame(xy, dummy)

    ktau = calculateKendalltaupb(expdata$xy, expdata$dummy)
    Family.Ktau[i] = ktau$cor
    Family.Ktauvar[i] = ktau$consistentvar

    #Analyse the generated data for each member of the family using the cid function to obtain Cliff's d and its variance
    Cliff = Cliffd.test(
      DataSet$AlternativeData,
      DataSet$BaselineData,
      alpha = alpha,
      alternative = alternative,
      sigfig = -1
    )
    Cliffd[i] = Cliff$d
    Cliffdvar[i] = Cliff$sqse.d
    Cliffd.sig[i] = Cliff$d.sig


    #Analyse the generated data for each member of the family using the bmp function to obtain phat and its variance

    PHat.res = PHat.test(
      DataSet$BaselineData,
      DataSet$AlternativeData,
      alpha = alpha,
      alternative = alternative,
      sigfig = -1
    )
    PHat[i] = PHat.res$phat
    PHatvar[i] = PHat.res$sqse.phat
    PHatdf[i] = PHat.res$phat.df
    PHat.sig[i] = PHat.res$phat.sig

    # Prepare to do a standard analysis for comparison
    # We can do both a standardized effect size or correlation. For the paper we use the standardized mean difference effect size

    Family.Var[i] = (stats::var(DataSet$BaselineData) + stats::var(DataSet$AlternativeData)) / 2
    Family.MD[i] = base::mean(DataSet$AlternativeData) - base::mean(DataSet$BaselineData)
    Family.StdMD[i] = (Family.MD[i]) / sqrt(Family.Var[i])


    tempttest = t.test(DataSet$AlternativeData,
                       DataSet$BaselineData,
                       alternative = alternative)
    df[i] = as.numeric(tempttest$parameter)
    tval[i] = tempttest$statistic

    t.sig[i] = tempttest$p.value < 0.05

    c[i] = reproducer::calculateSmallSampleSizeAdjustment(df[i])
    Family.StdMDAdj[i] = Family.StdMD[i] * c[i]

    # N in each group

    Family.StdMDvar.approx[i] = 2  / N + Family.StdMDAdj[i] ^ 2 / (2 * df[i])

    Family.StdMDAdjvar.approx[i] = c[i] ^ 2 * Family.StdMDvar.approx[i]

    Family.StdMDAdjvar.exact[i] = varStandardizedEffectSize(Family.StdMD[i], 2 / N, df[i], returnVarg =
                                                              TRUE)


    Family.StdMDvar.exact[i] =  varStandardizedEffectSize(Family.StdMD[i], 2 / N, df[i], returnVarg =
                                                            FALSE)
    Pearsonr = Family.StdMDAdj[i] / sqrt(Family.StdMDAdj[i] ^ 2 + 4)
    TempP = as.numeric(Pearsonr)
    TempP = reproducer::transformRtoZr(TempP)
    Family.PBScor[i] = TempP
    Family.PBScorvar[i] = 1 / (2 * N - 3)

  }

  if (!returnES)
  {
    NumExp = Exp

    # First aggregate as if the family is a planned distributed experiment, i.e. the "experiment" is a blocking factor

    # Aggregate tau_pb values


    AveKtau = base::mean(Family.Ktau)
    AveKtauctvar = sum(Family.Ktauvar) / NumExp ^ 2

    #Obtain confidence intervals using t-distribtion with 2*N-3*NumExp degrees of freedom

    AveKtausig = (
      calcEffectSizeConfidenceIntervals(
        effectsize = AveKtau,
        effectsize.variance = AveKtauctvar,
        #effectsize.df = (2 * N - 3) * (NumExp - 1), ERROR
        effectsize.df = (2 * N - 3) * NumExp,
        alpha = 0.05,
        alternative = alternative,
        UpperValue = 2 * N / (2 * (2 * N - 1)),
        LowerValue = -2 * N / (2 * (2 * N - 1))
      )
    )$ES.sig



    #Aggregate the Cliff's d values

    AveCliffd = base::mean(Cliffd)
    AveCliffdvar = sum(Cliffdvar) / NumExp ^ 2

    d.df = (2 * N - 2) * NumExp

    AveCliffdsig = (
      calcCliffdConfidenceIntervals(
        d.value = AveCliffd,
        d.variance = AveCliffdvar,
        # d.df = (2 * N - 2) * (NumExp - 1), ERROR
        d.df = d.df,
        alpha = alpha,
        alternative = alternative
      )
    )$d.sig


    #Aggregate the phat values

    Avephat = base::mean(PHat)
    Avephatadj = Avephat - 0.5 # The null effect for phat is 0.5
    Avephatvar = sum(PHatvar) / NumExp ^ 2
    # pdf = sum(PHatdf - NumExp + 1)  ERROR
    pdf = sum(PHatdf)

    Avephatsig = (
      calcPHatConfidenceIntervals(
        phat = Avephat,
        phat.variance = Avephatvar,
        phat.df = pdf,
        alpha = alpha,
        alternative = alternative
      )
    )$phat.sig



    # tau_pb meta-analysis

    meth = MAMethod

    KtauMA.res = metafor::rma(Family.Ktau, Family.Ktauvar, method = meth)

    # Extract the data from the analysis

    KtauMAResults = ExtractMAStatistics(KtauMA.res,
                                        N,
                                        N,
                                        Transform = FALSE,
                                        returnse = TRUE)
    KtauMA = KtauMAResults$mean

    KtauMAvar = KtauMAResults$se ^ 2

    KtauMAsig = (
      calcEffectSizeConfidenceIntervals(
        effectsize = KtauMA,
        effectsize.variance = KtauMAvar,
        # effectsize.df = (2 * N - 3) * (NumExp - 1), ERROR
        effectsize.df = (2 * N - 3) * NumExp,
        alpha = 0.05,
        alternative = alternative,
        UpperValue = 2 * N / (2 * (2 * N - 1)),
        LowerValue = -2 * N / (2 * (2 * N - 1))
      )
    )$ES.sig

    KtauMAQE = KtauMAResults$QE
    KtauMAQEp = KtauMAResults$QEp
    KtauMAHetsig = KtauMAQEp < 0.05


    # Cliff's d meta-analysis
    MAd.res = metafor::rma(Cliffd, Cliffdvar, method = meth)
    MAd.Results = ExtractMAStatistics(MAd.res, N, N, Transform = FALSE)
    Mean.d = MAd.Results$mean
    MAvar.d = as.numeric(MAd.res$se) ^ 2

    d.sig = (
      calcCliffdConfidenceIntervals(
        d.value = Mean.d,
        d.variance = MAvar.d,
        #d.df = (2 * N - 2) * (NumExp - 1), ERROR
        d.df = d.df,
        alpha = alpha,
        alternative = alternative
      )
    )$d.sig


    # Probability of Superiority phat meta-analysis
    PHatAdj = PHat - 0.5
    MAphat.res = metafor::rma(PHatAdj, PHatvar, method = meth) # meta-analyse phat-0.5
    MAphat.Results = ExtractMAStatistics(MAphat.res, N, N, Transform =
                                           FALSE)
    Mean.phat = MAphat.Results$mean + 0.5 # Add back the 0.5 value
    MAvar.phat = as.numeric(MAphat.res$se) ^ 2
    phat.sig = (
      calcPHatConfidenceIntervals(
        phat = Mean.phat,
        phat.variance = MAvar.phat,
        phat.df = pdf,
        alpha = alpha,
        alternative = alternative
      )
    )$phat.sig

    # Do a Pearson correlation analysis as a comparison with tau_pb

    MAPBScor.res = metafor::rma(Family.PBScor, Family.PBScorvar, method = meth)

    MAPBScorResults = ExtractMAStatistics(
      MAPBScor.res,
      N,
      N,
      Transform = TRUE,
      type = "r",
      returnse = TRUE
    )
    PBScorMA = MAPBScorResults$mean
    PBScorMAvar = MAPBScorResults$se ^ 2
    PBScorMAsig = (
      calcEffectSizeConfidenceIntervals(
        effectsize = PBScorMA,
        effectsize.variance = PBScorMAvar,
        #effectsize.df = (2 * N - 3) * (NumExp - 1), ERROR
        effectsize.df = (2 * N - 3) * NumExp,
        alpha = 0.05,
        alternative = alternative,
        UpperValue = 2 * N / (2 * N - 1),
        LowerValue = -2 * N / (2 * N - 1)
      )
    )$ES.sig
    PBScorMAhetsig = MAPBScorResults$QEp <= 0.05


    # Do an analysis of the small sample size standardized effect size as a comparison with Cliff's d

    MAgres.approx = metafor::rma(Family.StdMDAdj, Family.StdMDAdjvar.approx, method = meth)
    MAgresults.approx = ExtractMAStatistics(MAgres.approx, N, N, Transform =
                                              FALSE)
    StdMDAdjMAapprox =	MAgresults.approx$mean
    StdMDAdjMAapproxvar = as.numeric(MAgres.approx$se) ^ 2
    StdMDAdjMAapproxsig = (
      calcEffectSizeConfidenceIntervals(
        effectsize = StdMDAdjMAapprox,
        effectsize.variance = StdMDAdjMAapproxvar,
        effectsize.df = 0,
        alpha = 0.05,
        alternative = alternative
      )
    )$ES.sig


    # Note when N is the same for each study, the variance for the study based on the average of d is the same for each study. The variance of g is affected by the degrees of freedom which based on a Welch style analysis may be less than a simple ANOVA

    # Cohen d unweighted analysis


    #Family.df = sum(df) - NumExp + 1 ERROR
    Family.df = sum(df)

    AveMD = mean(Family.MD)
    AveVar = mean(Family.Var)

    AveMDVar = 2 * AveVar / (NumExp * GroupSize) # Estimate of variance for each MD

    Ave.Family.StdMD = AveMD / sqrt(AveVar)
    Family.tvalue = AveMD / sqrt(AveMDVar)

    Family.J = calculateSmallSampleSizeAdjustment(Family.df, exact = TRUE)
    varFamily.StdMD = 2 / (GroupSize * NumExp) + Family.J ^ 2 * Ave.Family.StdMD ^
      2 / (2 * Family.df)


    FamilyMDsig = (
      calcEffectSizeConfidenceIntervals(
        effectsize = Ave.Family.StdMD,
        effectsize.variance = varFamily.StdMD,
        effectsize.df = Family.df,
        alpha = 0.05,
        alternative = alternative
      )
    )$ES.sig

    SmallSampleAnalysis = metaanalyseSmallSampleSizeExperiments(Family.StdMD, df, 2 /
                                                                  N)


    StdMDAdjUnweighted = SmallSampleAnalysis$UnweightedMean


    StdMDAdjUnweightedvar = Family.J ^ 2 * (2 / (NumExp * N) + StdMDAdjUnweighted ^
                                              2 / (2 * Family.df))

    StdMDAdjUnweightedsig = (
      calcEffectSizeConfidenceIntervals(
        effectsize = StdMDAdjUnweighted,
        effectsize.variance = StdMDAdjUnweightedvar,
        effectsize.df = Family.df,
        alpha = 0.05,
        alternative = alternative
      )
    )$ES.sig

    StdMDUnweighted = mean(Family.StdMD)

    StdMDUnweightedvar = 2 / (NumExp * N) + Family.J ^ 2 * StdMDUnweighted ^
      2 / (2 * Family.df)
    StdMDUnweightedsig = (
      calcEffectSizeConfidenceIntervals(
        effectsize = StdMDUnweighted,
        effectsize.variance = StdMDUnweightedvar,
        effectsize.df = Family.df,
        alpha = 0.05,
        alternative = alternative
      )
    )$ES.sig

    HedgesMA.Weighted = SmallSampleAnalysis$WeightedMean
    HedgesMA.Weightedvar = SmallSampleAnalysis$VarWeightedMean

    HedgesMA.Weightedsig = (
      calcEffectSizeConfidenceIntervals(
        effectsize = HedgesMA.Weighted,
        effectsize.variance = HedgesMA.Weightedvar,
        effectsize.df = Family.df,
        alpha = 0.05,
        alternative = alternative
      )
    )$ES.sig


    MA.g.exact = metafor::rma(Family.StdMDAdj, Family.StdMDAdjvar.exact, method = meth)
    StdMDAdjMAexact = as.numeric(MA.g.exact$beta)
    StdMDAdjMAexactvar = as.numeric(MA.g.exact$se) ^ 2
    # Although exact variance used, metafor analysis produces CIs based on normal distribution, so for consistency with metafor we have also used the normal CIs

    StdMDAdjMAexactsig =  (
      calcEffectSizeConfidenceIntervals(
        effectsize = StdMDAdjMAexact,
        effectsize.variance = StdMDAdjMAexactvar,
        effectsize.df = 0,
        alpha = 0.05,
        alternative = alternative
      )
    )$ES.sig


    # Cohen's d metaanalysis
    StdMDUnweighted = mean(Family.StdMD)
    Cohend.res = metafor::rma(Family.StdMD, Family.StdMDvar.exact, method = meth)
    StdMDMAexact = as.numeric(Cohend.res$beta)
    StdMDMAexactvar = as.numeric(Cohend.res$se) ^ 2
    StdMDMAexactsig =  (
      calcEffectSizeConfidenceIntervals(
        effectsize = StdMDMAexact,
        effectsize.variance = StdMDMAexactvar,
        effectsize.df = 0,
        alpha = 0.05,
        alternative = alternative
      )
    )$ES.sig


    Cohend.approx.res = metafor::rma(Family.StdMD, Family.StdMDvar.approx, method = meth)
    StdMDMAapprox = as.numeric(Cohend.approx.res$beta)
    StdMDMAapproxvar = as.numeric(Cohend.approx.res$se) ^ 2

    StdMDMAapproxsig =  (
      calcEffectSizeConfidenceIntervals(
        effectsize = StdMDMAapprox,
        effectsize.variance = StdMDMAapproxvar,
        effectsize.df = 0,
        alpha = 0.05,
        alternative = alternative
      )
    )$ES.sig

    output = tibble::tibble(
      NumExp,
      GroupSize,
      #           AveKtau,
      #           AveKtauctvar,
      #           AveKtausig,
      AveCliffd,
      AveCliffdvar,
      AveCliffdsig,
      Avephat,
      Avephatvar,
      Avephatsig,
      AveMDStd = Ave.Family.StdMD,
      AveMDStdvar = varFamily.StdMD,
      AveMDStdsig = FamilyMDsig,
      #           KtauMA,
      #            KtauMAvar,
      #            KtauMAsig,
      #            KtauMAQE,
      #            KtauMAQEp,
      #            KtauMAHetsig,
      #            PBScorMA,
      #            PBScorMAvar,
      #            PBScorMAsig,
      #            PBScorMAhetsig,
      MAphat = Mean.phat,
      MAphatvar = MAvar.phat,
      MAphatsig = phat.sig,
      MACliffd = Mean.d,
      MACliffdvar = MAvar.d,
      MACliffdsig = d.sig,
      StdMDAdjUnweighted,
      StdMDAdjUnweightedvar,
      StdMDAdjUnweightedsig,
      StdMDUnweighted,
      StdMDUnweightedvar,
      StdMDUnweightedsig,
      HedgesMA.Weighted,
      HedgesMA.Weightedvar,
      HedgesMA.Weightedsig,
      StdMDAdjMAexact,
      StdMDAdjMAexactvar,
      StdMDAdjMAexactsig,
      StdMDAdjMAapprox,
      StdMDAdjMAapproxvar,
      StdMDAdjMAapproxsig,
      StdMDMAapprox,
      StdMDMAapproxvar,
      StdMDMAapproxsig,
      StdMDMAexact,
      StdMDMAexactvar,
      StdMDMAexactsig
    )
  }

  else  {
    output = tibble::as_tibble(
      dplyr::bind_cols(
        MeanExp = Family.MD,
        VarExp = Family.Var,
        StdMD = Family.StdMD,
        df = df,
        tval = tval,
        t.sig = t.sig,
        Cliffd = Cliffd,
        Cliffdvar = Cliffdvar,
        Cliffd.sig = Cliffd.sig,
        PHat = PHat,
        PHatvar = PHatvar,
        PHatdf = PHatdf,
        Phat.sig = PHat.sig,
        StdMDAdj = Family.StdMDAdj,
        StdMDAdjvar.exact =  Family.StdMDAdjvar.exact,
        StdMDAdjvar.approx = Family.StdMDAdjvar.approx,
        StdMDvar.exact = Family.StdMDvar.exact,
        StdMDvar.approx = Family.StdMDvar.approx

      )
    )
  }
  return(output)
}


#' @title NP4GMetaAnalysisSimulation
#' @description This function simulates data from a family of experiments, where the number of experiments in a family is defined by the parameter Exp. It simulates data from one of four distributions and uses the data to construct four of groups of equal size (GroupSize). Two groups are assigned as control groups and their distribution is based on the parameter, mean, and the parameter, spread. However, the mean and spread for the control group in Block 2 can be adjusted using the parameters BlockEffect and BlockStdAdj respectively. The other two groups are treatment groups and their distribution is based on the mean+diff and the spread parameter, but the distributions can be adjusted using the StdAdj, BlockEffect and BlockStdAdj parameters. The data from each experiment is analysed separately to estimate the non-parametric statistics P-hat, Cliff's d and their variances. In addition, the estimates of the standardized mean difference and the small sample size adjusted standardized mean difference are calculated. The effect size statistics are then meta-analysed using the method specified by the MAMethod parameter. We output both the average non-parametric effect statistics across the Exp experimet analysed as if they arose from a single large experiment and also the results of meta-analysising each non-parametric effect size. We use the standard parametric effect sizes and their meta-analysis as baselines.Tests of significance are one-sided if the mean difference is non-zero. If the mean difference is zero, two-sided tests are used. In addition, the user can force the use of two-sided tests using the parameter AlwaysTwoSidedTests. This should only be used for comparison with results reported in other simulation studies. The alpha parameter determines the significance level used in the tests.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export NP4GMetaAnalysisSimulation
#' @param mean The default value used for the group means in the simulated data. It can be any real number including zero.
#' @param sd The default value used for the spread of the control group and the spread of the treatment group in the simulated data. The value must be a real value greater than 0.
#' @param diff mean+diff is the value used for the mean of the treatment group. It can be zero.
#' @param GroupSize is the size of each of the 4 groups comprising one experiment. Groupsize should be an integer of 4 or more
#' @param Exp is the number of experiments being simulated. Exp should be an integer of 2 or more. It defaults to 5.
#' @param type specifies the distribution being simulated. The permitted values are "n" for the normal distribution,  "l" for the lognormal distribution, "g" for the gamma distribution and "lap" for the Laplace dsitribution. The parameter defaults to "n".
#' @param alpha the Type 1 error rate level use for statistical tests.
#' @param seed specifies the seed to be used to initiate the simulation, so the simulation is repeatable. It defaults to 123.
#' @param StdExp defines whether any additional heterogeneity is introduced between families. The value (set to 0 or 0.5 for our simulations) is used when we generate a deviation to be added to the control mean (control rate for gamma data) for each family. The deviation is generated from a Normal distribution with mean 0 and standard deviation=0.5. If StdExp=0 we do not add any deviations to the mean.
#' @param BlockEffect is the effect of having two different blocks
#' @param BlockStdAdj is the variance associated with the Block. If BlockStdAdj is zero it means we are treating the block effect as a fixed effect. If BlockStdAdj>0, we treat the block effect as a random effect and increase the variance of Block 2 data.
#' @param StdAdj The value used to introduce heterogeneity into the treatment groups variance if required.
#' @param MAMethod defines the method used for meta-analysis
#' @param returnES This determines the format of the output. If returnES=FALSE it returns the summary meta-analysis statistics otherwise it returns the effect sizes and their variances for each experiment in the family
#' @param AlwaysTwoSidedTests If this parameter is TRUE, the function always does two-sided tests. IF the parameter is FALSE, the function does two-sided statistical tests if the difference between treatment groups is 0, if the difference is not 0, it does one-sided tests
#' @return If returnES is FALSE, the function returns the summary meta-analysis summary statistics otherwise the function returns the effect sizes for each experiment
#' @examples

#' as.data.frame(NP4GMetaAnalysisSimulation(mean=0,sd=1,diff=0.8,GroupSize=5,Exp=5,type="n",
#' alpha=0.05,seed=457,StdAdj=0,BlockEffect=0.5,BlockStdAdj=0,StdExp=0,MAMethod="FE",returnES=TRUE))
#' #     MeanExp    VarExp       StdMD       df       tval t.sig Cliffd  Cliffdvar PHat PHatvar PH..
#' #1  1.0761565 1.3874542  0.91362108 14.42773  2.0429188  TRUE   0.52 0.05530667 0.76  0.0132 13..
#' #2  0.1012680 0.9779431  0.10240368 12.74930  0.2289816 FALSE   0.20 0.09048000 0.60  0.0224 10..
#' #3  1.2100986 0.9909894  1.21558760 11.16850  2.7181365  TRUE   0.64 0.04720000 0.82  0.0110 13..
#' #4 -0.1452027 2.3106703 -0.09552252 11.93764 -0.2135949 FALSE   0.04 0.09888000 0.52  0.0244 10..
#' #5  1.1701075 0.9623530  1.19277505 12.72802  2.6671261  TRUE   0.52 0.05048000 0.76  0.0124 15..
#' #     StdMDAdj StdMDAdjvar.exact StdMDAdjvar.approx StdMDvar.exact StdMDvar.approx
#' #1  0.86514731         0.2357247          0.2025998      0.2664156       0.2259389
#' #2  0.09623845         0.2098977          0.1769637      0.2377103       0.2003632
#' #3  1.13176658         0.2732955          0.2230773      0.3262821       0.2573441
#' #4 -0.08937076         0.2106627          0.1753619      0.2407210       0.2003345
#' #5  1.12084087         0.2623764          0.2201822      0.3050637       0.2493511

#' as.data.frame(NP4GMetaAnalysisSimulation(mean=0,sd=1,diff=0.8,GroupSize=5,Exp=5,type="n",
#'alpha=0.05,seed=457,StdAdj=0,BlockEffect=0.5,BlockStdAdj=0,StdExp=0,MAMethod="FE",returnES=FALSE))
#' #  NumExp GroupSize AveCliffd AveCliffdvar AveCliffdsig Avephat Avephatvar Avephatsig  AveMDStd..
#' #1      5         5     0.384   0.01369387         TRUE   0.692   0.003336       TRUE 0.5927084..
#' #  AveMDStdsig    MAphat  MAphatvar MAphatsig  MACliffd MACliffdvar MACliffdsig StdMDAdjUnweigh..
#' #1        TRUE 0.7253858 0.00300356      TRUE 0.4471125  0.01246219        TRUE          0.6249..
#' #  StdMDAdjUnweightedvar StdMDAdjUnweightedsig StdMDUnweighted StdMDUnweightedvar StdMDUnweight..
#' #1            0.04220968                  TRUE        0.665773         0.04366035            TRUE
#' #  HedgesMA.Weighted HedgesMA.Weightedvar HedgesMA.Weightedsig StdMDAdjMAexact StdMDAdjMAexactvar
#' #1         0.6250243           0.04574766                 TRUE       0.5709401         0.04711703
#' #  StdMDAdjMAexactsig StdMDAdjMAapprox StdMDAdjMAapproxvar StdMDAdjMAapproxsig StdMDMAapprox St..
#' #1               TRUE        0.5715637          0.03950437                TRUE     0.6090632 0...
#' #  StdMDMAapproxsig StdMDMAexact StdMDMAexactvar StdMDMAexactsig
#' #1             TRUE    0.6013198      0.05417894            TRUE

NP4GMetaAnalysisSimulation = function(mean,
                                      sd,
                                      diff,
                                      GroupSize,
                                      Exp = 5,
                                      type = "n",
                                      alpha = 0.05,
                                      seed = 123,
                                      StdAdj = 0,
                                      BlockEffect = 0,
                                      BlockStdAdj = 0,
                                      StdExp = 0,
                                      MAMethod,
                                      returnES = FALSE,
                                      AlwaysTwoSidedTests = FALSE) {
  # 14-03-2022  Parameter name BlockStdadj changed to BlockStdAdj for consisncy with other functions
  # 1-07-2022 Allow for one-sided tests
  # 3-07-2022 Restructured to improve modularization
  # 26-10-2022 Changed after finalizing the equations for var(SMD) and var(SMDAdj)

  N = GroupSize

  alpha.set = alpha

  if (AlwaysTwoSidedTests | (diff == 0)) {
    alternative = "two.sided"
  }
  else{
    # Find direction of MD

    if (type == "g")
    {
      # A negative value for diff decreases the rate & the MD is positive
      Positive.MD = (diff < 0)
    }
    else {
      Positive.MD = (diff > 0)
    }
    if (Positive.MD)
    {
      alternative = "greater"
    }
    else
    {
      alternative = "less"
    }
  }

  set.seed(seed)
  Family.Ktau = rep(NA, Exp) # Used to hold Kendall's tau for each experiment
  Family.Ktauvar = rep(NA, Exp) # Used to hold the consistent variance of Kendall's tau for each experiment
  Cliffd = rep(NA, Exp) # Used to hold Cliff's d for each experiment
  Cliffdvar = rep(NA, Exp) # Used to hold the variance Cliff's d for each experiment
  Cliffd.sig = rep(NA, Exp) #This holds significance of Cliff's d based on the Cliff d confidence interval


  PHat = rep(NA, Exp) #This holds the probability of superiority phat
  PHatvar = rep(NA, Exp) #This holds the variance of phat
  PHatdf = rep(NA, Exp) # This holds the degrees of freedom of the t-test for phat
  PHat.sig = rep(NA, Exp) # This holds the significance of the t-test for phat

  Family.StdMDAdj = rep(NA, Exp) # Used to hold the small sample size adjusted standardized effect size
  Family.StdMDAdjvar.approx = rep(NA, Exp) # Used to hold the large sample size approximate variance of the small sample size adjusted standardized effect size
  Family.StdMDAdjvar.exact = rep(NA, Exp) # Used to hold the exact variance of the small sample size adjusted standardized effect size
  Family.PBScor = rep(NA, Exp) # This holds the estimate of the point bi-serial correlation effect size
  Family.PBScorvar = rep(NA, Exp) # This holds the variance of normal transformation of the correlation effect size
  c = rep(NA, Exp) # This holds the small sample size adjustment factor for each experiment
  df = rep(NA, Exp) # This holds the degrees of freedom of the t test of the unstandardized effect size for each for each experiment

  Family.StdMD = rep(NA, Exp) # This holds standardized mean difference effect size of each experiment
  Family.StdMDvar.exact  = rep(NA, Exp)  # This holds the variance of the standardized mean difference effect size of each experiment
  Family.StdMDvar.approx = rep(NA, Exp) # This holds the large sample size approximate variance of the standardized mean difference effect size of each experiment
  tval = rep(NA, Exp) # This holds t test value for each experiment
  t.sig = rep(NA, Exp)  # This holds the significance of the t-test for each experiment

  Family.MD = rep(NA, Exp) # This holds the unstandardized effect size
  Family.Var = rep(NA, Exp) # This holds the pooled within group variance

  retry = 0
  trySimulation = TRUE
  NumRetries = 0

  while (trySimulation) {
    # Generate 4-group experiments for a family of experiments (the number in the family defined by the "Exp" parameter), with the underlying distribution defined by the "type" parameter

    for (i in 1:Exp) {
      DataError = TRUE
      Numerrs = 0

      while (DataError) {
        # 17-03-2022 Revised to stop bias induced by using only a positive
        if (StdExp == 0)
          ExpAdj = 0
        else {
          ExpAdj = rnorm(1, 0, StdExp)

          if ((type == "g") & (mean + ExpAdj + diff <= 0))
          {
            while (mean + ExpAdj + diff <= 0)
            {
              ExpAdj = rnorm(1, 0, StdExp)
            }
          }
        }

        DataSet = simulate4GExperimentData(
          mean = mean,
          sd = sd,
          diff = diff,
          GroupSize = GroupSize,
          type = type,
          ExpAdj = ExpAdj,
          StdAdj = StdAdj,
          BlockEffect = BlockEffect,
          BlockStdAdj = BlockStdAdj
        )


        NumNAsBD1 = sum(as.numeric(is.na(DataSet$BaselineData.B1)))
        NumNAsAD1 = sum(as.numeric(is.na(DataSet$AlternativeData.B1)))
        NumNAsBD2 = sum(as.numeric(is.na(DataSet$BaselineData.B2)))
        NumNAsAD2 = sum(as.numeric(is.na(DataSet$AlternativeData.B2)))


        if (NumNAsBD1 == 0 &
            NumNAsAD1 == 0 & NumNAsBD2 == 0 & NumNAsAD2 == 0)  {
          DataError = FALSE
        }
        else {
          Numerrs = Numerrs + 1
          if (Numerrs > 50) {
            stop()
          }
          seed = seed + 1
          set.seed(seed)
        }


      }


      newlist = list()
      newlist[[1]] = DataSet$BaselineData.B1
      newlist[[2]] = DataSet$AlternativeData.B1
      newlist[[3]] = DataSet$BaselineData.B2
      newlist[[4]] = DataSet$AlternativeData.B2

      # The data in each experiment is analysed to obtain the values of three non-parametric effect size: Kendall's tau, its consistent variance, Cliff's d with its consistent variance and phat with its variance plus t-test results based on the Welch method applied to a randomized blocks experiment.

      NPStats = NULL

      # Calc4GroupNPStats finds the average value of various non-parametric statistics for a randomised block experimental design
      NPStats = Calc4GroupNPStats(
        newlist[[2]],
        newlist[[1]],
        newlist[[4]],
        newlist[[3]],
        alternative = alternative,
        alpha = alpha
      )
      Family.Ktau[i] = as.numeric(NPStats$cor)
      Family.Ktauvar[i] = NPStats$ctvar
      Cliffd[i] = NPStats$d
      Cliffdvar[i] = NPStats$vard
      PHat[i] = NPStats$phat
      PHatvar[i] = NPStats$phat.var
      PHatdf[i] = NPStats$phat.df # Holds the degrees of freedom



      # Do a parametric analysis of the data
      Family.MD[i] = (
        base::mean(newlist[[2]]) + base::mean(newlist[[4]]) - base::mean(newlist[[1]]) - base::mean(newlist[[3]])
      ) /
        2

      Family.Var[i] = (
        stats::var(newlist[[2]]) + stats::var(newlist[[4]]) + stats::var(newlist[[1]]) + stats::var(newlist[[3]])
      ) /
        4

      Family.StdMD[i] = (Family.MD[i]) / sqrt(Family.Var[i])

      vec = c(-1, 1, -1, 1) / 2

      # RandomizedBlocksAnalysis does a Welch-based linear contrast analysis
      res.t = RandomizedBlocksAnalysis(newlist,
                                       con = vec,
                                       alpha = alpha,
                                       alternative = alternative)

      dataframe.ttest = base::data.frame(res.t$psihat)

      testres = base::data.frame(res.t$test)

      tval[i] = testres$test

      tpval = dataframe.ttest$p.value
      df[i] = testres$df

      t.sig[i] = as.logical(res.t$sig)


      # Find small sample-size adjusted standardized mean difference

      c[i] = reproducer::calculateSmallSampleSizeAdjustment(df[i])
      Family.StdMDAdj[i] = Family.StdMD[i] * c[i]



      Family.StdMDAdjvar.exact[i] = varStandardizedEffectSize(mean(Family.StdMDAdj[i]), 1 / N, df[i], returnVarg = TRUE)


      Family.StdMDvar.approx[i] = 1 / N + Family.StdMDAdj[i] ^ 2 / (2 * df[i])

      Family.StdMDAdjvar.approx[i] = c[i] ^ 2 * Family.StdMDvar.approx[i]

      Family.StdMDvar.exact[i] = varStandardizedEffectSize(mean(Family.StdMD[i]), 1 / N, df[i], returnVarg =
                                                             FALSE)

      # Calculate r from d

      Family.PBScor[i] = Family.StdMDAdj[i] / sqrt(Family.StdMDAdj[i] ^ 2 + 4)
      Family.PBScor[i] = reproducer::transformRtoZr(Family.PBScor[i])
      Family.PBScorvar[i] = 1 / (4 * N - 3)

    } # End of for statement. All data sets generated and effect sizes found


    if (returnES) {
      # Tabulate the basic experiment statistics for each experiment
      trySimulation = FALSE


      output = tibble::tibble(
        dplyr::bind_cols(
          MeanExp = Family.MD,
          VarExp = Family.Var,
          StdMD = Family.StdMD,
          df = df,
          tval = tval,
          t.sig = t.sig,
          Cliffd = Cliffd,
          Cliffdvar = Cliffdvar,
          PHat = PHat,
          PHatvar = PHatvar,
          PHatdf = PHatdf,
          StdMDAdj = Family.StdMDAdj,
          StdMDAdjvar.exact =  Family.StdMDAdjvar.exact,
          StdMDAdjvar.approx = Family.StdMDAdjvar.approx,
          StdMDvar.exact = Family.StdMDvar.exact,
          StdMDvar.approx = Family.StdMDvar.approx
        )
      )


    }
    else
    {
      # Try to do a meta-analysis but cater for rma iterations failing to converge


      NumExp = Exp
      # Analyse the data from the Exp experiments as a single distributed experiment
      AveKtau = base::mean(Family.Ktau)
      AveKtauctvar = sum(Family.Ktauvar) / NumExp ^ 2

      # Use t-based bounds

      AveKtausig =      (
        calcEffectSizeConfidenceIntervals(
          effectsize = AveKtau,
          effectsize.variance = AveKtauctvar,
          #     effectsize.df = (4 * N - 3) * (NumExp - 1), ERROR
          effectsize.df = (4 * N - 3) * NumExp	,
          alpha = 0.05,
          alternative = alternative,
          UpperValue = 4 * N / (2 * 4 * (N - 1)),
          LowerValue = -4 * N / (2 * 4 * (N - 1))
        )
      )$ES.sig

      AveCliffd = base::mean(Cliffd)
      AveCliffdvar = sum(Cliffdvar) / NumExp ^ 2
      Cliffd.df = 4 * (N - 1) * NumExp

      AveCliffdsig = (
        calcCliffdConfidenceIntervals(
          d.value = AveCliffd,
          d.variance = AveCliffdvar,
          # d.df = 4 * (N - 1) * (NumExp - 1), ERROR
          d.df = Cliffd.df,
          alpha = alpha,
          alternative = alternative
        )
      )$d.sig


      Avephat = base::mean(PHat)
      Avephatadj = Avephat - 0.5 # The null effect for phat is 0.5
      Avephatvar = sum(PHatvar) / NumExp ^ 2
      #phatdf = sum(PHatdf - NumExp + 1) ERROR
      phatdf = sum(PHatdf)

      Avephatsig = (
        calcPHatConfidenceIntervals(
          phat = Avephat,
          phat.variance = Avephatvar,
          phat.df = phatdf,
          alpha = alpha,
          alternative = alternative
        )
      )$phat.sig


      # Perform a standard meta-analysis for  all effect sizes but allow for pathalogical simulated data

      meth = MAMethod

      Ktau.res = NULL
      MA.dres = NULL
      MA.phat = NULL
      MA.g.exact = NULL
      MA.g.approx = NULL
      MAP.res = NULL
      Cohend.res = NULL
      Cohendapprox.res = NULL


      Ktau.res =  CatchError(metafor::rma(Family.Ktau, Family.Ktauvar, method = meth))

      if (is.character(Ktau.res)) {
        problem = "Ktau.res Problem"


      }
      MA.dres = CatchError(metafor::rma(Cliffd, Cliffdvar, meth =
                                          meth))

      if (is.character(MA.dres)) {
        problem = "MA.dres Problem"

      }

      PHatAdj = PHat - 0.5
      MA.phat =  CatchError(metafor::rma(PHatAdj, PHatvar, method = meth))
      if (is.character(MA.phat)) {
        problem = "MA.phat Problem"

      }


      MA.g.exact = CatchError(metafor::rma(Family.StdMDAdj,
                                           Family.StdMDAdjvar.exact,
                                           method = meth))
      if (is.character(MA.g.exact)) {
        problem = "MA.g.exact Problem"

      }

      MA.g.approx = CatchError(metafor::rma(Family.StdMDAdj,
                                            Family.StdMDAdjvar.approx,
                                            method = meth))
      if (is.character(MA.g.approx)) {
        problem = "MA.g.approx Problem"

      }

      MAPBScor.res = CatchError(metafor::rma(Family.PBScor, Family.PBScorvar, method = meth))
      if (is.character(MAP.res)) {
        problem = "MAP.res Problem"
      }


      Cohend.res = CatchError(metafor::rma(Family.StdMD, Family.StdMDvar.exact, method = meth))
      if (is.character(Cohend.res)) {
        problem = "Cohend.res Problem"
      }

      Cohend.approx.res = CatchError(metafor::rma(Family.StdMD, Family.StdMDvar.approx, method = meth))
      if (is.character(Cohendapprox.res)) {
        problem = "Cohendapprox.res Problem"
      }

      if (is.character(Ktau.res) |
          is.character(MA.dres) |
          is.character(MA.phat) |
          is.character(MA.g.exact) |
          is.character(MA.g.approx) |
          is.character(MAPBScor.res) |
          is.character(Cohend.res) |
          is.character(Cohend.approx.res))  		{
        # Problem with meta-analyis try another simulation
        retry = retry + 1
        NumRetries = NumRetries + 1
        trySimulation = TRUE
        if (retry == 100) {
          #Give up
          problem = paste(problem, as.character(seed))
          stop(problem)
        }

        else {
          seed = seed + 1
          set.seed(seed)

        }

      }
      else
      {
        # Meta-analysis successful. No need to try a different dataset

        trySimulation = FALSE
        # Extract the data from the meta-analysis for ktau

        KtauMA = as.numeric(Ktau.res$beta)
        KtauMAvar = as.numeric(Ktau.res$se ^ 2)
        KtauMAQE = as.numeric(Ktau.res$QE)
        KtauMAQEp = as.numeric(Ktau.res$QEp)

        KtauMAsig = (
          calcEffectSizeConfidenceIntervals(
            effectsize = KtauMA,
            effectsize.variance = KtauMAvar,
            #effectsize.df = 4 * (N - 3) * (NumExp - 1),ERROR
            effectsize.df = 4 * (N - 3) * NumExp,
            alpha = 0.05,
            alternative = alternative ,
            UpperValue = 2 * N / (4 * N - 1),
            LowerValue = -2 * N / (4 * N - 1)
          )
        )$ES.sig

        KtauMAHetsig = KtauMAQEp < 0.05

        #	Meta-analyse Cliff'd d result

        Mean.d = as.numeric(MA.dres$beta)

        MAvar.d = as.numeric(MA.dres$se) ^ 2

        d.sig = (
          calcCliffdConfidenceIntervals(
            d.value = Mean.d,
            d.variance = MAvar.d,
            # d.df = 4 * (N - 1) * (NumExp - 1), ERROR
            d.df = Cliffd.df,
            alpha = alpha,
            alternative = alternative
          )
        )$d.sig


        #	Meta-analysis of phat

        Mean.phat = as.numeric(MA.phat$beta) + 0.5 # Add the 0.5 back

        MAvar.phat = as.numeric(MA.phat$se) ^ 2

        phat.sig = (
          calcPHatConfidenceIntervals(
            phat = Mean.phat,
            phat.variance = MAvar.phat,
            phat.df = phatdf,
            alpha = alpha,
            alternative = alternative
          )
        )$phat.sig

        # Meta-analysis of StdDMAdj values using both exact and approx variance.

        # Note when N is the same for each study the variance for the study based on the average of d is the same for each study. The variance of g is affected by the degrees of freedom which based on a Welch style analysis may be less than a simple ANOVA

        StdMDAdjMAexact = as.numeric(MA.g.exact$beta)

        StdMDAdjMAexactvar = as.numeric(MA.g.exact$se) ^ 2

        StdMDAdjMAexactsig =
          (
            calcEffectSizeConfidenceIntervals(
              effectsize = StdMDAdjMAexact,
              effectsize.variance = StdMDAdjMAexactvar,
              effectsize.df = 0,
              alpha = 0.05,
              alternative = alternative
            )
          )$ES.sig



        StdMDAdjMAapprox = as.numeric(MA.g.approx$beta)

        StdMDAdjMAapproxvar = as.numeric(MA.g.approx$se) ^ 2

        StdMDAdjMAapproxsig = (
          calcEffectSizeConfidenceIntervals(
            effectsize = StdMDAdjMAapprox,
            effectsize.variance = StdMDAdjMAapproxvar,
            effectsize.df = 0,
            alpha = 0.05,
            alternative = alternative
          )
        )$ES.sig


        # Do a Pearson correlation analysis


        MAPBScorResults = ExtractMAStatistics(
          MAPBScor.res,
          N,
          N,
          Transform = TRUE,
          type = "r",
          returnse = TRUE
        )
        PBScorMA = MAPBScorResults$mean
        PBScorMAvar = MAPBScorResults$se ^ 2

        PBScorMAsig = (
          calcEffectSizeConfidenceIntervals(
            effectsize = PBScorMA,
            effectsize.variance = PBScorMAvar,
            #effectsize.df = (4 * N - 3) * (NumExp - 1), ERROR
            effectsize.df = (4 * N - 3) * NumExp,
            alpha = 0.05,
            alternative = alternative,
            UpperValue = 4 * N / (4 * N - 1),
            LowerValue = -4 * N / (4 * N - 1)
          )
        )$ES.sig

        PBScorMAhetsig = MAPBScorResults$QEp <= 0.05

        # Meta-analysis of Cohen's d

        StdMDMAexact = as.numeric(Cohend.res$beta)
        StdMDMAexactvar = as.numeric(Cohend.res$se) ^ 2
        StdMDMAexactsig =  (
          calcEffectSizeConfidenceIntervals(
            effectsize = StdMDMAexact,
            effectsize.variance = StdMDMAexactvar,
            effectsize.df = 0,
            alpha = 0.05,
            alternative = alternative
          )
        )$ES.sig

        StdMDMAapprox = as.numeric(Cohend.approx.res$beta)
        StdMDMAapproxvar = as.numeric(Cohend.approx.res$se) ^ 2

        StdMDMAapproxsig =  (
          calcEffectSizeConfidenceIntervals(
            effectsize = StdMDMAapprox,
            effectsize.variance = StdMDMAapproxvar,
            effectsize.df = 0,
            alpha = 0.05,
            alternative = alternative
          )
        )$ES.sig



        # Family.df = sum(df) - NumExp + 1 ERROR
        Family.df = sum(df)

        AveMD = mean(Family.MD)
        AveVar = mean(Family.Var)

        AveMDVar = AveVar / (NumExp * GroupSize) # Estimate of variance for each MD

        Ave.Family.StdMD = AveMD / sqrt(AveVar)

        Family.J = calculateSmallSampleSizeAdjustment(Family.df, exact = TRUE)
        varFamily.StdMD = 1 / (GroupSize * NumExp) + Family.J ^
          2 * Ave.Family.StdMD ^ 2 / (2 * Family.df)


        FamilyMDsig = (
          calcEffectSizeConfidenceIntervals(
            effectsize = Ave.Family.StdMD,
            effectsize.variance = varFamily.StdMD,
            effectsize.df = Family.df,
            alpha = 0.05,
            alternative = alternative
          )
        )$ES.sig


        SmallSampleAnalysis = metaanalyseSmallSampleSizeExperiments(Family.StdMD, df, 1 /
                                                                      N)

        StdMDAdjUnweighted = SmallSampleAnalysis$UnweightedMean

        StdMDAdjUnweightedvar = Family.J ^ 2 * (1 / (NumExp * N) +
                                                  StdMDAdjUnweighted ^ 2 / (2 * Family.df))
        StdMDAdjUnweightedsig = (
          calcEffectSizeConfidenceIntervals(
            effectsize = StdMDAdjUnweighted,
            effectsize.variance = StdMDAdjUnweightedvar,
            effectsize.df = Family.df,
            alpha = 0.05,
            alternative = alternative
          )
        )$ES.sig

        StdMDUnweighted = mean(Family.StdMD)

        StdMDUnweightedvar = 1 / (NumExp * N) + Family.J ^ 2 * StdMDUnweighted ^
          2 / (2 * Family.df)

        StdMDUnweightedsig = (
          calcEffectSizeConfidenceIntervals(
            effectsize = StdMDUnweighted,
            effectsize.variance = StdMDUnweightedvar,
            effectsize.df = Family.df,
            alpha = 0.05,
            alternative = alternative
          )
        )$ES.sig

        HedgesMA.Weighted = SmallSampleAnalysis$WeightedMean
        HedgesMA.Weightedvar = SmallSampleAnalysis$VarWeightedMean

        HedgesMA.Weightedsig = (
          calcEffectSizeConfidenceIntervals(
            effectsize = HedgesMA.Weighted,
            effectsize.variance = HedgesMA.Weightedvar,
            effectsize.df = Family.df,
            alpha = 0.05,
            alternative = alternative
          )
        )$ES.sig

        #Tabulate the results for output.

        output = tibble::tibble(
          NumExp,
          GroupSize,
          #           AveKtau,
          #           AveKtauctvar,
          #           AveKtausig,
          AveCliffd,
          AveCliffdvar,
          AveCliffdsig,
          Avephat,
          Avephatvar,
          Avephatsig,
          AveMDStd = Ave.Family.StdMD,
          AveMDStdvar = varFamily.StdMD,
          AveMDStdsig = FamilyMDsig,
          #           KtauMA,
          #            KtauMAvar,
          #            KtauMAsig,
          #            KtauMAQE,
          #            KtauMAQEp,
          #            KtauMAHetsig,
          #            PBScorMA,
          #            PBScorMAvar,
          #            PBScorMAsig,
          #            PBScorMAhetsig,
          MAphat = Mean.phat,
          MAphatvar = MAvar.phat,
          MAphatsig = phat.sig,
          MACliffd = Mean.d,
          MACliffdvar = MAvar.d,
          MACliffdsig = d.sig,
          StdMDAdjUnweighted,
          StdMDAdjUnweightedvar,
          StdMDAdjUnweightedsig,
          StdMDUnweighted,
          StdMDUnweightedvar,
          StdMDUnweightedsig,
          HedgesMA.Weighted,
          HedgesMA.Weightedvar,
          HedgesMA.Weightedsig,
          StdMDAdjMAexact,
          StdMDAdjMAexactvar,
          StdMDAdjMAexactsig,
          StdMDAdjMAapprox,
          StdMDAdjMAapproxvar,
          StdMDAdjMAapproxsig,
          StdMDMAapprox,
          StdMDMAapproxvar,
          StdMDMAapproxsig,
          StdMDMAexact,
          StdMDMAexactvar,
          StdMDMAexactsig
        )
      }

    }
  }
  return(output)

}



############################################################################################################
#' @title MetaAnalysisSimulations
#' @description This function simulates data from many families of experiments. The number of families simulated is defined by the Replications parameter. The parameter Exp determines the number of experiments in each family. The function simulates data from one of four distributions and uses the data to construct two of groups of equal size (GroupSize). The experimental design of individual experiments in each family is determined by the FourGroup parameter. If FourGroup=FALSE, the basic experimental design is a balanced two group randomized experiment, otherwise the experimental design is a balanced four group experiment corresponding to a randomized blocks experiment. The function calls either NP2GMetaAnalysisSimulation or NP2GMetaAnalysisSimulation to generate and analyse data for each individual family. The function either returns the meta-analysed data from each experiment or provides summary statistics.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export MetaAnalysisSimulations
#' @param mean the value used for the mean of control group in the simulated data. It can be any real number including zero.
#' @param sd the value used for the spread of the control group and the spread of the treatment group in the simulated data. The value must be a real value greater than 0.
#' @param diff mean+diff is the value used for the mean of the treatment group. It can be zero.
#' @param GroupSize is the size of each of the groups comprising one experiment. Groupsize should be an integer of 4 or more
#' @param Exp is the number of experiments in each family of experiments being
#' simulated. Exp should be an integer of 2 or more (default 5).
#' @param Replications The number of times the set of experiments is simulated.
#' @param type specifies the distribution being simulated. The permitted values
#' are 'n' for the normal distribution (default), 'l' for the lognormal
#' distribution, 'g' for the gamma distribution and 'lap' for the Laplace
#' distribution.
#' @param seed specifies the seed to be used to initiate the simulation, so the simulation is repeatable.
#' @param alpha The significance level used for tests and confidence intervals
#' (default 0.05)
#' @param FourGroup is a Boolean variable that determines whether the experiment
#' is a two group experiments or a 4-Group randomised block experiment. It
#' defaults to FALSE which means a two-group experiment is the default condition
#' @param StdAdj If non-zero that can be used to introduced variability into
#' the treatment spread/variance (default 0). Not appropriate for gamma data.
#' @param BlockEffect A factor used to change the mean difference between
#' blocks (default 0.5)
#' @param BlockStdAdj if non-zero this can be used to change the BlockEffect
#' from a fixed to random effect (default 0).
#' @param StdExp if non-zero it simulates a random effect between experiments
#' in the same family (default 0).
#' @param MAMethod specifies the model to be used when experimental effect sizes ar aggregated using the R metafor package.
#' @param returnES if TRUE the function outputs the summary statistics otherwise it outputs the meta-analysis results for each family (default FALSE)
#' @param AlwaysTwoSidedTests This parameter can be used to override the
#' one-sided tests used as default if the diff parameter is non-zero
#' (default FALSE). This should only be set to TRUE to check simulation
#' reported in other papers that seem to have used two-sided tests.
#' @return The parameter either returns the meta-analysis values obtained from
#' each family or the average values of the meta-analysis over all replications.
#' @examples
#' as.data.frame(
#'   MetaAnalysisSimulations(
#'     mean=0, sd=1, diff=0.5, GroupSize=10, type='n', Replications=5, Exp=5,
#'     seed=456, alpha=0.05, FourGroup=FALSE, StdAdj=0, BlockEffect=0.5,
#'     BlockStdAdj=0,StdExp=0,MAMethod='PM',returnES=FALSE))
#' #  AverageCliffd AverageCliffdvar AverageCliffdsig Averagephat Averagephatvar
#' #1        0.3336        0.0132419              0.8      0.6668    0.003214756
#' #Averagephatsig  AveMDStd AveMDStdvar AveMDStdsig MAMean.phat  MAphat.var
#' #1          0.9 0.6176206  0.04278117         0.9    0.689908 0.003888047
#' #MAphat.sig MAMean.Cliffd MACliffd.var MACliffd.sig Mean.StdMDUnweighted
#' #1      0.9       0.37984   0.01575063          0.9            0.6449963
#' #StdMDUnweighted.var StdMDUnweighted.sig Mean.StdMDAdjUnweighted
#' #1        0.04299001                 0.9               0.6145034
#' #StdMDAdjUnweighted.var StdMDAdjUnweighted.sig Mean.HedgesMA Hedges.var
#' #1           0.04192908                    0.9     0.6150575 0.04455833
#' #Hedges.sig Mean.StdMDAdjMA.exact StdMDAdjMA.exact.var StdMDAdjMA.exact.sig
#' #1      0.9             0.5834754           0.05171067                  0.8
#' #Mean.StdMDAdjMA.approx StdMDAdjMA.approx.var StdMDAdjMA.approx.sig
#' #1              0.58643            0.04749064                   0.9
#' #Mean.StdMDMA.exact StdMDMA.exact.var StdMDMA.exact.sig Mean.StdMDMA.approx
#' #1        0.6134374        0.05711235               0.8           0.6165884
#' #StdMDMA.approx.var StdMDMA.approx.sig
#' #1       0.05242339                0.9
#' #as.data.frame(
#'  # MetaAnalysisSimulations(
#'  #   mean=0, sd=1, diff=0.5, GroupSize=10, type='n', Replications=50, Exp=5,
#'  #   seed=456, alpha=0.05, FourGroup=FALSE, StdAdj=0, BlockEffect=0.5,
#'  #   BlockStdAdj=0,StdExp=0,MAMethod='PM',returnES=FALSE))
#' # AverageCliffd AverageCliffdvar AverageCliffdsig Averagephat Averagephatvar
#' #1      0.29808       0.01333744             0.74     0.64904    0.003236444
#' # Averagephatsig  AveMDStd   AveMDStdvar AveMDStdsig MAMean.phat  MAphat.var
#' #           0.78 0.5450377    0.04217901        0.78   0.6677884 0.004538661
#' #  MAphat.sig MAMean.Cliffd MACliffd.var MACliffd.sig  Mean.StdMDUnweighted
#' #1       0.72     0.3356298   0.01833956         0.72             0.5686653
#' #  StdMDUnweighted.var StdMDUnweighted.sig Mean.StdMDAdjUnweighted
#' #1          0.04237386                0.82               0.5419554
#' #StdMDAdjUnweighted.var StdMDAdjUnweighted.sig Mean.HedgesMA Hedges.var
#' #            0.04138573                   0.78     0.5420552 0.04388383
#' #  Hedges.sig Mean.StdMDAdjMA.exact StdMDAdjMA.exact.var StdMDAdjMA.exact.sig
#' #1       0.76             0.5163304           0.05874152                 0.72
#' #Mean.StdMDAdjMA.approx StdMDAdjMA.approx.var StdMDAdjMA.approx.sig
#' #1            0.5203279            0.05591752                  0.74
#' # Mean.StdMDMA.exact StdMDMA.exact.var
#' #        0.5418705        0.06468786
#' # StdMDMA.exact.sig Mean.StdMDMA.approx StdMDMA.approx.var StdMDMA.approx.sig
#' #              0.72           0.5461255         0.06159257               0.74
#'
#' #as.data.frame(
#' #   MetaAnalysisSimulations(
#' #     mean=0, sd=1, diff=0.5, GroupSize=10, type='n', Replications=50, Exp=5,
#' #     seed=456, alpha=0.05, FourGroup=TRUE, StdAdj=0, BlockEffect=0.5,
#' #     BlockStdAdj=0, StdExp=0, MAMethod='PM', returnES=FALSE))
#' #  AverageCliffd AverageCliffdvar AverageCliffdsig Averagephat ...
#' #1       0.27968       0.00683327             0.92     0.63984 ...
#' # as.data.frame(
#' #   MetaAnalysisSimulations(
#' #     mean=0, sd=1, diff=0.5, GroupSize=10, type='n', Replications=10, Exp=5,
#' #     seed=456, alpha=0.05, FourGroup=TRUE, StdAdj=0, BlockEffect=0.5,
#' #     BlockStdAdj=0, StdExp=0, MAMethod='PM', returnES=TRUE))
#' #Family NumExp GroupSize AveCliffd AveCliffdvar AveCliffdsig Avephat  ...
#' #     1      1         5        10        0.252  0.007423693    TRUE  ...
#' # Family NumExp GroupSize AveCliffd AveCliffdvar AveCliffdsig Avephat ...
#' #1     1      5        10     0.252  0.007423693         TRUE   0.626 ...

MetaAnalysisSimulations <-
  function(mean = 0,
           sd = 1,
           diff = 0.5,
           GroupSize = 10,
           type = "n",
           Replications = 50,
           Exp = 5,
           seed = 456,
           alpha = 0.05,
           FourGroup = FALSE,
           StdAdj = 0,
           BlockEffect = 0,
           BlockStdAdj = 0,
           StdExp = 0,
           MAMethod = "PM",
           returnES = FALSE,
           AlwaysTwoSidedTests = FALSE) {
    # Set up variables to hold the outcome of the simulations

    AverageCliffd <- 0 # Holds the average of Cliff's d
    AverageCliffdvar <-
      0 # Holds the average of the variance of Cliff's d
    AverageCliffdsig <-
      0 # Holds the average number of times the Cliff's d was assessed as significant using average hypothesis testing variance
    Averagephat <- 0 # Holds the average of phat
    Averagephatvar <- 0 # Holds the average of the variance of phat
    Averagephatsig <-
      0 # Holds the average number of times phat was assessed as signifciant using average hypothesis testing variance
    AveMDStd <-
      0 # Holds the mean of standardized effect sizes based in aggregating the MD
    AveMDStdvar <- 0
    AveMDStdsig <- 0


    Mean.StdMDAdjMA.exact <-
      0 # Holds the standardized effect size based on meta-analysis with the exact variance
    StdMDAdjMA.exact.var <- 0
    StdMDAdjMA.exact.sig <-
      0 # Holds the number of times the standardized effect size was significant with a meta-analysis using the exact variance

    Mean.StdMDAdjMA.approx <-
      0 # Holds the standardized effect size based on meta-analysis with the approximate variance
    StdMDAdjMA.approx.var <- 0
    StdMDAdjMA.approx.sig <-
      0 # Holds the number of times the standardized effect size was significant with a meta-analysis using the approximate variance

    MAMean.phat <-
      0 # Holds the mean phat effect size  based on a meta-analysis
    MAphat.sig <-
      0 # Holds the number of times the phat effect size was significant with a meta-analysis
    MAphat.var <- 0

    MAMean.Cliffd <-
      0 # Holds the mean Cliff's d effect size based on meta-analysis
    MACliffd.sig <-
      0 # Holds the number of times the Cliff's d effect size was significant with a meta-analysis
    MACliffd.var <- 0

    Mean.StdMDUnweighted <-
      0 # Holds the mean of the unweighted family level estimate of Cohen's d
    StdMDUnweighted.sig <-
      0 # Holds the proportion of significant unweighted family level estimates of Cohen's d
    StdMDUnweighted.var <- 0

    Mean.StdMDAdjUnweighted <-
      0 # Holds the mean of the unweighted family level, small-sample size adjusted,  estimate of the Cohen's d
    StdMDAdjUnweighted.sig <-
      0 # Holds the proportion of significant unweighted family level, small-sample size adjusted,  estimate of the Cohen's d
    StdMDAdjUnweighted.var <- 0

    Mean.StdMDMA.exact <-
      0 # Holds the mean SMD effect size based on meta-analysis
    StdMDMA.exact.var <-
      0 # Holds the mean variance of the StdMDMA.exact effect size
    StdMDMA.exact.sig <-
      0 # Holds the number of times the SMD effect size was significant with a meta-analysis

    Mean.StdMDMA.approx <-
      0 # Holds the mean SMD effect size based on meta-analysis
    StdMDMA.approx.var <-
      0 # Holds the nvariance of the StdMDMA.approx effect size
    StdMDMA.approx.sig <-
      0 # Holds the number of times the SMD effect size was significant with a meta-analysis


    # HedgesMA.Mean.Unweighted=0 # Holds the average of the small sample size adjusted StdMD
    Mean.HedgesMA <-
      0 # Holds the weighted average of the small sample size adjusted StdMD based on Hedges spcial analysis process
    Hedges.sig <-
      0 # Holds the proportion of times the small sample size adjusted StdMD was significant
    Hedges.var <- 0


    DataTable <- NULL

    FourGroupExp <- c(rep(FourGroup, Replications))



    for (i in 1:Replications) {
      if (!FourGroupExp[i]) {
        # 14-03-2022 parameter BlockStdadj changed to BlockStdAdj for consistency
        res <- NP2GMetaAnalysisSimulation(
          mean,
          sd,
          diff,
          GroupSize,
          Exp = Exp,
          type = type,
          seed = seed + i,
          alpha = alpha,
          StdAdj = StdAdj,
          StdExp = StdExp,
          MAMethod = MAMethod,
          AlwaysTwoSidedTests = AlwaysTwoSidedTests
        )
      }
      if (FourGroupExp[i]) {
        res <- NP4GMetaAnalysisSimulation(
          mean,
          sd,
          diff,
          GroupSize,
          Exp = Exp,
          type = type,
          seed = seed + i,
          alpha = alpha,
          StdAdj = StdAdj,
          BlockEffect = BlockEffect,
          BlockStdAdj = BlockStdAdj,
          StdExp = StdExp,
          MAMethod = MAMethod,
          AlwaysTwoSidedTests = AlwaysTwoSidedTests
        )
      }


      if (!returnES) {
        AverageCliffd <- AverageCliffd + res$AveCliffd
        AverageCliffdvar <- AverageCliffdvar + res$AveCliffdvar
        if (res$AveCliffdsig) {
          AverageCliffdsig <- AverageCliffdsig + 1
        }

        Averagephat <- Averagephat + res$Avephat
        Averagephatvar <- Averagephatvar + res$Avephatvar
        if (res$Avephatsig) {
          Averagephatsig <- Averagephatsig + 1
        }

        AveMDStd <- AveMDStd + res$AveMDStd
        AveMDStdvar <- AveMDStdvar + res$AveMDStdvar
        if (res$AveMDStdsig) {
          AveMDStdsig <- AveMDStdsig + 1
        }


        MAMean.phat <- MAMean.phat + res$MAphat
        MAphat.var <- MAphat.var + res$MAphatvar
        if (res$MAphatsig) {
          MAphat.sig <- MAphat.sig + 1
        }

        MAMean.Cliffd <- MAMean.Cliffd + res$MACliffd
        MACliffd.var <- MACliffd.var + res$MACliffdvar
        if (res$MACliffdsig) {
          MACliffd.sig <- MACliffd.sig + 1
        }

        Mean.StdMDUnweighted <-
          res$StdMDUnweighted + Mean.StdMDUnweighted
        StdMDUnweighted.var <-
          StdMDUnweighted.var + res$StdMDUnweightedvar
        if (res$StdMDUnweightedsig) {
          StdMDUnweighted.sig <- StdMDUnweighted.sig + 1
        }

        Mean.StdMDAdjUnweighted <-
          res$StdMDAdjUnweighted + Mean.StdMDAdjUnweighted
        if (res$StdMDAdjUnweightedsig) {
          StdMDAdjUnweighted.sig <- StdMDAdjUnweighted.sig + 1
        }
        StdMDAdjUnweighted.var <-
          StdMDAdjUnweighted.var + res$StdMDAdjUnweightedvar

        Mean.StdMDAdjMA.exact <-
          Mean.StdMDAdjMA.exact + res$StdMDAdjMAexact
        StdMDAdjMA.exact.var <-
          StdMDAdjMA.exact.var + res$StdMDAdjMAexactvar
        if (res$StdMDAdjMAexactsig) {
          StdMDAdjMA.exact.sig <- StdMDAdjMA.exact.sig + 1
        }

        Mean.StdMDAdjMA.approx <-
          Mean.StdMDAdjMA.approx + res$StdMDAdjMAapprox
        StdMDAdjMA.approx.var <-
          StdMDAdjMA.approx.var + res$StdMDAdjMAapproxvar
        if (res$StdMDAdjMAapproxsig) {
          StdMDAdjMA.approx.sig <- StdMDAdjMA.approx.sig + 1
        }


        Mean.StdMDMA.exact <- Mean.StdMDMA.exact + res$StdMDMAexact
        StdMDMA.exact.var <- StdMDMA.exact.var + res$StdMDMAexactvar
        if (res$StdMDMAexactsig) {
          StdMDMA.exact.sig <- StdMDMA.exact.sig + 1
        }


        Mean.StdMDMA.approx <- Mean.StdMDMA.approx + res$StdMDMAapprox
        StdMDMA.approx.var <-
          StdMDMA.approx.var + res$StdMDMAapproxvar
        if (res$StdMDMAapproxsig) {
          StdMDMA.approx.sig <- StdMDMA.approx.sig + 1
        }

        Mean.HedgesMA <- Mean.HedgesMA + res$HedgesMA.Weighted
        if (res$HedgesMA.Weightedsig) {
          Hedges.sig <- Hedges.sig + 1
        }
        Hedges.var <- Hedges.var + res$HedgesMA.Weightedvar
      } else {
        # Store the outcome from each replication
        DataTable <-
          tibble::as_tibble(dplyr::bind_rows(DataTable, dplyr::bind_cols(Family = i, res)))
      }
    }
    if (!returnES) {
      # Calculate averages.

      AverageCliffd <- AverageCliffd / Replications
      AverageCliffdvar <- AverageCliffdvar / Replications
      AverageCliffdsig <- AverageCliffdsig / Replications

      Averagephat <- Averagephat / Replications
      Averagephatvar <- Averagephatvar / Replications
      Averagephatsig <- Averagephatsig / Replications

      AveMDStd <- AveMDStd / Replications
      AveMDStdvar <- AveMDStdvar / Replications
      AveMDStdsig <- AveMDStdsig / Replications

      MAMean.phat <- MAMean.phat / Replications
      MAphat.var <- MAphat.var / Replications
      MAphat.sig <- MAphat.sig / Replications

      MAMean.Cliffd <- MAMean.Cliffd / Replications
      MACliffd.var <- MACliffd.var / Replications
      MACliffd.sig <- MACliffd.sig / Replications

      Mean.StdMDUnweighted <- Mean.StdMDUnweighted / Replications
      StdMDUnweighted.var <- StdMDUnweighted.var / Replications
      StdMDUnweighted.sig <- StdMDUnweighted.sig / Replications

      Mean.StdMDAdjUnweighted <-
        Mean.StdMDAdjUnweighted / Replications
      StdMDAdjUnweighted.var <- StdMDAdjUnweighted.var / Replications
      StdMDAdjUnweighted.sig <- StdMDAdjUnweighted.sig / Replications


      Mean.StdMDAdjMA.exact <- Mean.StdMDAdjMA.exact / Replications
      StdMDAdjMA.exact.var <- StdMDAdjMA.exact.var / Replications
      StdMDAdjMA.exact.sig <- StdMDAdjMA.exact.sig / Replications


      Mean.StdMDAdjMA.approx <- Mean.StdMDAdjMA.approx / Replications
      StdMDAdjMA.approx.var <- StdMDAdjMA.approx.var / Replications
      StdMDAdjMA.approx.sig <- StdMDAdjMA.approx.sig / Replications

      Mean.StdMDMA.exact <- Mean.StdMDMA.exact / Replications
      StdMDMA.exact.var <- StdMDMA.exact.var / Replications
      StdMDMA.exact.sig <- StdMDMA.exact.sig / Replications

      Mean.StdMDMA.approx <- Mean.StdMDMA.approx / Replications
      StdMDMA.approx.var <- StdMDMA.approx.var / Replications
      StdMDMA.approx.sig <- StdMDMA.approx.sig / Replications

      Mean.HedgesMA <- Mean.HedgesMA / Replications
      Hedges.var <- Hedges.var / Replications
      Hedges.sig <- Hedges.sig / Replications
    }

    if (!returnES) {
      outcome <- tibble::tibble(
        AverageCliffd,
        AverageCliffdvar,
        AverageCliffdsig,
        Averagephat,
        Averagephatvar,
        Averagephatsig,
        AveMDStd,
        AveMDStdvar,
        AveMDStdsig,
        MAMean.phat,
        MAphat.var,
        MAphat.sig,
        MAMean.Cliffd,
        MACliffd.var,
        MACliffd.sig,
        Mean.StdMDUnweighted,
        StdMDUnweighted.var,
        StdMDUnweighted.sig,
        Mean.StdMDAdjUnweighted,
        StdMDAdjUnweighted.var,
        StdMDAdjUnweighted.sig,
        Mean.HedgesMA,
        Hedges.var,
        Hedges.sig,
        Mean.StdMDAdjMA.exact,
        StdMDAdjMA.exact.var,
        StdMDAdjMA.exact.sig,
        Mean.StdMDAdjMA.approx,
        StdMDAdjMA.approx.var,
        StdMDAdjMA.approx.sig,
        Mean.StdMDMA.exact,
        StdMDMA.exact.var,
        StdMDMA.exact.sig,
        Mean.StdMDMA.approx,
        StdMDMA.approx.var,
        StdMDMA.approx.sig
      )
    } else {
      outcome <- DataTable
    }

    return(outcome)
  }

####################################################################################################################################################

# Function used to find the expected means,variance and effect sizes from different distributions

#' @title calculatePopulationStatistics
#' @description This helper function constructs the theoretical effect sizes and distribution statistics four (normal, lognormal, Laplace & gamma) given specific parameter values for the distributions and is used to support the calculation of population statistics for two and four group experiments.
#' @author Barbara Kitchenham and Lech Madeyski
#' @param mean The theoretical central location parameter for the distribution specified by the type parameter.
#' @param std The theoretical spread parameter for the distribution specified by the type parameter.
#' @param type String identifying the distribution, 'n' for normal, 'ln' for lognormal, 'lap' for Laplace, 'g' for Gamm
#' @return dataframe containing the expected standardized effect size, mean, variance,skewness and kurtosis statistics for samples from the specific distribution
#' @examples
#' reproducer:::calculatePopulationStatistics(mean=0, std=1, type='l')
#' #   RawMean RawVariance RawEffectSize RawSkewness RawKurtosis
#' #1 1.648721    4.670774      0.762874    6.184877    88.54343
#' reproducer:::calculatePopulationStatistics(mean=0, std=1, type='n')
#' #   RawMean RawVariance RawEffectSize RawSkewness RawKurtosis
#' # 1       0           1             0           0           3

calculatePopulationStatistics <- function(mean, std, type = "n") {
  if (type == "n") {
    # The expected values of a sample from the normal distribution
    RawMean <- mean
    RawVariance <- std ^ 2
    RawSkewness <- 0
    RawKurtosis <- 3
  }
  if (type == "l") {
    # The expected values of a sample from the lognormal distribution
    RawMean <- exp(mean + std ^ 2 / 2)
    RawVariance <- (exp(std ^ 2) - 1) * exp(2 * mean + std ^ 2)
    RawSkewness <- (exp(std ^ 2) + 2) * sqrt(exp(std ^ 2) - 1)
    RawKurtosis <-
      exp(4 * std ^ 2) + 2 * exp(2 * std ^ 2) + 3 * exp(2 * std ^ 2) - 3
  }
  if (type == "g") {
    # The expected values of a sample from the gamma distribution
    shape <- std
    rate <- mean
    RawMean <- shape / rate
    RawVariance <- shape / rate ^ 2
    RawSkewness <- 2 / sqrt(shape)
    RawKurtosis <- 6 / shape + 3
  }
  if (type == "lap") {
    # The expected values of a sample from the Laplace distribution

    location <- mean
    scale <- std
    RawMean <- location
    RawVariance <- 2 * scale ^ 2
    RawSkewness <- 0
    RawKurtosis <- 6
  }
  RawEffectSize <- RawMean / sqrt(RawVariance)
  output <-
    tibble::tibble(RawMean,
                   RawVariance,
                   RawEffectSize,
                   RawSkewness,
                   RawKurtosis)
  return(output)
}

#' @title RandomizedDesignEffectSizes
#' @description This function creates the theoretical effect sizes for data from one of four different distributions for specified parameter values for the distribution specified by the type parameter. It assumes there are two samples, one corresponding to a control group and the other to the treatment group. It returns the theoretical effect sizes for a fully randomized experiment.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export RandomizedDesignEffectSizes
#' @param m1 The theoretical mean for the control group
#' @param std1 The theoretical variance for the control group
#' @param m2 The theoretical mean for the treatment group
#' @param std2 The theoretical variance for the treatment group
#' @param type String identifying the distribution, 'n' for normal, 'ln' for lognormal, 'lap' for Laplace, 'g' for Gamma
#' @return dataframe containing the expected values of the unstandardized mean difference effect size, the pooled within group variance, the standardized mean difference effect size and the point bi-serial correlation.
#' @examples
#' RandomizedDesignEffectSizes(m1=0, std1=1, m2=1, std2=3, type = 'n')
#' #  ES Var     StdES      rPBS
#' #1  1   5 0.4472136 0.2182179
#' RandomizedDesignEffectSizes(m1=0, std1=1, m2=1, std2=3, type = 'l')
#' #        ES       Var     StdES        rPBS
#' #1 243.0432 242552663 0.0156056 0.007802562
#'  RandomizedDesignEffectSizes(m1=0, std1=1, m2=0.266, std2=1, type = 'l')
#' #          ES      Var     StdES       rPBS
#' # 1 0.5024232 6.310995 0.1999957 0.09950162
RandomizedDesignEffectSizes <-
  function(m1, std1, m2, std2, type = "n") {
    G1.results <- calculatePopulationStatistics(m1, std1, type = type)
    G2.results <- calculatePopulationStatistics(m2, std2, type = type)
    ES <- G2.results$RawMean - G1.results$RawMean
    Var <- (G2.results$RawVariance + G1.results$RawVariance) / 2
    StdES <- ES / sqrt(Var)
    rPBS <- StdES / sqrt(StdES ^ 2 + 4)
    output <- tibble::tibble(ES, Var, StdES, rPBS)
    return(output)
  }

#' @title RandomizedBlockDesignEffectSizes
#' @description This function finds the theoretical effect sizes for a four-group randomized block experiments assuming one of four different underlying distributions specified by the type parameter. The design assumes two blocks each comprising a control and treatment group. If required a fixed Blocking effect is added to the mean for Block 2.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export RandomizedBlockDesignEffectSizes
#' @param m1 The theoretical mean for the control group in Block 1
#' @param std1 The theoretical variance for the control group in Block 1
#' @param m2 The theoretical mean for the treatment group in Block 1
#' @param std2 The theoretical variance for the treatment group in Block 1
#' @param m3 The theoretical mean for the control group in Block 2
#' @param std3 The theoretical variance for the control group in Block 2
#' @param m4 The theoretical mean for the treatment group in Block 2
#' @param std4 The theoretical variance for the treatment group in Block 2
#' @param BE A fixed block effect to be added to the Block 2 mean values.
#' @param type String identifying the distribution, 'n' for normal, 'ln' for lognormal, 'lap' for Laplace, 'g' for Gamma
#' @return dataframe holing the expected unstandardized mean difference effect size, the pooled within group variance, the standardized effect size and the point bi-serial correlation.
#' @examples
#' RandomizedBlockDesignEffectSizes(m1=0,std1=1,m2=1,std2=1,m3=0,std3=1,m4=1,
#'   std4=1,BE = 1,type = 'n')
#' # ES Var StdES      rPBS
#' #1  1   1     1 0.4472136
#' RandomizedBlockDesignEffectSizes(m1=0,std1=1,m2=1,std2=1,m3=0,std3=1,m4=1,
#'   std4=1,BE = 1,type = 'l')
#' #        ES      Var     StdES      rPBS
#' #1 5.266886 82.17791 0.5810004 0.2789675
#' RandomizedBlockDesignEffectSizes(
#'   m1=0,std1=1,m2=0.266,std2=1,m3=0,std3=1,m4=0.266,std4=1,BE = 0,type = 'l')
#' #        ES      Var     StdES       rPBS
#' #1 0.5024232 6.310995 0.1999957 0.09950162
RandomizedBlockDesignEffectSizes <-
  function(m1,
           std1,
           m2,
           std2,
           m3,
           std3,
           m4,
           std4,
           BE = 0,
           type = "n") {
    # Find the basic statistics for each group
    G1.results <- calculatePopulationStatistics(m1, std1, type = type)
    G2.results <- calculatePopulationStatistics(m2, std2, type = type)

    # 12-03-2022 Block effect incorrectly handled for all types. Revised so that the block effect is added to the shape for Gamma variables and
    # the mean for all other data types

    if (type == "g") {
      G3.results <-
        calculatePopulationStatistics(m3, std3 + BE, type = type)
      G4.results <-
        calculatePopulationStatistics(m4, std4 + BE, type = type)
    } else {
      G3.results <-
        calculatePopulationStatistics(m3 + BE, std3, type = type)
      G4.results <-
        calculatePopulationStatistics(m4 + BE, std4, type = type)
    }

    # Calculate the expected unstandardized effect size allowing for the experimental design
    ES <-
      (G2.results$RawMean - G1.results$RawMean + G4.results$RawMean - G3.results$RawMean) / 2
    # Calculate the within groups pooled variance
    Var <-
      (
        G2.results$RawVariance + G1.results$RawVariance + G3.results$RawVariance + G4.results$RawVariance
      ) / 4
    # Calculate the standardized mean difference effect size
    StdES <- ES / sqrt(Var)
    # Calculate the point bi-serial effect size
    rPBS <- StdES / sqrt(StdES ^ 2 + 4)
    output <- tibble::tibble(ES, Var, StdES, rPBS)
    return(output)
  }


#' @title crossoverResidualAnalysis
#' @description This function analyses one or more crossover experiments where each experiment can be either a two group or four group experiment as specified by Type parameter. The file parameter includes a dataset comprising one or more experiments as defined by the ExperimentNames parameter and each experiment includes values for every output variable defined in the Metrics parameter. After being analysed using the linear modeling lmer function of the lme4 package, the residuals are assessed for normality based on  the Anderson-Darling test. Warning 1. This function should only be used with data sets that include one or more individual crossover experiments in the format used in the datasets reported in the reproducer package that were used in the paper B. Kitchenham, L. Madeyski, G. Scanniello, and C. Gravino, “The importance of the correlation between results from individual participants in crossover studies,” IEEE Transactions in SoftwareEngineering, 2021 (Accepted for Publication). [Online]. Available: https://doi.org/10.1109/TSE.2021.30. Warning 2. The lmer function assumes that when experiments include multiple data from the same participant, the correlation between measures for the same participant will be positive. If there is no positive correlation, the function will deliver the warning: boundary (singular) fit: see ?isSingular. This does not mean the analysis has failed. It means that the within participant and between participant variance are set to the same value.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export crossoverResidualAnalysis
#' @param file The dataset to be analysed.
#' @param StudyID A character string used to identify the origin of the dataset
#' @param ExperimentNames A vector of one or more strings variables identifying each experiment in the file.
#' @param Type A vector of string variables identifying the design for each experiment. Each element should have the value '2G' or '4G'
#' @param Metrics A vector of string variables identifying the variables to be analysed.
#' @return The results of analysing the residuals for each experiment and metric using the Anderson-Darling test (ADpval) and the number of outliers (NUmOut).
#' @examples
#' File=reproducer::KitchenhamEtAl.CorrelationsAmongParticipants.Scanniello15EMSE

#' crossoverResidualAnalysis(
#'   File,StudyID='S1',ExperimentNames=c('USB2'),Type=c('4G'),
#'   Metrics=c('Correctness','Time','Efficiency'))

#' #  Study  Exp     Metrics  N ADpval NumOut
#' #1    S1 USB2 Correctness 24 0.0846      2
#' #2    S1 USB2        Time 24 0.0448      1
#' #3    S1 USB2  Efficiency 24  0.365      4

crossoverResidualAnalysis <-
  function(file,
           StudyID,
           ExperimentNames,
           Type,
           Metrics) {
    NumMets <- length(Metrics)
    NumExp <- length(ExperimentNames)

    ExpData <-
      reproducer::ExtractExperimentData(
        file,
        ExperimentNames = ExperimentNames,
        idvar = "ParticipantID",
        timevar = "Period",
        ConvertToWide = FALSE
      )

    table <- NULL
    for (j in 1:NumMets) {
      # Perform a linear mixed model analysis for each experiment and each metric

      for (i in 1:NumExp) {
        # Construct a row for the outcome of the linear model analysis depending on the experiment type and analyse the residuals

        LinearModel <- doLM(ExpData[[i]], Metrics[j], Type[i])

        Residuals <- summary(LinearModel)$residuals

        N <- length(Residuals) / 2

        row <-
          data.frame(AnalyseResiduals(Residuals, ExperimentNames[i]))


        table <-
          data.frame(rbind(
            table,
            cbind(
              Study = StudyID,
              Exp = ExperimentNames[i],
              Metrics = Metrics[j],
              N,
              ADpval = signif(row$AndersonDarling,
                              3),
              NumOut = row$NumOut
            )
          ))
      }
    }

    return(table)
  }

#' @title doLM
#' @description This helper function is called by the function 'crossoverResidualAnalysis' to perform either an AB/BA crossover analysis or a four-group crossover on the data set defined by the parameter DataSet depending on the value of the Type parameter.
#' @author Barbara Kitchenham and Lech Madeyski
#' @param DataSet The dataset to be analysed.
#' @param Metric The name of the variable to be analysed
#' @param Type Defines the experimental design
#' @return The data analysis results provided by the lmer function of the lme4 package

doLM <- function(DataSet, Metric, Type) {
  form1 <-
    paste(Metric,
          "~Period+Treatment+System+CrossOverID+(1|ParticipantID)")
  form2 <- paste(Metric, "~Period+Treatment+(1|ParticipantID)")

  if (Type == "4G") {
    LinearModel <- lme4::lmer(form1, data = DataSet)
  }

  if (Type == "2G") {
    LinearModel <- lme4::lmer(form2, data = DataSet)
  }

  return(LinearModel)
}


#
# 6 functions to generate results in our paper "Recommendations for Analyzing Software Engineering Experiments":
#

#' @title calculateMAType1Error
#' @description The function simulates multiple five group families of either two-group or four-group experiments and estimates the Type1 Error rate obtained after synthesizing  the analysis results obtained from the experiments in each family. The Type1 Error is estimated as the percentage of families for which the overall mean of the five experiments was significantly different from zero. The experiment data may be one of four different type: Normal, Log-normal, Gamma or Laplace. The simulations can be repeated for different sample sizes depending on the parameter N. The output is a table of values identifying three observed effect size estimates (Cliff's d, PHat and StdMD) and their related type 1 error rates for each set of simulated families. The synthesis method for all three effect sizes is based on calculating the overall mean and variance for the family of experiments, and then using those values to calculate the effect size variance and its variance. This function supports the production of the values reported in data tables in the paper "Recommendations for Analyzing Small Sample Size Software Engineering Experiments" and its Supplementary Material.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export calculateMAType1Error
#' @param mean This is the mean value of the control and treatment group(s) used in the simulations of each experiment of each family for simulations of a specified sample size (default 0).
#' @param sd This is the standard deviation  value of the control group(s) and treatment group(s) used in the simulations of each experiment of each family for simulations of a specified sample size (default 1).
#' @param N This specifies the sample sizes per group that will be used in each set of simulations (default c(5,10,15,20,30,40)).
#' @param reps The number of families simulated for each sample size.
#' @param type This specifies the distribution of the data samples that will be simulated. Options ae "n" for Normal, "l", for Log-normal,'g" for Gamma, "lap" for LaPlace (default "n").
#' @param seed A seed for the simulations (default 123).
#' @param Experiments The number of experiments in each family (default 5).
#' @param FourG If FourG is FALSE the individual experiments in each family will be two-group experiments, otherwise the individual experiments will be four-group families (default FALSE).
#' @param StdAdj Used to introduce variance heterogeneity for Laplace and Normal samples (default 0).
#' @param Blockmean Used to set a fixed block effect for four-group experiments (default 0).
#' @param BlockStdAdj Not used (default 0).
#' @param StdExp Used to introduce heterogeneity among families of experiments (default 0).
#' @param MAMethod Not used (default "PM").
#' @param alpha The significance level for statistical tests (default 0.05).
#' @return Design. Specifies the type of experiment 2G or 4G, the sample distribution (n,l,g,lap), and whether variance heterogeneity was added (het)
#' @return BEIncluded. Specifies whether or not a block effect was introduced. Always set to "No" for two-group experiments.
#' @return GrpSize. Specifies the size of each group in the individual experiments.
#' @return ObsPHat. The average of the average Phat value found for each family in the set of simulations.
#' @return ObsCliffd. The average of the average Cliffd value found for each family in the set of simulations.
#' @return ObsStdES. The average of StdMD calculated for each family in the set of simulations.
#' @return PHatType1ER. The percentage of the simulations, for a specific group size, for which the overall Phat estimate was significantly different from zero at the nominated alpha level.
#' @return CliffdType1ER. The percentage of the simulations, for a specific group size, for which the overall Cliff's d estimate was significantly different from zero at the nominated alpha level.
#' @return StdMDType1ER. The percentage of the simulations, for a specific group size, for which the overall StdMD estimate was significantly different from zero at the nominated alpha level.
#' @examples
#'# as.data.frame(calculateMAType1Error(mean=0,sd=1,N=c(5,10),reps=10,type="n",Experiments=5,
#'#  FourG=FALSE,StdAdj=0,Blockmean=0,seed=123))
#' #  Design BEIncluded GrpSize ObsPHat ObsCliffd     ObsStdES PHatType1ER CliffdType1ER StdMDType1ER
#' #1   2G_n         No       5  0.4848   -0.0304 -0.054156883           0             0          0.0
#' #2   2G_n         No      10  0.5036    0.0072  0.002888142           0             0          0.1
#'
#'#as.data.frame(calculateMAType1Error(mean=0,sd=1,N=c(5,10),reps=10,type="l",Experiments=5,
#'#  FourG=FALSE,StdAdj=0,Blockmean=0,seed=123))
#' #   Design BEIncluded GrpSize ObsPHat ObsCliffd    ObsStdES PHatType1ER CliffdType1ER StdMDType1ER
#' #1   2G_l         No       5  0.4848   -0.0304 -0.02789656           0             0          0.0
#' #2   2G_l         No      10  0.5036    0.0072  0.06473696           0             0          0.2
#'
#'#as.data.frame(calculateMAType1Error(mean=0,sd=1,N=c(5,10),reps=10,type="n",Experiments=5,
#'#  FourG=TRUE,StdAdj=0.5,Blockmean=0.5,seed=123))
#' #   Design BEIncluded GrpSize ObsPHat ObsCliffd   ObsStdES PHatType1ER CliffdType1ER StdMDType1ER
#' #  1 4G_n_het        Yes       5  0.5108    0.0216 0.01361820           0             0          0.1
#' #  2 4G_n_het        Yes      10  0.5069    0.0138 0.01700672           0             0          0.0
#'
#'as.data.frame(calculateMAType1Error(mean=0,sd=1,N=c(5,10),reps=5,type="l",Experiments=5,
#'  FourG=TRUE,StdAdj=0,Blockmean=0.5,seed=123))
#'  #Results for reps=10
#' #   Design BEIncluded GrpSize ObsPHat ObsCliffd   ObsStdES PHatType1ER CliffdType1ER StdMDType1ER
#' #1   4G_l        Yes       5  0.5108    0.0216 0.07578257           0             0          0.2
#' #2   4G_l        Yes      10  0.5072    0.0144 0.04839936           0             0          0.0

calculateMAType1Error = function(mean = 0,
                                 sd = 1,
                                 N = c(5, 10, 15, 20, 30, 40),
                                 reps,
                                 type = "n",
                                 seed = 123,
                                 Experiments = 5,
                                 FourG = FALSE,
                                 StdAdj = 0,
                                 Blockmean = 0,
                                 BlockStdAdj = 0,
                                 StdExp = 0,
                                 MAMethod = "PM",
                                 alpha = 0.05) {
  Type1ErrorTable = NULL

  Design = "2G"

  if (FourG) {
    Design = "4G"
  }

  Design = paste(Design, type, sep = "_")

  if (StdAdj > 0) {
    Design = paste(Design, "het", sep = "_")
  }

  BEIncluded = "No"
  if (Blockmean > 0 & FourG) {
    BEIncluded = "Yes"
  }


  NumSamples = length(N)

  for (i in 1:NumSamples) {
    Out = reproducer::MetaAnalysisSimulations(
      mean = mean,
      sd = sd,
      diff = 0,
      GroupSize = N[i],
      type = type,
      Replications = reps,
      Exp = Experiments,
      FourGroup = FourG,
      StdAdj = StdAdj,
      BlockEffect = Blockmean,
      BlockStdAdj = BlockStdAdj,
      StdExp = StdExp,
      MAMethod = MAMethod,
      returnES = TRUE,
      seed = seed,
      alpha = 0.05
    )


    PHatType1ER = mean(as.numeric(Out$Avephatsig))

    CliffdType1ER = mean(as.numeric(Out$AveCliffdsig))

    AveMDStdType1ER = mean(as.numeric(Out$AveMDStdsig))


    temp = dplyr::bind_cols(
      Design = Design,
      BEIncluded = BEIncluded,
      GrpSize = as.character(N[i]),
      ObsPHat = mean(Out$Avephat),
      ObsCliffd = mean(Out$AveCliffd),
      ObsStdES = mean(Out$AveMDStd),
      PHatType1ER =  PHatType1ER,
      CliffdType1ER = CliffdType1ER,
      StdMDType1ER = AveMDStdType1ER
    )

    Type1ErrorTable = dplyr::bind_rows(Type1ErrorTable, temp)
  }

  return(Type1ErrorTable)
}

#################################################################################################################################

#' @title calculateMABias
#' @description The function simulates multiple five group families of either two-group or four-group experiments and estimates the power, individual estimate error, and the small sample bias obtained after synthesizing the analysis results obtained from the experiments in each family. The power is estimated as the percentage of families for which the overall mean of the five experiments was significantly different from zero. The experiment data may be one of four different type: Normal, Log-normal, Gamma or Laplace. The simulations can be repeated for different mean differences between the control mean and treatment mean depending on the parameter diff. The output is a table of values identifying the observed values of three effect sizes: Cliff's d, PHat and StdMD, estimate error and their related small sample bias and power for each set of simulated families. The synthesis method for all the effect sizes is based on calculating the overall mean and variance for experiments in each family and then using those values to calculate the overall effect size and its variance. This function supports the production of the values reported in data tables in the paper "Recommendations for Analyzing Small Sample Size Software Engineering Experiments" and its Supplementary Material.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export calculateMABias
#' @param mean This is the mean value of the control and treatment group(s) used in the simulations of each experiment of each family for simulations of a specified sample size (default 0).
#' @param sd This is the standard deviation  value of the control group(s) and treatment group(s) used in the simulations of each experiment of each family for simulations of a specified sample size (default 1).
#' @param N This specifies the sample size per group that will be used in each set of simulations.
#' @param reps The number of families simulated for each sample size.
#' @param diff This specifies the difference between the control and treatment that will be used in each set of simulations. It must always have three values representing small, medium and large values (default c(0.2, 0.5, 0.8)).
#' @param Experiments The number of experiments in each family (default 5).
#' @param Expected.StdMD This defines the expected value of the overall average StdMD for each mean difference (default c(0.2,0.5,0.8)).
#' @param Expected.PHat This defines the expected population value of the overall average Phat for each mean difference (default c(0.556,0.638,0.714)).
#' @param type This specifies the distribution of the data samples that will be simulated. Options ae "n" for Normal, "l", for Log-normal,'g" for Gamma, "lap" for LaPlace (default "n").
#' @param FourG If FourG is FALSE (default) the individual experiments in each family will be two-group experiments, otherwise the individual experiments will be four-group families.
#' @param seed A seed for the simulations (default 123).
#' @param StdAdj Used to introduce variance heterogeneity for Laplace and Normal samples (default 0).
#' @param Blockmean Used to set a fixed block effect for four-group experiments (default 0).
#' @param StdExp Used to introduce heterogeneity among families of experiments (default 0).
#' @param MAMethod Not used (default "PM").
#' @param alpha The significance level for statistical tests (default 0.05).
#' @return Design. Specifies the type of experiment 2G or 4G, the sample distribution (n,l,g,lap), and whether variance heterogeneity was added (het)
#' @return BEIncluded. Specifies whether or not a block effect was introduced. Always set to "No" for two-group experiments.
#' @return GrpSize. Specifies the size of each group in the individual experiments.
#' @return Diff. The size of the difference between the control and treatment converted to an ordinal scale (Small, Medium, Large)
#' @return NPBias The relative difference between the average of the observed values of either Cliff's d or centralised PHat and the population value
#' @return StdMDBias. The relative difference between the average of the observed values of StdMDBias and the theoretical value
#' @return NPMdMRE The median of the absolute relative difference between the observed values of either Cliff's d or centralised PHat and the theoretical value for each experiment.
#' @return StdMDMdMRE The median of the absolute relative difference between the observed values of StdMD and the population value for each experiment.
#' @return ObsPHat. The average of the average Phat value found for each family in the set of simulations.
#' @return ObsCliffd. The average of the average Cliffd value found for each family in the set of simulations.
#' @return ObsStdES. The average of StdMD calculated for each family in the set of simulations.
#' @return PHatPower. The percentage of the simulations, for a specific mean difference, for which the overall Phat estimate was significantly different from zero at the nominated alpha level using one-sided tests.
#' @return CliffdPower. The percentage of the simulations, for a specific mean difference, for which the overall Cliff's d estimate was significantly different from zero at the nominated alpha level using one-sided tests.
#' @return StdMDPower. The percentage of the simulations, for a specific mean difference, for which the overall StdMD estimate was significantly different from zero at the nominated alpha level using one-sided tests.
#' @examples
#'# as.data.frame(calculateMABias(mean=0,sd=1,N=10,diff=c(0.2,0.5,0.8), Experiments=5,reps=10,
#'# Expected.StdMD=c(0.2,0.5,0.8), Expected.PHat=c(0.556,0.638,0.714), type="n",FourG=FALSE,
#'# seed= 123, StdAdj = 0, Blockmean=0, StdExp=0))
#' #  Design Blockmean GrpSize   Diff     NPBias  StdMDBias   NPMdMRE StdMDMdMRE ObsPHat ObsCliffd..
#' #1   2G_n        No      10  Small 0.09285714 0.02606704 0.8928571  1.0741432  0.5612    0.1224..
#' #2   2G_n        No      10 Medium 0.03768116 0.01740262 0.2391304  0.4171896  0.6432    0.2864..
#' #3   2G_n        No      10  Large 0.03738318 0.01523651 0.2009346  0.2490287  0.7220    0.4440..
#' #  PHatPower CliffdPower StdESPower
#' #1       0.2         0.2        0.3
#' #2       0.7         0.7        0.7
#' #3       1.0         1.0        1.0
#'as.data.frame(calculateMABias(mean=0,sd=1,N=10,diff=c(0.2,0.5,0.8), Experiments=5,reps=4,
#'  Expected.StdMD=c(0.2,0.5,0.8), Expected.PHat=c(0.556,0.638,0.714), type="n",FourG=TRUE,
#'  seed= 123,StdAdj = 0.5,Blockmean=0.5,StdExp=0))
#'  #Results for reps=10
#' #    Design Blockmean GrpSize   Diff     NPBias  StdMDBias   NPMdMRE StdMDMdMRE ObsPHat ObsClif..
#' #1 4G_n_het       Yes      10  Small -0.1321429 -0.1372277 0.6696429  0.4698935  0.5486  0.0972..
#' #2 4G_n_het       Yes      10 Medium -0.1869565 -0.1882479 0.2318841  0.1472392  0.6122  0.2244..
#' #3 4G_n_het       Yes      10  Large -0.1864486 -0.2010029 0.1612150  0.1531253  0.6741  0.3482..
#' #  PHatPower CliffdPower StdESPower
#' #1       0.4         0.4        0.4
#' #2       0.9         0.9        0.8
#' #3       1.0         1.0        1.0


calculateMABias = function(mean = 0,
                           sd = 1,
                           N,
                           reps,
                           diff = c(0.2, 0.5, 0.8),
                           Experiments = 5,
                           Expected.StdMD = c(0.2, 0.5, 0.8),
                           Expected.PHat = c(0.556, 0.638, 0.714),
                           type = "n",
                           FourG = FALSE,
                           seed = 223,
                           StdAdj = 0,
                           Blockmean = 0,
                           StdExp = 0,
                           MAMethod = "PM",
                           alpha = 0.05) {
  MdMRETable = NULL

  NumESizes = length(diff)


  Design = "2G"

  if (FourG) {
    Design = "4G"
  }

  Design = paste(Design, type, sep = "_")

  DiffGen = c("Small", "Medium", "Large")

  if (StdAdj > 0) {
    Design = paste(Design, "het", sep = "_")
  }

  BlockMeanInfo = "No"
  if (Blockmean > 0 & FourG)
    BlockMeanInfo = "Yes"

  for (i in 1:NumESizes) {
    Out = reproducer::MetaAnalysisSimulations(
      mean = mean,
      sd = sd,
      diff = diff[i],
      GroupSize = N,
      type = type,
      Replications = reps,
      Exp = Experiments,
      FourGroup = FourG,
      seed = seed,
      alpha = alpha,
      StdAdj = StdAdj,
      BlockEffect = Blockmean,
      BlockStdAdj = 0,
      StdExp = StdExp,
      MAMethod = MAMethod,
      returnES = TRUE
    )

    CentralPhat = Out$Avephat - 0.5

    ExpectedCentralPhat = Expected.PHat[i] - 0.5

    Expected.Cliffd = (Expected.PHat[i] - 0.5) * 2

    NPBias = (mean(CentralPhat) - ExpectedCentralPhat) / ExpectedCentralPhat

    StdMDBias = (mean(Out$AveMDStd) - Expected.StdMD[i]) / Expected.StdMD[i]

    NPMdMRE = median(abs((Out$AveCliffd - Expected.Cliffd) / Expected.Cliffd))

    StdMDMdMRE = median(abs((Out$AveMDStd - Expected.StdMD[i]) / Expected.StdMD[i]))

    temp = dplyr::bind_cols(
      Design = Design,
      Blockmean = BlockMeanInfo,
      GrpSize = as.character(N),
      Diff = DiffGen[i],
      NPBias = NPBias,
      StdMDBias =  StdMDBias,
      NPMdMRE = NPMdMRE,
      StdMDMdMRE = StdMDMdMRE,
      ObsPHat = mean(Out$Avephat),
      ObsCliffd = mean(Out$AveCliffd),
      ObsStdES = mean(Out$AveMDStd),
      PHatPower = mean(Out$Avephatsig),
      CliffdPower = mean(Out$AveCliffdsig),
      StdESPower = mean(Out$AveMDStdsig)
    )

    MdMRETable = tibble::tibble(dplyr::bind_rows(MdMRETable, temp))
  }
  return(MdMRETable)

}

######################################################################################################################################
#' @title calculate2GType1Error
#' @description The function simulates  multiple two-group experiments and estimates the Type1 Error rate obtained from the set of simulated experiments. The Type1 Error is estimated as the percentage of experiments for which the mean the experiment was significantly different from zero at the 0.05 significance level using two-sided tests. The experiment data may be one of four different type: Normal, Log-normal, Gamma or Laplace. The output is a set of values identifying three observed effect size estimates (Cliff's d, PHat and StdMD) and their related type 1 error rates. This function supports the production of the values reported in data tables in the paper "Recommendations for Analyzing Small Sample Size Software Engineering Experiments" and its Supplementary Material.

#' @author Barbara Kitchenham and Lech Madeyski
#' @export  calculate2GType1Error
#' @param mean This is the mean value of the control and treatment group(s) used in the simulations of each experiment for simulations of a specified sample size (default 0).
#' @param sd This is the standard deviation  value of the control group(s) and treatment group(s) used in the simulations of each experiment of each family for simulations of a specified sample size (default 1).
#' @param N This specifies the sample size per group that will be used in each set of simulations (default 5).
#' @param reps The number of experiments to simulated.
#' @param type This specifies the distribution of the data samples that will be simulated. Options ae "n" for Normal, "l", for Log-normal,'g" for Gamma, "lap" for LaPlace (default "n").
#' @param seed A seed for the simulations (default 123).
#' @param StdAdj Used to introduce variance heterogeneity for Laplace and Normal samples(default 0).
#' @return Design. Specifies the type of experiment 2G or 4G, the sample distribution (n,l,g,lap), and whether variance heterogeneity was added (het)
#' @return GrpSize. Specifies the size of each group in the individual experiments.
#' @return ObsPHat. The average Phat values found in the set of simulations.
#' @return ObsCliffd. The average Cliffd values found  in the set of simulations.
#' @return ObsStdES. The average of StdMD values found in the set of simulations.
#' @return PHatType1ER. The proportion of the simulations for which the Phat estimate was significantly different from zero at the nominated alpha level.
#' @return CliffdType1ER. The  proportion of the simulations for which the  Cliff's d estimate was significantly different from zero at the nominated alpha level.
#' @return StdMDType1ER. The  proportion of the simulations for which the StdMD estimate was significantly different from zero at the nominated 0.05 significance level.

#' @examples

#'calculate2GType1Error(mean=1,sd=3,N=10,reps=100,type="g",seed=3256,StdAdj = 0)
#'# A tibble: 1 x 8
#'#   Design GrpSize ObsPHat ObsCliffd ObsStdES PHatType1ER CliffdType1ER StdESType1ER
#'#   <chr>  <chr>     <dbl>     <dbl>    <dbl>       <dbl>         <dbl>        <dbl>
#'# 1 2G_g   10        0.498   -0.0034 -0.00464        0.02          0.01         0.02


calculate2GType1Error = function(mean = 0,
                                 sd = 1,
                                 N = 10,
                                 reps,
                                 type = "n",
                                 seed = 123,
                                 StdAdj = 0) {
  Type1ErrorTable = NULL

  Out1 = reproducer::RandomExperimentSimulations(
    mean = mean,
    sd = sd,
    diff = 0,
    N = N,
    reps = reps,
    type = type,
    seed = seed,
    StdAdj = StdAdj,
    returnData = TRUE
  )

  Design = "2G"
  Design = paste(Design, type, sep = "_")

  if (StdAdj > 0) {
    Design = paste(Design, "het", sep = "_")
  }

  Type1ErrorTable = tibble::tibble(
    dplyr::bind_cols(
      Design = Design,
      GrpSize = as.character(N),
      ObsPHat = mean(Out1$PHat),
      ObsCliffd = mean(Out1$Cliffd),
      ObsStdES = mean(Out1$StdES),
      PHatType1ER = mean(Out1$PHatSig),
      CliffdType1ER = mean(Out1$CliffdSig),
      StdESType1ER = mean(Out1$ESSig)
    )
  )


  return(Type1ErrorTable)
}

###############################################################################################################################
#' @title calculate2GBias
#' @description The function simulates two-group experiments and estimates the power, individual estimate error, and the small sample bias obtained obtained from the set of simulated experiments. The set of simulations for a specific mean difference are repeated for three different values of the difference between the treatment and control groups specified by the parameter "diff". The power is estimated as the percentage of experiments for which the  mean of the experiment was significantly different from zero. The experiment data may be one of four different type: Normal, Log-normal, Gamma or Laplace. The output is a table of values identifying the observed values of three effect sizes: Cliff's d, PHat and StdMD, estimate error and their related small sample bias and power for each set of simulated experiments. This function supports the production of the values reported in data tables in the paper "Recommendations for Analyzing Small Sample Size Software Engineering Experiments" and its Supplementary Material.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export calculate2GBias
#' @param mean This is the mean value of the control and treatment group(s) used in the simulations of each experiment for simulations of a specified sample size and mean difference (default 0).
#' @param sd This is the standard deviation  value of the control group(s) and treatment group(s) used in the simulations of each experiment of each family for simulations of a specified sample size (default 1).
#' @param N This specifies the sample size per group that will be used in each set of simulations.
#' @param reps The number of experiments simulated for each mean difference.
#' @param diff This specifies the mean difference between the control and treatment that will be used in each set of simulations. It must always have three values representing small, medium and large differences (default c(0.2, 0.5, 0.8)).
#' @param Expected.StdMD This defines the theoretical value of the average StdMD obtained from the simulations for each mean difference. (default c(0.2, 0.5, 0.8))
#' @param Expected.PHat This defines the expected population value of the average Phat obtained from the simulations for each mean difference (default c(0.556,0.638,0.714)).
#' @param type This specifies the distribution of the data samples that will be simulated. Options ae "n" for Normal, "l", for Log-normal,'g" for Gamma, "lap" for LaPlace (default "n").
#' @param seed A seed for the simulations (default 123).
#' @param StdAdj Used to introduce variance heterogeneity for Laplace and Normal samples (default 0).
#' @return Design. Specifies the type of experiment, the sample distribution (n,l,g,lap), and whether variance heterogeneity was added (het)
#' @return GrpSize. Specifies the size of each group in the simulated experiments.
#' @return Diff. The size of the difference between the control and treatment converted to an ordinal scale (Small, Medium, Large)
#' @return NPBias The relative difference between the average of the observed values of either Cliff's d or centralised PHat and the population value
#' @return StdMDBias. The relative difference between the average of the observed values of StdMDBias and the theoretical value
#' @return NPMdMRE The median of the absolute relative difference between the observed values of either Cliff's d or centralised PHat and the theoretical value for each experiment.
#' @return StdMDMdMRE The median of the relative difference between the observed values of StdMD and the population value for each experiment.
#' @return ObsPHat. The average of the  Phat values found  in the set of simulations.
#' @return ObsCliffd. The average of the  Cliffd values found  in the set of simulations.
#' @return ObsStdES. The average of StdMD values found in the set of simulations.
#' @return PHatPower. The percentage of the simulations, for a specific mean difference, for which the Phat estimate was significantly different from zero at the 0.05 alpha level based on one-sided tests.
#' @return CliffdPower. The percentage of the simulations, for a specific mean difference, for which the  Cliff's d estimate was significantly different from zero at the 0.05 alpha level based on one-sided tests.
#' @return StdMDPower. The percentage of the simulations, for a specific mean difference, for which the StdMD estimate was significantly different from zero at the 0.05 alpha level based on one-sided tests.
#' @examples
#'# as.data.frame(calculate2GBias(mean=0,sd=1,diff=c(0.2,0.5,0.8),Expected.StdMD=c(0.157,0.392,0.628),
#'#  Expected.PHat=c(0.544,0.609,0.671), N=5,reps=50, type="n", seed=523, StdAdj =0.5 ))
#'# Results for reps=100 (due to NOTE "Examples with CPU (user + system) or elapsed time > 5s"):
#'#    Design GrpSize   Diff        NPBias  StdMDBias  NPMdMRE StdMDMdMRE ObsPHat ObsCliffd  ObsSt..
#'# 1 2G_n_het       5  Small -6.308085e-16 0.07088601 3.272727  3.2700082  0.5440    0.0880 0.168..
#'# 2 2G_n_het       5 Medium  3.486239e-02 0.09914637 1.385321  1.3502057  0.6128    0.2256 0.430..
#'# 3 2G_n_het       5  Large  2.222222e-02 0.10446123 0.754386  0.8626523  0.6748    0.3496 0.693..
#'as.data.frame(calculate2GBias(mean=0,sd=1,diff=c(0.283,0.707104,1.131374),
#'  Expected.StdMD=c(0.157,0.392,0.628),Expected.PHat=c(0.556,0.636,0.705),N=10, reps=20,
#'  type="lap",seed=1423,StdAdj=0.5 ))
#'  #Parameter reps changed due to NOTE "Examples with CPU (user + system) or elapsed time > 5s"
#'  #Results for reps=100:
#'#      Design GrpSize   Diff      NPBias    StdMDBias   NPMdMRE StdMDMdMRE ObsPHat ObsCliffd  Ob..
#'#1 2G_lap_het      10  Small -0.11071429 -0.080855612 1.8928571  2.1256888  0.5498    0.0996 0.1..
#'#2 2G_lap_het      10 Medium -0.07426471  0.003940804 0.6323529  0.8170856  0.6259    0.2518 0.3..
#'#3 2G_lap_het      10  Large -0.05756098  0.023696619 0.4146341  0.5447941  0.6932    0.3864 0.6..

calculate2GBias = function(mean = 0,
                           sd = 1,
                           N,
                           reps,
                           diff = c(0.2, 0.5, 0.8),
                           Expected.StdMD = c(0.2, 0.5, 0.8),
                           Expected.PHat = c(0.556, 0.638, 0.714),
                           type = "n",
                           seed = 223,
                           StdAdj = 0) {
  MdMRETable = NULL

  NumESizes = length(diff)


  Design = "2G"
  Design = paste(Design, type, sep = "_")

  DiffGen = c("Small", "Medium", "Large")

  if (StdAdj > 0) {
    Design = paste(Design, "het", sep = "_")
  }

  for (i in 1:NumESizes) {
    Out1 = reproducer::RandomExperimentSimulations(mean,
                                                   sd,
                                                   diff[i],
                                                   N,
                                                   reps,
                                                   type,
                                                   seed,
                                                   StdAdj,
                                                   returnData = TRUE)

    CentralPhat = Out1$PHat - 0.5

    ExpectedCentralPhat = Expected.PHat[i] - 0.5

    Expected.Cliffd = (Expected.PHat[i] - 0.5) * 2

    NPBias = (mean(CentralPhat) - ExpectedCentralPhat) / ExpectedCentralPhat

    StdMDBias = (mean(Out1$StdES) - Expected.StdMD[i]) / Expected.StdMD[i]

    NPMdMRE = median(abs((Out1$Cliffd - Expected.Cliffd) / Expected.Cliffd))

    StdMDMdMRE = median(abs((Out1$StdES - Expected.StdMD[i]) / Expected.StdMD[i]))


    MdMRETable = tibble::tibble(dplyr::bind_rows(
      MdMRETable,
      dplyr::bind_cols(
        Design = Design,
        GrpSize = as.character(N),
        Diff = DiffGen[i],
        NPBias = NPBias,
        StdMDBias =  StdMDBias,
        NPMdMRE = NPMdMRE,
        StdMDMdMRE = StdMDMdMRE,
        ObsPHat = mean(Out1$PHat),
        ObsCliffd = mean(Out1$Cliffd),
        ObsStdES = mean(Out1$StdES),
        PHatPower = mean(Out1$PHatSig),
        CliffdPower = mean(Out1$CliffdSig),
        StdESPower = mean(Out1$ESSig)
      )
    ))

  }
  return(MdMRETable)

}

###########################################################################################################################################
#' @title calculate4GBias
#' @description The function simulates four-group experiments and estimates of the power, individual estimate error and small sample bias obtained from a set of simulated experiments. The function produces three set of simulations obtained using three different values of the mean difference between the treatment and control groups as specified by the parameter "diff". The power is estimated as the percentage of simulated experiments for which the mean of the experiment was significantly different from zero using one-sided tests. The experiment data may be one of four different type: Normal, Log-normal, Gamma or Laplace. The output is a table of values identifying the observed values of three effect sizes: Cliff's d, PHat and StdMD, their relted estimate error, small sample bias and power for each set of simulated experiments. This function supports the production of the values reported in data tables in the paper "Recommendations for Analyzing Small Sample Size" and its Supplementary Material.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export calculate4GBias
#' @param mean This is the mean value of the control group(s) used in the simulations of each experiment for simulations of a specified mean difference (default 0).
#' @param sd This is the standard deviation  value of the control group(s) and treatment group(s) used in the simulations of each experiment of each family for simulations of a specified sample size (default 1).
#' @param N This specifies the sample size per group that will be used in each set of simulations.
#' @param reps The number of families simulated for each sample size.
#' @param diff This specifies the difference between the control and treatment that will be used in each set of simulations. It must always have three values representing small, medium and large mean differences (default c(0.2, 0.5, 0.8)).
#' @param Expected.StdMD This defines the expected value of the overall average StdMD for each mean difference (default c(0.2, 0.5, 0.8)).
#' @param Expected.PHat This defines the expected population value of the overall average Phat for each mean difference (default c(0.556,0.638,0.714)).
#' @param type This specifies the distribution of the data samples that will be simulated. Options ae "n" for Normal, "l", for Log-normal,'g" for Gamma, "lap" for LaPlace (default "n").
#' @param seed A seed for the simulations (default 123).
#' @param StdAdj Used to introduce variance heterogeneity for Laplace and Normal samples (default 0).
#' @param Blockmean Specifies he value of the block effect (default 0).
#' @return Design. Specifies the type of experiment 2G, the sample distribution (n,l,g,lap), and whether variance heterogeneity was added (het)
#' @return BEIncluded. Specifies whether or not a block effect was introduced.
#' @return GrpSize. Specifies the size of each group in the individual experiments.
#' @return Diff. The size of the difference between the control and treatment converted to an ordinal scale (Small, Medium, Large)
#' @return NPBias The relative difference between the average of the observed values of either Cliff's d or centralised PHat and the population value
#' @return StdMDBias. The relative difference between the average of the observed values of StdMDBias and the theoretical value
#' @return NPMdMRE The median of the absolute relative difference between the observed values of either Cliff's d or centralised PHat and the theoretical value for each experiment.
#' @return StdMDMdMRE The median of the relative difference between the observed values of StdMD and the population value for each experiment.
#' @return ObsPHat. The average Phat value found for each simulation.
#' @return ObsCliffd. The  average Cliffd value found for each simulation.
#' @return ObsStdES. The average of StdMD calculated for each simulation.
#' @return PHatPower. The proportion of the simulations, for a given mean difference, for which the  Phat estimate was significantly different from zero at the 0.05 alpha level based on one-sided tests.
#' @return CliffdPower. The proportion of the simulations, for a given mean difference, for which the Cliff's d estimate was significantly different from zero at the 0.05 alpha level based on one-sided tests.
#' @return StdMDPower. The proportion of the simulations, for a given mean difference, for which the StdMD estimate was significantly different from zero at the 0.05 alpha level based on one-sided tests.
#' @examples
#'#as.data.frame(calculate4GBias(mean=0,sd=1,diff=c(0.266,0.72375,1.43633),
#'#  Expected.StdMD=c(0.2,0.5,0.8),Expected.PHat=c(0.575,0.696,0.845),N=10,reps=200,type="l",
#'#  seed=17+1823,StdAdj=0,Blockmean=0))
#'#  Design BEIncluded GrpSize   Diff      NPBias StdMDBias   NPMdMRE StdMDMdMRE  ObsPHat ObsCliffd.
#'#  1 4G_l         No      10  Small -0.05933333 0.1247408 0.8666667  1.2047848 0.570550   0.1411..
#'#  2 4G_l         No      10 Medium -0.01760204 0.1565643 0.3112245  0.4426859 0.692550   0.3851..
#'#  3 4G_l         No      10  Large -0.00326087 0.2273638 0.1594203  0.2924361 0.843875   0.6877..
#' as.data.frame(calculate4GBias(mean=1,sd=3,diff=c(0.1225,0.3415,0.6224),
#'  Expected.StdMD=c(-0.208,-0.52,-0.833),Expected.PHat=c(0.444,0.360,0.277),N=20,reps=30,type="g",
#'  seed=17+977,StdAdj=0 ,Blockmean=0.5))
#'# Results for reps=200:
#'#  Design BEIncluded GrpSize   Diff     NPBias  StdMDBias   NPMdMRE StdMDMdMRE   ObsPHat  ObsCli..
#'#1   4G_g        Yes      20  Small 0.04274554 0.02242895 0.8370536  0.7960052 0.4416062 -0.1167..
#'#2   4G_g        Yes      20 Medium 0.01959821 0.01585829 0.3348214  0.3210435 0.3572562 -0.2854..
#'#3   4G_g        Yes      20  Large 0.01303251 0.01515967 0.1905830  0.1871956 0.2740938 -0.4518..


calculate4GBias = function(mean = 0,
                           sd = 1,
                           N,
                           reps,
                           diff = c(0.2, 0.5, 0.8),
                           Expected.StdMD = c(0.2, 0.5, 0.8),
                           Expected.PHat = c(0.556, 0.638, 0.714),
                           type = "n",
                           seed = 223,
                           StdAdj = 0,
                           Blockmean = 0) {
  MdMRETable = NULL

  NumESizes = length(diff)


  Design = "4G"
  Design = paste(Design, type, sep = "_")

  DiffGen = c("Small", "Medium", "Large")

  if (StdAdj > 0) {
    Design = paste(Design, "het", sep = "_")
  }

  BEIncluded = "No"
  if (Blockmean > 0)
    BEIncluded = "Yes"

  for (i in 1:NumESizes) {
    Out1 = reproducer::RandomizedBlocksExperimentSimulations(
      mean = mean,
      sd = sd,
      diff = diff[i],
      N = N,
      reps = reps,
      type = type,
      seed = seed,
      StdAdj = StdAdj,
      Blockmean = 0,
      returnData = TRUE
    )

    CentralPhat = Out1$PHat - 0.5

    ExpectedCentralPhat = Expected.PHat[i] - 0.5

    Expected.Cliffd = (Expected.PHat[i] - 0.5) * 2

    NPBias = (mean(CentralPhat) - ExpectedCentralPhat) / ExpectedCentralPhat

    StdMDBias = (mean(Out1$StdES) - Expected.StdMD[i]) / Expected.StdMD[i]

    NPMdMRE = median(abs((Out1$Cliffd - Expected.Cliffd) / Expected.Cliffd))

    StdMDMdMRE = median(abs((Out1$StdES - Expected.StdMD[i]) / Expected.StdMD[i]))


    MdMRETable = tibble::tibble(dplyr::bind_rows(
      MdMRETable,
      dplyr::bind_cols(
        Design = Design,
        BEIncluded = BEIncluded,
        GrpSize = as.character(N),
        Diff = DiffGen[i],
        NPBias = NPBias,
        StdMDBias =  StdMDBias,
        NPMdMRE = NPMdMRE,
        StdMDMdMRE = StdMDMdMRE,
        ObsPHat = mean(Out1$PHat),
        ObsCliffd = mean(Out1$Cliffd),
        ObsStdES = mean(Out1$StdES),
        PHatPower = mean(Out1$PHatSig),
        CliffdPower = mean(Out1$CliffdSig),
        StdESPower = mean(Out1$ESSig)
      )
    ))

  }
  return(MdMRETable)

}

#################################################################################################################################
#' @title calculate4GType1Error
#' @description The function simulates multiple four-group experiments and estimates the Type1 Error rate obtained from the set of simulated experiments. The Type1 Error is estimated as the percentage of experiments for which the mean the experiment was significantly different from zero at the 0.05 significance level using two-sided tests. The experiment data may be one of four different type: Normal, Log-normal, Gamma or Laplace. The output is a set of values identifying three observed effect size estimates (Cliff's d, PHat and StdMD) and their related type 1 error rates. This function supports the production of the values reported in data tables in the paper "Recommendations for Analyzing Small Sample Size Software Engineering Experiments" and its Supplementary Material.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export  calculate4GType1Error
#' @param mean This is the mean value of the control and treatment group(s) used in the simulations (default 0).
#' @param sd This is the standard deviation  value of the control group(s) and treatment group(s) used in the simulations (default 1).
#' @param N This specifies the sample size per group that will be used in each simulation (default 5).
#' @param reps The number of experiments to simulated.
#' @param type This specifies the distribution of the data samples that will be simulated. Options are "n" for Normal, "l", for Log-normal,'g" for Gamma, "lap" for LaPlace (default "n").
#' @param seed A seed for the simulations (default 123).
#' @param StdAdj Used to introduce variance heterogeneity for Laplace and Normal samples (default 0).
#' @param Blockmean Used to specify the block effect (default 0).
#' @return Design. Specifies the type of experiment 2G or 4G, the sample distribution (n,l,g,lap), and whether variance heterogeneity was added (het)
#' @return GrpSize. Specifies the size of each group in the simulations.
#' @return BEIncluded. Specifies whether or not a block effect was introduced.
#' @return ObsPHat. The average of the average Phat values found in the set of simulations.
#' @return ObsCliffd. The average of the average Cliffd values found  in the set of simulations.
#' @return ObsStdES. The average of StdMD values found in the set of simulations.
#' @return PHatType1ER. The percentage of the simulations for which the Phat estimate was significantly different from zero at the 0.05 alpha level.
#' @return CliffdType1ER. The percentage of the simulations for which the overall Cliff's d estimate was significantly different from zero at the 0.05 alpha level.
#' @return StdMDType1ER. The percentage of the simulations for which the overall StdMD estimate was significantly different from zero at the 0.05 significance level.
#' @examples
#'as.data.frame(calculate4GType1Error(mean=0,sd=1,N=40,reps=100,type="n",seed=17+1056,StdAdj = 0.5,
#'  Blockmean=0.5))
#'  # Results for reps=300
#' #    Design GrpSize BEIncluded   ObsPHat   ObsCliffd   ObsStdES PHatType1ER CliffdType1ER StdES..
#' #1 4G_n_het      40        Yes 0.5034729 0.006945833 0.01316457        0.03    0.02333333 0.046..
#'
#'#as.data.frame(calculate4GType1Error(mean=0,sd=1,N=40,reps=300,type="lap",seed=17+2056,
#'#  StdAdj = 0.5,Blockmean=0.5))
#' #      Design GrpSize BEIncluded   ObsPHat    ObsCliffd  ObsStdES PHatType1ER CliffdType1ER Std..
#' #1 4G_lap_het      40        Yes 0.4992708 -0.001458333 0.0014446  0.04333333          0.04  0.06
calculate4GType1Error = function(mean = 0,
                                 sd = 1,
                                 N = 10,
                                 reps = 10,
                                 type = "n",
                                 seed = 123,
                                 StdAdj = 0,
                                 Blockmean = 0) {
  Type1ErrorTable = NULL

  Out1 = reproducer::RandomizedBlocksExperimentSimulations(
    mean = mean,
    sd = sd,
    diff = 0,
    N = N,
    reps = reps,
    type = type,
    seed = seed,
    StdAdj = StdAdj,
    Blockmean = Blockmean,
    returnData = TRUE
  )

  Design = "4G"
  Design = paste(Design, type, sep = "_")

  if (StdAdj > 0) {
    Design = paste(Design, "het", sep = "_")
  }

  BEIncluded = "No"
  if (Blockmean > 0) {
    BEIncluded = "Yes"
  }


  Type1ErrorTable = tibble::tibble(
    dplyr::bind_cols(
      Design = Design,
      GrpSize = as.character(N),
      BEIncluded = BEIncluded,
      ObsPHat = mean(Out1$PHat),
      ObsCliffd = mean(Out1$Cliffd),
      ObsStdES = mean(Out1$StdES),
      PHatType1ER = mean(Out1$PHatSig),
      CliffdType1ER = mean(Out1$CliffdSig),
      StdESType1ER = mean(Out1$ESSig)
    )
  )


  return(Type1ErrorTable)
}


#
# The functions for meta-analysis of Cliff's d and PHat
#

#' @title calcCliffdTestStatistics
#' @description This function is a helper function for meta-analysis of experiments using Cliff's d as an effect size. It returns the 100*(1-alpha/2)% confidence intervals based on the normal distribution probability values, the value of the t-test, the probability asssociated with the null hypothesis and the significance of the test. The pvalue and the significance vary according to the value of the alternative parameter and whether or not degrees of freedom are specifed.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export calcCliffdTestStatistics
#' @param d.value The overall estimate of Cliff's d from a group of effect sizes to be meta-analysed
#' @param d.variance The estimate of the variance of the overall estimate of Cliff's d
#' @param d.df The total degrees of freedom for the set of effect sizes. If d.df>0, the pvalues and significance test use the t-distribution probability values. If d.df=0 (default) the pvalues and significance test use the normal distribution probability values. The confidence intervals are always based on the normal probability values.
#' @param alpha The significance level used to control the significance tests and calculation of confidence limits (default 0.05).
#' @param alternative Specifies the type of significance test and can take the values "two.sided", "less" or "greater" (default "two.sided").
#' @return d.tvalue The value of the t-statistic
#' @return d.pvalue The p-value of the t-test if the parameter d.df>0, or the normal probability value if d.df=0
#' @return d.ci.lower The lower 100*(1-alpha/2)% confidence interval based on the normal probability value
#' @return d.ci.upper The upper 100*(1-alpha/2)% confidence interval based on the normal probability value
#' @return d.sig The significance of the statistical test of the d.tvalue return value at the alpha level for one sided tests and aplha/2 for two sided tests as specified by the input parameter alternative
#' @examples
#' aveCliffd=mean(c(0.84,0.2,-0.04,0.44,0.76))
#' aveCliffdvar=sum(c(0.04,0.18,0.21,0.15,0.06))/25
#' df=45
#' calcCliffdTestStatistics(d.value=aveCliffd,d.variance=aveCliffdvar,d.df=df)
#' # A tibble: 1 x 5
#'#   d.tvalue d.pvalue d.ci.lower d.ci.upper d.sig
#'#      <dbl>    <dbl>      <dbl>      <dbl> <lgl>
#' # 1     2.75  0.00855     0.0923      0.692 TRUE
calcCliffdTestStatistics = function(d.value,
                                    d.variance,
                                    d.df = 0,
                                    alpha = 0.05,
                                    alternative = "two.sided") {
  d <- d.value
  vard <- d.variance
  d.tvalue <- d / sqrt(vard)

  useTTest = d.df > 0


  if (alternative == "two.sided") {
    zv <- stats::qnorm(alpha / 2)
  }
  else {
    zv <- stats::qnorm(alpha)
  }


  d.cu <-
    (d - d ^ 3 - zv * sqrt(vard) * sqrt((1 - d ^ 2) ^ 2 + zv ^ 2 * vard)) /
    (1 - d ^ 2 + zv ^ 2 * vard)
  d.cl <-
    (d - d ^ 3 + zv * sqrt(vard) * sqrt((1 - d ^ 2) ^ 2 + zv ^ 2 * vard)) /
    (1 - d ^ 2 + zv ^ 2 * vard)
  d.sig <- d.cu < 0 | d.cl > 0

  if (useTTest) {
    d.pvalue <- 2 * (1 - stats::pt(abs(d.tvalue), d.df))
  }

  else {
    d.pvalue <- 2 * (1 - stats::pnorm(abs(d.tvalue)))
  }


  # If a one-sided test correct pvalue and significance
  if (alternative == "greater") {
    temp.cl <-
      (d - d ^ 3 + zv * sqrt(vard) * sqrt((1 - d ^ 2) ^ 2 + zv ^ 2 * vard)) /
      (1 - d ^ 2 + zv ^ 2 * vard)
    d.sig <- temp.cl > 0

    if (useTTest) {
      d.pvalue <- 1 - stats::pt(d.tvalue, d.df)
    }
    else {
      d.pvalue <- 1 - stats::pnorm(d.tvalue)
    }
  }

  if (alternative == "less") {
    temp.cu <-
      (d - d ^ 3 - zv * sqrt(vard) * sqrt((1 - d ^ 2) ^ 2 + zv ^ 2 * vard)) / (1 - d ^ 2 + zv ^ 2 * vard)

    d.sig <- temp.cu < 0

    if (useTTest) {
      d.pvalue <- stats::pt(d.tvalue, d.df)
    }
    else {
      d.pvalue <- stats::pnorm(d.tvalue)
    }
  }



  out = tibble::tibble(
    d.tvalue = d.tvalue,
    d.pvalue = d.pvalue,
    d.ci.lower = d.cl,
    d.ci.upper = d.cu,
    d.sig = d.sig
  )

  return(out)

}



#' @title metaanalyse.Cliffd
#' @description This function provides a simple meta-analysis of experiments using Cliff's d as an effect size. It returns the 100*(1-alpha/2)% confidence intervals based on the normal distribution probability values, the value of the t-test, the probability asssociated with the null hypothesis and the significance of the test. The pvalue and the significance vary according to the value of the parameter "alternative" and whether or not degrees of freedom are specifed. It also return the heterogeneity statistics Q and I-squared.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export metaanalyse.Cliffd
#' @param Cliffd A vector of one or more numerical values, identifying the effect sizes to be meta-analysed
#' @param  Cliffdvar A vector of the estimates variance of each of the effect sizes
#' @param df The total degrees of freedom for the set of effect sizes. If df>0, the pvalues and significance test use the t-distribution probability values. If df=0 (default) the pvalues and significance test use the normal distribution probability values. The confidence intervals are always based on the normal probability values, as recommended by Cliff.
#' @param alternative Specifies the type of significance test and can take the values "two.sided" (default), "less" or "greater".
#' @param alpha The significance level used to control the significance tests and calculation of confidence limits (default 0.05).
#' @return Estimate The overall estimate of Cliff's d obtained from the set of experiments
#' @return UpperCI The upper 100*(1-alpha/2)% confidence interval based on the normal probability value
#' @return LowerCI The lower 100*(1-alpha/2)% confidence interval based on the normal probability value
#' @return The variance of the Estimate
#' @return tvalue The value of the t-statistic
#' @return df The supplied degrees of freedom or NA if the input parameter df was set to zero
#' @return AltHyp Defines the alternative hypothesis used for significance testing and depends on the value of the input parameter alternative. It takes the values "Not=0", ">0", or "<0"
#' @return NullHyp Defines the null hypothesis and depends on the value of the input parameter alternative. It takes the values "~0", "<0", or ">0"
#' @return pvalue The p-value of the t-test if the parameter df>0, or the normal probability value if d=0
#' @return RejectNullHyp "Yes" or "No" depending on whether or not the null hypothesis should be rejected at the alpha/2 level for two-sided tests and alpha level for one-sided tests
#' @return The Q homogeneity statistic
#' @return The I-squared estimate of the extent of heterogeneity
#' @return ProbQHomogeneous. The probability that the set of Cliff's d values come from a set of homogeneous experiments.
#' @examples

#'Cliffd=c(0.84,0.2,-0.04,0.44,0.76)
#'CliffdvarInvalid=c(0.04,0.18,0.21,0.15)
#'Cliffdvar=c(0.04,0.18,0.21,0.15,0.06)
#'CliffdvarInvalid=c(0.04,0.18,0.21,0.15)
#'df=45

#'as.data.frame(metaanalyse.Cliffd(Cliffd=Cliffd,Cliffdvar=Cliffdvar,df=df,alternative="greater",
#'  alpha=0.05))
#' #  Estimate   UpperCI   LowerCI Variance tvalue df AltHyp NullHyp      pvalue
#' #1     0.44 0.6601381 0.1502568   0.0256   2.75 45     >0     <=0 0.004275955
#' #  RejectNullHyp    Q I.square ProbQHomogeneous
#' #1           Yes 21.5 81.39535     0.0002519835
#'as.data.frame(metaanalyse.Cliffd(Cliffd=Cliffd,Cliffdvar=Cliffdvar,df=df,alternative="less",
#'  alpha=0.05))
#' #  Estimate   UpperCI   LowerCI Variance tvalue df AltHyp NullHyp   pvalue RejectNullHyp
#' #1     0.44 0.6601381 0.1502568   0.0256   2.75 45     <0     >=0 0.995724            No
#' #     Q I.square ProbQHomogeneous
#' #1 21.5 81.39535     0.0002519835
#'as.data.frame(metaanalyse.Cliffd(Cliffd=Cliffd,Cliffdvar=Cliffdvar,df=df,alternative="two.sided",
#'  alpha=0.05))
#' #  Estimate  UpperCI    LowerCI Variance tvalue df AltHyp NullHyp      pvalue
#' #1     0.44 0.692073 0.09227496   0.0256   2.75 45  Not=0      ~0 0.008551911
#' #  RejectNullHyp    Q I.square ProbQHomogeneous
#' #1           Yes 21.5 81.39535     0.0002519835
#'as.data.frame(metaanalyse.Cliffd(Cliffd=Cliffd,Cliffdvar=Cliffdvar,df=df,alpha=0.05))
#' #  Estimate  UpperCI    LowerCI Variance tvalue df AltHyp NullHyp      pvalue
#' #1     0.44 0.692073 0.09227496   0.0256   2.75 45  Not=0      ~0 0.008551911
#' #  RejectNullHyp    Q I.square ProbQHomogeneous
#' #1           Yes 21.5 81.39535     0.0002519835

#' metaanalyse.Cliffd(Cliffd=Cliffd,Cliffdvar=Cliffdvar,df=0,alternative="two.sided",alpha=0.05)
#' #Error in testfunctionParameterChecks(alternative = alternative, alpha = alpha,  :
#' #  Invalid alternative parameter, choose one of two.sided, greater or less
#'# metaanalyse.Cliffd(Cliffd=Cliffd,Cliffdvar=CliffdvarInvalid,df=df,alternative="greater",
#'# alpha=0.05)
#' #Error in metaanalyse.Cliffd(Cliffd = Cliffd, Cliffdvar = CliffdvarInvalid,  :
#' #  Length of Cliffdvar parameter must equal the length of the Cliffd parameter

metaanalyse.Cliffd <-
  function(Cliffd,
           Cliffdvar,
           df = 0,
           alternative = "two.sided",
           alpha = 0.05) {
    NumExp <- length(Cliffd)

    if (NumExp != length(Cliffdvar)) {
      stop("Length of Cliffdvar parameter must equal the length of the Cliffd parameter")
    }


    AveCliffd <- base::mean(Cliffd)
    AveCliffdvar <- sum(Cliffdvar) / NumExp ^ 2

    testfunctionParameterChecks(
      alternative = alternative,
      alpha = alpha,
      stderr = sqrt(AveCliffdvar)
    )


    AveCliffCI <-
      calcCliffdTestStatistics(
        d.value = AveCliffd,
        d.variance = AveCliffdvar,
        d.df = df,
        alpha = alpha,
        alternative = alternative
      )

    Alternative = "Not=0"

    NullHyp = "~0"


    if (alternative == "greater") {
      Alternative = ">0"
      NullHyp = "<=0"
    }

    if (alternative == "less") {
      Alternative = "<0"
      NullHyp = ">=0"
    }


    #Heterogeneity analysis
    Q <- sum((Cliffd - AveCliffd) ^ 2) / AveCliffdvar

    q.df <- NumExp - 1

    ProbQHomogeneous = 1 - pchisq(Q, q.df)

    I.square <- 100 * (Q - (NumExp - 1)) / Q

    if (I.square < 0) {
      I.square = 0
    }

    if (df == 0) {
      df = "NA"
    }

    if (AveCliffCI$d.sig) {
      RejectNullHyp = "Yes"
    } else {
      RejectNullHyp = "No"
    }

    UpperCI = AveCliffCI$d.ci.upper
    if (UpperCI > 1) {
      UpperCI = 1
    }

    LowerCI = AveCliffCI$d.ci.lower
    if (LowerCI < -1) {
      LowerCI = -1
    }

    out <-
      tibble::tibble(
        Estimate = AveCliffd,
        UpperCI = UpperCI,
        LowerCI = LowerCI,
        Variance = AveCliffdvar,
        tvalue = AveCliffCI$d.tvalue,
        df = df,
        AltHyp = Alternative,
        NullHyp = NullHyp,
        pvalue = AveCliffCI$d.pvalue,
        RejectNullHyp = RejectNullHyp,
        Q = Q,
        I.square = I.square,
        ProbQHomogeneous = ProbQHomogeneous
      )


    return(out)

  }

#' @title PHatonesidedTestStatistics
#' @description This function is a helper function for meta-analysis of experiments using PHat as an effect size. It returns the 100*(1-alpha)% confidence intervals, the value of the t-test, the probability asssociated with the null hypothesis and the significance for a one-sided test. The direction of the one-sided test is determined by the parameter "alternative" which takes the values "greater" or "less".
#' @author Barbara Kitchenham and Lech Madeyski
#' @param effectsize The overall estimate of the centralized PHat (i.e. Phat-0.5) from a group of effect sizes to be meta-analysed
#' @param effectsize.variance The estimate of the variance of the overall estimate ofPHat
#' @param effectsize.df The total degrees of freedom for the set of effect sizes. If effectsize.df>0, the confidence intervals, pvalues and significance test use the t-distribution probability values. If effectsize.df=0 (default), the confidence intervals, the pvalues and significance test use the normal distribution probability values.
#' @param alpha The significance level (default 0.05) used to control the significance tests and calculation of confidence limits.
#' @param alternative Specifies the type of significance test and can take the values "less" or "greater" (default).
#' @return ES.test The value of the t-statistic
#' @return ES.pvalue The p-value of the two-sided t-test if the parameter d.df>0, or the normal probability value if d.df=0
#' @return ES.sig The significance of the statistical test of the d.tvalue return value at the alpha level for one sided tests and aplha/2 for two sided tests as specified by the input parameter alternative.
#' @return ES.ci.lower The lower 100*(1-alpha/2)% confidence interval of the average centralized PHat based on the t-distribution probability values if effectsize.df>0 or normal probability values if effectsize.df=0
#' @return ES.ci.upper The upper 100*(1-alpha/2)% confidence interval of the average centralized PHat based on the t-distribution probability values if effectsize.df>0 or normal probability values if effectsize.df=0
#' @examples
#'PHatES=mean(c(0.92,0.6,0.48,0.72,0.88))-0.5
#'PHatESvar=sum(c(0.01,0.04,0.05,0.04,0.01))/25
#'PHatdf=sum(c(6.63,6.63,5.08,5.61,8))

#'#PHatonesidedTestStatistics(effectsize=PHatES,effectsize.variance=PHatESvar,effectsize.df=PHatdf)
#'# A tibble: 1 x 5
#'#  ES.test ES.pvalue ES.sig ES.ci.lower ES.ci.upper
#'#    <dbl>     <dbl> <lgl>        <dbl>       <dbl>
#'#1    2.84   0.00389 TRUE        0.0888       0.351
#'#PHatonesidedTestStatistics(effectsize=PHatES,effectsize.variance=PHatESvar,effectsize.df=0,
#'# alternative="less")
#'# A tibble: 1 x 5
#'#  ES.test ES.pvalue ES.sig ES.ci.lower ES.ci.upper
#'#    <dbl>     <dbl> <lgl>        <dbl>       <dbl>
#'#1    2.84     0.998 FALSE       0.0926       0.347

PHatonesidedTestStatistics = function(effectsize,
                                      effectsize.variance,
                                      effectsize.df = 0,
                                      alpha = 0.05,
                                      alternative = "greater")
{
  ES.se <- sqrt(effectsize.variance)
  ES.test <- effectsize / ES.se
  useTTest <- (effectsize.df > 0)


  if (useTTest) {
    vv <- stats::qt(alpha, effectsize.df)
  }
  else {
    vv <- stats::qnorm(alpha)
  }

  ES.ci.lower <- effectsize + vv * ES.se
  ES.ci.upper <- effectsize - vv * ES.se

  if (alternative == "greater") {
    ES.sig <- ES.ci.lower > 0
    if (useTTest) {
      ES.pvalue <- (1 - stats::pt(ES.test, effectsize.df))
    }
    else {
      ES.pvalue <- (1 - stats::pnorm(ES.test))
    }
  }
  else
  {
    ES.sig <- ES.ci.upper < 0

    if (useTTest) {
      ES.pvalue <- (stats::pt(ES.test, effectsize.df))
    }
    else {
      ES.pvalue <- (stats::pnorm(ES.test))
    }

  }


  out <- tibble::tibble(
    ES.test = ES.test,
    ES.pvalue = ES.pvalue,
    ES.sig = ES.sig,
    ES.ci.lower = ES.ci.lower,
    ES.ci.upper = ES.ci.upper
  )

  return(out)

}


#' @title PHattwosidedTestStatistics
#' @description This function is a helper function for meta-analysis of experiments using PHat as an effect size. It returns the 100*(1-alpha/2)% confidence intervals, the value of the t-test, the probability asssociated with the null hypothesis and the significance for a two-sided test.
#' @author Barbara Kitchenham and Lech Madeyski
#' @param effectsize The overall estimate of the centralized PHat (ie.Phat-0.5) from a group of effect sizes to be meta-analysed
#' @param effectsize.variance The estimate of the variance of the overall estimate ofPHat
#' @param effectsize.df The total degrees of freedom for the set of effect sizes. If effectsize.df>0, the confidence intervals, pvalues and significance test use the t-distribution probability values. If effectsize.df=0 (default), the confidence intervals, the pvalues and significance test use the normal distribution probability values.
#' @param alpha The significance level (default 0.05) used to control the significance tests and calculation of confidence limits.
#' @return ES.test The value of the t-statistic
#' @return ES.pvalue The p-value of the two-sided t-test if the parameter d.df>0, or the normal probability value if d.df=0
#' @return ES.sig The significance of the statistical test of the d.tvalue return value at the alpha level for one sided tests and aplha/2 for two sided tests as specified by the input parameter alternative.
#' @return ES.ci.lower The lower 100*(1-alpha/2)% confidence interval of the average centralized PHat based on the t-distribution probability values if effectsize.df>0 or normal probability values if effectsize.df=0
#' @return ES.ci.upper The upper 100*(1-alpha/2)% confidence interval of the average centralized PHat based on the t-distribution probability values if effectsize.df>0 or normal probability values if effectsize.df=0
#' @examples
#'PHatES=mean(c(0.92,0.6,0.48,0.72,0.88))-0.5
#'PHatESvar=sum(c(0.01,0.04,0.05,0.04,0.01))/25
#'PHatdf=sum(c(6.63,6.63,5.08,5.61,8))

#'#PHattwosidedTestStatistics(effectsize=PHatES,effectsize.variance=PHatESvar)
#' # A tibble: 1 x 5
#' # ES.test ES.pvalue ES.sig ES.ci.lower ES.ci.upper
#' #     <dbl>     <dbl> <lgl>        <dbl>       <dbl>
#' # 1    2.84   0.00451 TRUE        0.0682       0.372
#' # PHattwosidedTestStatistics(effectsize=PHatES,effectsize.variance=PHatESvar,effectsize.df=PHatdf)
#' #  A tibble: 1 x 5
#' #   ES.test ES.pvalue ES.sig ES.ci.lower ES.ci.upper
#' #     <dbl>     <dbl> <lgl>        <dbl>       <dbl>
#' # 1    2.84   0.00778 TRUE        0.0622       0.378

PHattwosidedTestStatistics = function(effectsize,
                                      effectsize.variance,
                                      effectsize.df = 0,
                                      alpha = 0.05) {
  ES.se = sqrt(effectsize.variance)
  ES.test = effectsize / ES.se
  useTTest = (effectsize.df > 0)



  if (useTTest) {
    vv = stats::qt(alpha / 2, effectsize.df)
    ES.pvalue = 2 * (1 - stats::pt(abs(ES.test), effectsize.df))
  }
  else {
    vv = stats::qnorm(alpha / 2)
    ES.pvalue = 2 * (1 - stats::pnorm(abs(ES.test)))
  }
  ES.ci.lower = effectsize + vv * ES.se
  ES.ci.upper = effectsize - vv * ES.se
  ES.sig = (ES.ci.upper < 0 | ES.ci.lower > 0)



  out = tibble::tibble(
    ES.test = ES.test,
    ES.pvalue = ES.pvalue,
    ES.sig = ES.sig,
    ES.ci.lower = ES.ci.lower,
    ES.ci.upper = ES.ci.upper
  )

  return(out)

}



#' @title calcPHatMATestStatistics
#' @description This function is a helper function for meta-analysis of experiments using PHat as an effect size. It returns the 100*(1-alpha/2)% confidence intervals,the value of the t-test, the probability associated with the null hypothesis and the significance of the test. The confidence intervals, pvalue and the significance of the statistical test vary according to the value of the "alternative" parameter and whether or not degrees of freedom are specifed.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export calcPHatMATestStatistics
#' @param effectsize The overall estimate of the centralized PHat (ie.Phat-0.5) from a group of effect sizes to be meta-analysed
#' @param effectsize.variance The estimate of the variance of the overall estimate ofPHat
#' @param effectsize.df The total degrees of freedom for the set of effect sizes. If effectsize.df>0, the confidence intervals, pvalues and significance test use the t-distribution probability values. If effectsize.df=0 (default), the confidence intervals, the pvalues and significance test use the normal distribution probability values.
#' @param alpha The significance level used to control the significance tests and calculation of confidence limits (default 0.05).
#' @param alternative Specifies the type of significance test and can take the values "two.sided" (default), "less" or "greater"
#' @return ES.test The value of the t-statistic
#' @return ES.pvalue The p-value of the t-test if the parameter d.df>0, or the normal probability value if d.df=0
#' @return ES.sig The significance of the statistical test of the d.tvalue return value at the alpha level for one sided tests and aplha/2 for two sided tests as specified by the input parameter alternative.
#' @return ES.ci.lower The lower 100*(1-alpha/2)% confidence interval of the average centralized PHat based on the t-distribution probability values if effectsize.df>0 or normal probability values if effectsize.df=0
#' @return ES.ci.upper The upper 100*(1-alpha/2)% confidence interval of the average centralized PHat based on the t-distribution probability values if effectsize.df>0 or normal probability values if effectsize.df=0
#' @examples
#'avePHat=mean(c(0.92,0.6,0.48,0.72,0.88))
#'avePHatvar=sum(c(0.01,0.04,0.05,0.04,0.01))/25
#'PHatdf=sum(c(6.63,6.63,5.08,5.61,8))
#'calcPHatMATestStatistics(effectsize=avePHat-0.5,effectsize.variance=avePHatvar,effectsize.df=PHatdf)
#'# A tibble: 1 x 5
#'#   ES.test ES.pvalue ES.sig ES.ci.lower ES.ci.upper
#'#     <dbl>     <dbl> <lgl>        <dbl>       <dbl>
#'# 1    2.84   0.00778 TRUE        0.0622       0.378
calcPHatMATestStatistics = function(effectsize,
                                    effectsize.variance,
                                    effectsize.df = 0,
                                    alpha = 0.05,
                                    alternative = "two.sided") {
  ES.se <- sqrt(effectsize.variance)
  ES.test <- effectsize / ES.se
  useTTest <- (effectsize.df > 0)


  testfunctionParameterChecks(alternative = alternative,
                              alpha = alpha,
                              stderr = ES.se)


  CIValues <-
    PHattwosidedTestStatistics(
      effectsize = effectsize,
      effectsize.variance = effectsize.variance,
      effectsize.df = effectsize.df
    )


  if (alternative == "two.sided") {
    #Always use two-sided confidence intervals
    testStats <- CIValues
  }
  else {
    # For one sided tests adjust the p-value and significance assessment
    testStats <-
      PHatonesidedTestStatistics(
        effectsize = effectsize,
        effectsize.variance = effectsize.variance,
        effectsize.df = effectsize.df,
        alternative = alternative
      )
  }

  ES.ci.lower = CIValues$ES.ci.lower
  ES.ci.upper = CIValues$ES.ci.upper

  if (ES.ci.lower < -0.5) {
    ES.ci.lower = -0.5
  }
  if (ES.ci.upper > 0.5) {
    ES.ci.upper = 0.5
  }


  out = tibble::tibble(
    ES.test = ES.test,
    ES.pvalue = testStats$ES.pvalue,
    ES.sig = testStats$ES.sig,
    ES.ci.lower ,
    ES.ci.upper
  )

  return(out)

}


#' @title metaanalyse.PHat
#' @description This function performs a meta-analysis of experiments using PHat as an effect size. It returns the 100*(1-alpha/2)% confidence interval limits, the value of the t-test, the probability asssociated with the null hypothesis and the significance of the test. The confidence interval limits, the pvalue and the significance of the statistical test vary according to the value of the alternative parameter and whether or not degrees of freedom are specifed. It also reurns the Q and I.squared statistics to assess the heterogeneity among the experiments.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export metaanalyse.PHat
#' @param PHat The estimates of PHat obtained from a group of experiments to be meta-analysed
#' @param PHatvar The estimate of the variance of each PHat estimate
#' @param DFUnknown If DFUnknown=FALSE the degrees of freedom for each experiment is known, and the df parameter must be a vector specifying the effect size of each experiment, otherwise the df parameter is ignored.
#' @param df If DFUnknown is TRUE, this parameter is a vector of numerical values specifying the degrees of freedom for each experiment, and the confidence intervals, pvalues and significance test use the t-distribution probability values. If the parameter DFUNknown is FALSE, the confidence intervals, pvalues and significance test use the normal distribution probability values.
#' @param alternative Specifies the type of significance test and can take the values "two.sided" (default), "less" or "greater".
#' @param alpha The significance level (default 0.05) used to control the significance tests and calculation of confidence limits.
#' @return Estimate. The simple average of the PHat values recommended by Kromrey as the best estimator for meta-analysis.
#' @return UpperCI The upper 100*(1-alpha/2)% confidence interval of the average PHat based on the t-distribution probability values if DFUnknown is FALSE, or normal probability values if DFUnknown is TRUE.
#' @return LowerCI The lower 100*(1-alpha/2)% confidence interval of the average PHat based on the t-distribution probability values ifDFUnknown is FALSE, or normal probability values if DFUnknown is TRUE.
#' @return Variance The variance of the Estimate output
#' @return tvalue The value of the t-statistic
#' @return df Either NA if the parameter DFUnknown is TRUE, or sum of the degrees of freedom for each experiment.
#' @return AltHyp Defines the alternative hypothesis used for significance testing and depends on the value of the input parameter alternative. It takes the values "Not=0.5", ">0.5", or "<0.5".
#' @return NullHyp Defines the null hypothesis and depends on the value of the input parameter alternative. It takes the values "~0.5", "<0.5", or ">0.5".
#' @return pvalue The p-value of the t-test if the parameter DFUnknown is FALSE, otherwise the normal probability value.
#' @return RejectNullHyp "Yes" or "No" depending on whether or not the null hypothesis should be rejected at the alpha/2 level for two-sided tests and alpha level for one-sided tests
#' @return The I-squared estimate of the extent of heterogeneity
#' @return The Q homogeneity statistic
#' @return ProbQHomogeneous. The probability that the set of Phat values come from a set of homogeneous experiments.
#' @examples
#'PHat=c(0.92,0.6,0.48,0.72,0.88)
#'PHatvar=c(0.01,0.04,0.05,0.04,0.01)
#'PHatdf=c(6.63,6.63,5.08,5.61,8)
#'PHatInvalid=c(0.92,0.6,0.48,0.72)

#' as.data.frame(metaanalyse.PHat(PHat=PHat,PHatvar=PHatvar,DFUnknown=FALSE,df=PHatdf,
#'  alternative="greater",alpha=0.05))
#'#  Estimate   UpperCI   LowerCI Variance   tvalue    df AltHyp NullHyp      pvalue RejectNullHyp..
#'#  1     0.72 0.8777899 0.5622101    0.006 2.840188 31.95   >0.5   <=0.5 0.003890609         Yes..
#' as.data.frame(metaanalyse.PHat(PHat=PHat,PHatvar=PHatvar,DFUnknown=TRUE,df=PHatdf,
#'   alternative="greater",alpha=0.05))
#'#  Estimate   UpperCI   LowerCI Variance   tvalue df AltHyp NullHyp      pvalue RejectNullHyp..
#'# 1     0.72 0.8718182 0.5681818    0.006 2.840188 NA   >0.5   <=0.5 0.002254349           Yes..
#' as.data.frame(metaanalyse.PHat(PHat=PHat,PHatvar=PHatvar,DFUnknown=FALSE,df=PHatdf,
#'  alternative="two.sided",alpha=0.05))
#' #  Estimate   UpperCI   LowerCI Variance   tvalue    df  AltHyp NullHyp      pvalue RejectNullH..
#' #1     0.72 0.8777899 0.5622101    0.006 2.840188 31.95 Not=0.5    ~0.5 0.007781218         Yes..
#' as.data.frame(metaanalyse.PHat(PHat=PHat,PHatvar=PHatvar,DFUnknown=TRUE,df=PHatInvalid,
#'  alpha=0.05))
#' # Estimate   UpperCI   LowerCI Variance   tvalue df  AltHyp NullHyp      pvalue RejectNullHyp I..
#' #1     0.72 0.8718182 0.5681818    0.006 2.840188 NA Not=0.5    ~0.5 0.004508698         Yes 82..

metaanalyse.PHat <-
  function(PHat,
           PHatvar,
           DFUnknown,
           df,
           alternative = "two.sided",
           alpha = 0.05) {
    NumExp <- length(PHat)

    if (NumExp != length(PHatvar)) {
      stop("The length of the PHatvar parameter must equal the length  of the PHat parameter")
    }

    if (DFUnknown == FALSE) {
      if (NumExp != length(df)) {
        stop("The length of the df parameter must equal the length  of the PHat parameter")
      }
      else {
        effectsize.df <- sum(df)
      }

    }
    else {
      effectsize.df <- 0
    }



    CentralPHat <- PHat - 0.5

    AvePHat <- base::mean(CentralPHat)
    AvePHatvar <- sum(PHatvar) / NumExp ^ 2
    AvePHatCI <-
      calcPHatMATestStatistics(
        effectsize <- AvePHat,
        effectsize.variance <- AvePHatvar,
        effectsize.df <- effectsize.df,
        alpha <- alpha,
        alternative <- alternative
      )

    Alternative = "Not=0.5"

    NullHyp = "~0.5"

    if (alternative == "greater") {
      Alternative = ">0.5"
      NullHyp = "<=0.5"
    }

    if (alternative == "less") {
      Alternative = "<0.5"
      NullHyp = ">=0.5"
    }

    if (DFUnknown) {
      effectsize.df = "NA"
    }

    if (AvePHatCI$ES.sig) {
      RejectNullHyp = "Yes"
    } else {
      RejectNullHyp = "No"
    }


    #Heterogeneity analysis
    Q <- sum((CentralPHat - AvePHat) ^ 2) / AvePHatvar

    q.df <- NumExp - 1

    ProbQHomogeneous <- 1 - pchisq(Q, q.df)


    I.square <- 100 * (Q - (NumExp - 1)) / Q

    if (I.square < 0) {
      I.square <- 0
    }

    out <-
      tibble::tibble(
        Estimate = AvePHat + 0.5,
        UpperCI = 0.5 + AvePHatCI$ES.ci.upper,
        LowerCI = 0.5 + AvePHatCI$ES.ci.lower,
        Variance = AvePHatvar,
        tvalue = AvePHatCI$ES.test,
        df = effectsize.df,
        AltHyp = Alternative,
        NullHyp = NullHyp,
        pvalue = AvePHatCI$ES.pvalue,
        RejectNullHyp = RejectNullHyp,
        I.square = I.square,
        Q = Q,
        ProbQHomogeneous = ProbQHomogeneous
      )


    return(out)

  }
