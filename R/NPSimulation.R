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
#' df=c(5,10,17)
#' adjexact=calculateSmallSampleSizeAdjustment(df)
#' # > adjexact
#' # [1] 0.8407487 0.9227456 0.9551115
#' # Hedges and Olkin values 0.8408, 0.9228,0.9551
#' adjapprox=calculateSmallSampleSizeAdjustment(df,FALSE)
#' # > adjapprox
#' # [1] 0.8421053 0.9230769 0.9552239
#' # Another example:
#' df=c(10,25,50)
#' calculateSmallSampleSizeAdjustment(df,exact=TRUE)
#' # [1] 0.9227456 0.9696456 0.9849119
#' calculateSmallSampleSizeAdjustment(df,exact=FALSE)
#' # [1] 0.9230769 0.9696970 0.9849246
calculateSmallSampleSizeAdjustment = function(df, exact = TRUE) {
  exactvec = c(rep(exact, length(df)))
  # If exact is TRUE but the df is too large gamma cannot be calculated and the approximate value is used
  c = ifelse(exactvec &
               df < 340,
             sqrt(2 / df) * gamma(df / 2) / gamma((df - 1) / 2) ,
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
#' d=0.5
#' varStandardizedEffectSize(d,sqrt(1/30),29,returnVarg=FALSE)
#' # [1] 0.2007699
#' varStandardizedEffectSize(d,sqrt(1/30),29,returnVarg=TRUE)
#' # [1] 0.1904167
varStandardizedEffectSize = function(d, A, f, returnVarg = TRUE) {
  c = reproducer::calculateSmallSampleSizeAdjustment(f)
  g = d * c # g is a better estimate of the population standardized effect size delta than d
  var = (f / (f - 2)) * (A + g ^ 2) - g ^ 2 / c ^ 2 # best estimate of the variance of d
  if (returnVarg)
    var = c ^ 2 * var # best estimate of the variance of g
  return(var)
}



#' @title RandomizedBlocksAnalysis
#' @description The function performs a heteroscedastic test of a two treatment by J blocks randomized blocks effect size. The data are assumed to be stored in $x$ in list mode, a matrix or a data frame. If in list mode, length(x) is assumed to correspond to the total number of groups. All groups are assumed to be independent. Missing values are automatically removed.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export RandomizedBlocksAnalysis
#' @param x the structure holding the data. In list format, for a 2 treatment by J block randomized blocks experiments, there are  2J list elements each one specifying the outcome for a specific block and a specific treatment.
#' @param con is a 2J list containing the contrast coefficients that are used to calculate the mean effet size.
#' @param alpha default to 0.05 is the Type 1 error level used for the test of significance
#' @return The t-test and its associated metrics (i.e., critical value stansard error and degrees of freedom) and the estimate of the contrast with its upper and lower confidence interval bounds and p-value.
#' @examples
#' set.seed(123)
#' x=list()
#' x[[1]]=stats::rnorm(10,0,1)
#' x[[2]]=stats::rnorm(10,0.8,1)
#' x[[3]]=stats::rnorm(10,0.5,1)
#' x[[4]]=stats::rnorm(10,1.3,1)
#' vec=c(-1,1,-1,1)/2
#' RandomizedBlocksAnalysis(x,con=vec,alpha=0.05)
#' # $n
#' # [1] 10 10 10 10
#' # $test
#' #      test     crit        se       df
#' # [1,] 4.432644 2.038622 0.2798104 31.33793
#' # $psihat
#' #      psihat  ci.lower ci.upper      p.value
#' # [1,] 1.2403 0.6698721 1.810728 0.0001062952
#' # dat=c(x[[1]],x[[2]],x[[3]],x[[4]])
#' # matx=matrix(dat,nrow=10,ncol=4)
#' # RandomizedBlocksAnalysis(matx,con=c(-1,1,-1,1)/2,alpha=0.05)
#' #$n
#' #[1] 10 10 10 10

#' #$test
#' #         test     crit        se       df
#' #[1,] 4.432644 2.038622 0.2798104 31.33793

#' #$psihat
#' #     psihat  ci.lower ci.upper      p.value
#' #[1,] 1.2403 0.6698721 1.810728 0.0001062952

RandomizedBlocksAnalysis <-
  function(x,
           con = c(-0.5, 0.5, -0.5, 0.5),
           alpha = .05) {
    if (base::is.data.frame(x))
      x = base::as.matrix(x)

    flag <- TRUE
    if (alpha != .05 && alpha != .01)
      flag <- FALSE

    if (base::is.matrix(x)) {
      # Corrected to cope with matrix format
      y = list()
      for (j in 1:ncol(x))
        y[[j]] <- x[, j]
      x = list()
      x = y
    }

    if (!base::is.list(x))
      stop("Data must be stored in a matrix or in list mode.")
    con <- base::as.matrix(con)
    if (ncol(con) > 1)
      stop("Only one linear contrast permitted for a standard randomized blocks experiment")
    J <- length(x)
    sam = NA
    h <- base::vector("numeric", J)
    w <- base::vector("numeric", J)
    xbar <- base::vector("numeric", J)
    for (j in 1:J) {
      xx <- !is.na(x[[j]])
      val <- x[[j]]
      x[[j]] <- val[xx]  # Remove missing values
      sam[j] = length(x[[j]])
      h[j] <- length(x[[j]])
      # h is the number of observations in the jth group.
      w[j] <-
        ((length(x[[j]]) - 1) * stats::var(x[[j]])) / (h[j] * (h[j] - 1)) # The variance of the jth group
      xbar[j] <-
        base::mean(x[[j]]) # The mean of the jth group
    }
    #if(sum(con^2)>0){
    if (nrow(con) != length(x)) {
      stop("The number of groups does not match the number of contrast coefficients.")
    }
    psihat <- base::matrix(0, 1, 4)
    dimnames(psihat) <-
      list(NULL, c("psihat", "ci.lower", "ci.upper",
                   "p.value"))
    test <- matrix(0, 1, 4)
    dimnames(test) <- list(NULL, c("test", "crit", "se", "df"))
    df <- 0


    psihat[1, 1] <- sum(con[, 1] * xbar)
    sejk <-
      sqrt(sum(con[, 1] ^ 2 * w)) # The pooled standard error of the contrast

    test[1, 1] <-
      sum(con[, 1] * xbar) / sejk # The value of the t-test
    df <-
      (sum(con[, 1] ^ 2 * w)) ^ 2 / sum(con[, 1] ^ 4 * w ^ 2 / (h - 1)) # Degrees of freeedom allowing for heterogeneity
    vv = (1 - alpha / 2)
    crit <-
      stats::qt(vv, df) #The critical value of the t-test for the degrees of freedom
    test[1, 2] <- crit
    test[1, 3] <- sejk
    test[1, 4] <- df
    psihat[1, 2] <- psihat[1, 1] - crit * sejk
    psihat[1, 3] <- psihat[1, 1] + crit * sejk
    psihat[1, 4] <- 2 * (1 - stats::pt(abs(test[1, 1]), df))


    list(n = sam,
         test = test,
         psihat = psihat)
  }



#' @title Kendalltaupb
#' @description  Computes point bi-serial version of Kendall's tau plus a 1-alpha confidence interval using the method recommended by Long and Cliff (1997).  The algorithm is based on Wilcox's code but was extended to return the consistent variance and the confidence intervals based on the t-distribution. Also added a Diagnostic parameter to output internal calculations.
#' @author Rand Wilcox, Barbara Kitchenham and Lech Madeyski
#' @export Kendalltaupb
#' @param x either a matrix with two columns containg two correlated variables or a vector of variables
#' @param y if y=NULL, assume x is a matrix with two columns, otherwise y is a vector of variables with x[i] and x[i] being from the same experimental unit
#' @param alpha = 0.05, the Type 1 error level used for statistical tests
#' @return list containing the estimate of Kendall's tau, it hypothesis testing variance, and the t-test value obtained from it, the significance of the t-test, the consistent variance of tau and its confidence intervals based on both the normal dstribution and the t-test (recommended by Long and Cliff)
#' @examples
#' x=c(1.2,3,1.8,2,2,0.5,0.5,1,3,1)
#' y=c(1,1,1,1,1,0,0,0,0,0)

#' Kendalltaupb(x,y,alpha=.05)
#' # $cor
#' # [1] 0.3555556
#' # $ci
#' # [1] -0.04198026  0.55555556
#' # $cit
#' # [1] -0.1240567  0.5555556
#' # $test
#' # [1] 1.431084
#' # $sqse
#' # [1] 0.0617284
#' # $consistentvar
#' # [1] 0.04113925
#' # $siglevel
#' # [1] 0.1524063
Kendalltaupb <- function(x, y = NULL, alpha = .05) {
  if (length(x) <= 3)
    stop("Too few data points")
  if (length(x) != length(y))
    stop("Invalid input vectors")
  if (length(x) != length(x[!is.na(x)]))
    stop("Missing values not permitted")
  if (length(y) != length(y[!is.na(y)]))
    stop("Missing values not permitted")

  # Needs a test to ensure one of the variables contains only 1 or 0 values

  m = cbind(x, y)

  x = m[, 1]
  y = m[, 2]
  xdif <- base::outer(x, x, FUN = "-")
  ydif <- base::outer(y, y, FUN = "-")
  tv <- sign(xdif) * sign(ydif)

  #Corrects error in Wilcox's algorithm
  n <- length(x)
  dbar <- base::apply(tv, 1, sum) / (n - 1)

  tau <- sum(tv) / (n * (n - 1))

  A <- sum((dbar - tau) ^ 2) / (n - 1)
  B <-
    (n * (n - 1) * (-1) * tau ^ 2 + sum(tv ^ 2)) / (n ^ 2 - n - 1)
  C <-
    (4 * (n - 2) * A + 2 * B) / (n * (n - 1)) # C is the consistent variance

  # Confidence interval based on normal distribution - not recommended by Long and Cliff
  crit <- stats::qnorm(alpha / 2)
  cilow <- tau + crit * sqrt(C)
  cihi <- tau - crit * sqrt(C)


  if (cilow < (-n / (2 * (n - 1))))
    cilow = -n / (2 * (n - 1)) # Applies limits assuming a point bi-serial tau

  if (cihi > n / (2 * (n - 1)))
    cihi = n / (2 * (n - 1)) # Applies limits assuming a point bi-serial tau


  # Confidence interval based on t distribution - recommended by Long and Cliff
  vv = stats::qt(alpha / 2, n - 3)

  cilowt = tau + vv * sqrt(C)
  if (cilowt < (-n / (2 * (n - 1))))
    cilowt = -n / (2 * (n - 1))
  cihit = tau - vv * sqrt(C)
  if (cihit > n / (2 * (n - 1)))
    cihit = n / (2 * (n - 1))

  # t-test based on hypothesis test variance - not recommended by Long and Cliff
  se = sqrt((2 * (2 * n + 5)) / (9 * n * (n - 1)))
  test <- tau / se

  siglevel <- 2 * (1 - stats::pnorm(abs(test)))


  list(
    cor = tau,
    ci = c(cilow, cihi),
    cit = c(cilowt, cihit),
    test = test,
    sqse = se ^ 2,
    consistentvar = C,
    siglevel = siglevel
  )

}


#' @title Cliffd
#' @description This function implements finds Cliff's d and its confidence intervals. The null hypothesis is that for two independent group, P(X<Y)=P(X>Y). The function reports a 1-alpha confidence interval for P(X>Y)-P(X<Y). The algorithm computes a confidence interval for Cliff's d using the method in Cliff, 1996, p. 140, eq 5.12. The function is based on code produce by Rand Wilcox but has been amended. The plotting function has been removed and the dependency on Wilcox's binomci function has been removed. Construction of confidence intervals if values in one group are all larger than values in the other group has been amended to use the smallest non-zero variance method. Upper and lower confidence interval bounds cannot assume invalid values, i.e. values <-1 or >1.
#' @author Rand Wilcox, amendments Barbara Kitchenham and Lech Madeyski
#' @export Cliffd
#' @param x is a vector of values from group 1
#' @param y is a vector of values from group 2
#' @param alpha is the Type 1 error level for statistical tests
#' @param sigfig is the number of significant digit. If sigfig>0 the data in x and y is truncated to the specified value.
#' @return list including the value of Cliffs d its consistent variance and confidence intervals and the equivalent probability of superiority value and its confidence intervals.
# ' @examples
#' x=c(1.2,3,2.2,4,2.5,3)
#' y=c(3,4.2,4,6,7,5.9)
#' # Cliffd(x,y)
#' #  $n1
#' # [1] 6
#' # $n2
#' # [1] 6
#' # $cl
#' # [1] -0.9772519
#' # $cu
#' # [1] -0.3476461
#' # $d
#' # [1] -0.8611111
#' # $sqse.d
#' # [1] 0.02017931
#' # $phat
#' # [1] 0.06944444
#' # $summary.dvals
#' #        P(X<Y)     P(X=Y)     P(X>Y)
#' # [1,] 0.8888889 0.08333333 0.02777778
#' # $p.cl
#' # [1] 0.01137405
#' # $p.cu
#' # [1] 0.326177
#' #
#' z=c(1,2,3,4)
#' y=c(5,6,7,8)
#' Cliffd(z,y)
#' # $n1
#' # [1] 4
#' # $n2
#' # [1] 4
#' # $cl
#' # [1] -1
#' # $cu
#' # [1] -0.4368172
#' # $d
#' # [1] -1
#' # $sqse.d
#' # [1] 0.009765625
#' # $phat
#' # [1] 0
#' # $summary.dvals
#' #      P(X<Y) P(X=Y) P(X>Y)
#' # [1,]      1      0      0
#' # $p.cl
#' # [1] 0
#' # $p.cu
#' # [1] 0.2815914

Cliffd <- function(x,
                   y,
                   alpha = .05,
                   sigfig = -1) {
  # Check that the data is valid
  if (length(x) <= 1)
    stop("Too few data points")
  if (length(y) <= 1)
    stop("Too few data points")
  if (length(x) != length(x[!is.na(x)]))
    stop("Missing values not permitted")
  if (length(y) != length(y[!is.na(y)]))
    stop("Missing values not permitted")

  # Truncate the data if necessary
  if (sigfig > 0) {
    x = signif(x, sigfig)
    y = signif(y, sigfig)
  }

  m <- base::outer(x, y, FUN = "-")
  msave <- m
  m <- sign(m)

  d <- base::mean(m)

  phat <- (1 + d) / 2

  flag = TRUE
  if (phat == 0 || phat == 1)
    flag = FALSE
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
      di[i] <-
      sum(x[i] > y) / length(y) - sum(x[i] < y) / length(y)

    dh <- NA
    for (i in 1:length(y))
      dh[i] <-
      sum(y[i] < x) / length(x) - sum(y[i] > x) / length(x)


    sdi <- stats::var(di)
    sdh <- stats::var(dh)
    # sh is the consistent variance of d
    sh <-
      ((length(y) - 1) * sdi + (length(x) - 1) * sdh + sigdih) / (length(x) *
                                                                    length(y))


    zv <- stats::qnorm(alpha / 2)
    cu <-
      (d - d ^ 3 - zv * sqrt(sh) * sqrt((1 - d ^ 2) ^ 2 + zv ^ 2 * sh)) /
      (1 - d ^ 2 + zv ^ 2 * sh)
    cl <-
      (d - d ^ 3 + zv * sqrt(sh) * sqrt((1 - d ^ 2) ^ 2 + zv ^ 2 * sh)) /
      (1 - d ^ 2 + zv ^ 2 * sh)
  }
  if (!flag) {
    # phat=0 and d=-1 or phat=1 and d=1
    # Cannot calculate standard error of d, use alternative minimum change approach to calculate a "close" slightly too large estimate. Calculates the smallest non-zero consistent variance
    sh = 0
    di = NA
    dh = NA
    nx = length(x)
    ny = length(y)

    tempn = nx + ny - 1

    if (phat == 1) {
      cu <- 1
      tempn = nx + ny - 1
      disturbedx = c(nx:tempn)
      disturbedy = c(1:nx)
      disturbedanalysis = Cliffd(disturbedx, disturbedy)
      cl = disturbedanalysis$cl

    }

    if (phat == 0) {
      cl = -1
      disturbedy = c(nx:tempn)
      disturbedx = c(1:nx)
      disturbedanalysis = Cliffd(disturbedx, disturbedy)
      cu = disturbedanalysis$cu
    }

    sh = disturbedanalysis$sqse.d
  }
  # Need to construct confidence intervals on phat equivalent to the d confidence intervals
  pci = c((1 + cu) / 2, (1 + cl) / 2)

  cid.results = list(
    n1 = length(x),
    n2 = length(y),
    cl = cl,
    cu = cu,
    d = d,
    sqse.d = sh,
    phat = phat,
    summary.dvals = c.sum,
    p.cl = pci[2],
    p.cu = pci[1]
  )
  return(cid.results)
}

#' @title calculatePhat
#' @description This function extract the probability of superiority (i.e., Phat) and its confidence interval based on Brunner and Munzel (2000) heteroscedastic analog of WMW test. It is based on Wilcox'x bmp function with some amendments. It does not include a plotit facility. It uses the smallest non-zero variance to idetify confidence intervals and statistical significance for values of Phat=0 and Phat=1. It ensure that confidence intervals do not take on invalid values such as values <0 or >1.
#' @author Rand Wilcox amendments by Barbara Kitchenham and Lech Madeyski
#' @export calculatePhat
#' @param x is a vector of values from group 1
#' @param y is a vector of values from group 2
#' @param alpha is the Type 1 error level for statistical tests
#' @param sigfig If sigfig>0 the data in x and y is truncated to the specified number of significant digits.
#' @return list including the value of the t-test for PHat, the estimate of PHat and Cliff's d, and the confidence intervals for PHat.
#' @examples
#' x=c(1.2, 3.0, 2.2, 4.0, 2.5, 3.0)
#' y=c(3,4.2,4,6,7,5.9)
#' calculatePhat(x,y)
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
#' z=c(1,2,3,4)
#' y=c(5,6,7,8)
#' calculatePhat(z,y)
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
                          alpha = .05,
                          sigfig = -1) {
  if (sigfig > 0) {
    x = signif(x, sigfig)
    y = signif(y, sigfig)
  }

  x <- x[!is.na(x)]  # Remove any missing values
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

  flag = TRUE # TRUE if phat ne 0 and phat ne 1
  if (phat == 0 | phat == 1)
    flag = FALSE

  if (flag) {
    Rg1 <- base::rank(x)
    Rg2 <- base::rank(y)
    S1sq <-
      sum((R[flag1] - Rg1 - R1 + (n1 + 1) / 2) ^ 2) / (n1 - 1)
    S2sq <-
      sum((R[flag2] - Rg2 - R2 + (n2 + 1) / 2) ^ 2) / (n2 - 1)
    sig1 <- S1sq / n2 ^ 2
    sig2 <- S2sq / n1 ^ 2
    se <- sqrt(N) * sqrt(N * (sig1 / n1 + sig2 / n2))
    bmtest <- (R2 - R1) / se

    df <-
      (S1sq / n2 + S2sq / n1) ^ 2 / ((S1sq / n2) ^ 2 / (n1 - 1) + (S2sq /
                                                                     n1) ^ 2 / (n2 - 1))
    sig <- 2 * (1 - stats::pt(abs(bmtest), df))
    vv <- stats::qt(alpha / 2, df)



    ci.p <- c(phat + vv * se / N, phat - vv * se / N)

  }
  else {
    # Calculate the smallest non-negative variance and use results to approximate variance and confidence intervals of phat. Gives a more realistic confidence interval than other methods
    Nl1 = N - 1
    newx = c(1:n1)
    newy = c(n1:Nl1)
    newres = calculatePhat(newx, newy)
    se = newres$s.e. * N
    df = newres$df
    sig = newres$sig.level
    vv <- stats::qt(alpha / 2, df)


    ci.p = c(phat + vv * se / N, phat - vv * se / N)
    bmtest = newres$test.stat
    if (sum(x) > sum(y))
      bmtest = bmtest * -1
  }


  # Ensure confidence intervals dont assume impossible values
  if (ci.p[1] < 0)
    ci.p[1] = 0
  if (ci.p[2] > 1)
    ci.p[2] = 1


  list(
    test.stat = bmtest,
    phat = phat,
    dhat = dhat,
    sig.level = sig,
    s.e. = se / N,
    ci.p = ci.p,
    df = df
  )
}

#' @title Calc4GroupNPStats
#' @description This function does a non-parametric analysis of a randomized blocks experiment assuming 2 blocks and 2 treatment conditions.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export Calc4GroupNPStats
#' @param x1 is the data associated with treatment A in one block 1
#' @param x2 is the data associated with treatment B in block 1
#' @param x3 is the data associated with treatment A in block 2
#' @param x4  is the data associated with treatment B in block 2
#' @param sigfig is the number of significant digits in the data. If >0 the datav will be appropriately truncated
#' @param alpha is the significance level for all statistical tests
#' @return The function returns the point biserial version of Kendall's tau and its variance, Cliff's d and its variance, the probability of superiority, phat, and its variance, for the 4 group experiment experiment.
#' @examples
#' set.seed(123)
#' x=list()
#' x[[1]]=stats::rnorm(10,0,1)
#' x[[2]]=stats::rnorm(10,0.8,1)
#' x[[3]]=stats::rnorm(10,0.5,1)
#' x[[4]]=stats::rnorm(10,1.3,1)
#' Calc4GroupNPStats(x[[1]],x[[2]],x[[3]],x[[4]],sigfig=-1,alpha=0.05)
#'# A tibble: 1 x 17
#'#       N  phat phat.var phat.df phat.test phat.pvalue phat.sig      d   vard d.sig    cor   sqse
#'# <int> <dbl>    <dbl>   <dbl>     <dbl>       <dbl> <lgl>     <dbl>  <dbl> <lgl>  <dbl>  <dbl>
#'#   1    40  0.17  0.00497    31.0     -4.68   0.0000532 TRUE     -0.660 0.0206 TRUE  -0.347 0.0132
#'# … with 5 more variables: ctvar <dbl>, n1 <int>, n2 <int>, sigCVt <lgl>, sigCVn <lgl>
Calc4GroupNPStats = function(x1,
                             x2,
                             x3,
                             x4,
                             sigfig = -1,
                             alpha = 0.05) {
  #     Check the significant digits to ensure that equal values are properly detected
  if (sigfig > 0) {
    x1 = signif(x1, sigfig)
    x2 = signif(x2, sigfig)
    x3 = signif(x3, sigfig)
    x4 = signif(x4, sigfig)
  }
  # Set up a dummy variable such that the observations using treatment A in block 1 are associated with the value 1 and observations using treatment B in block 1 are associated with the value 0
  dummy1 = c(rep(1, length(x1)), rep(0, length(x2)))
  # Concatenate the observations in block 1
  xCO1 = c(x1, x2)
  #	Use Wilcox's function to find tau and its two variances for block 1
  tau1 = Kendalltaupb(xCO1, dummy1)
  n1 = length(x1) + length(x2)
  # Set up a dummy variable such that the observations using treatment A in block 2 are associated with the value 1 and observations using treatment B in block 2 are associated with the value 0

  dummy2 = c(rep(1, length(x3)), rep(0, length(x4)))
  # Concatenate the observations in block 12
  xCO2 = c(x3, x4)

  #	Use Wilcox's function to find tau and its two variances for block 2

  tau2 = Kendalltaupb(xCO2, dummy2)
  n2 = length(x3) + length(x4)
  N = n1 + n2

  average.tau = (tau1$cor + tau2$cor) / 2
  combinedsqse = (tau1$sqse + tau2$sqse) / 4

  ctvar = (tau1$consistentvar + tau2$consistentvar) / 4

  # Find the confidence limits on the combined tau using t-distribution
  vv = stats::qt(alpha / 2, N - 6)
  ci.t <-
    c(average.tau + vv * sqrt(combinedsqse),
      average.tau - vv * sqrt(combinedsqse))

  sigCVt = ci.t[1] > 0 | ci.t[2] < 0

  # Find the confidence limits on the combined tau using normal distribution

  ci.n <-
    c(
      average.tau + stats::qnorm(alpha / 2) * sqrt(combinedsqse),
      average.tau - stats::qnorm(alpha / 2) * sqrt(combinedsqse)
    )

  sigCVn = ci.n[1] > 0 | ci.n[2] < 0


  # Find the average d and the combined variances for the full experiment
  d = (Cliffd(x1, x2)$d + Cliffd(x3, x4)$d) / 2
  vard = (Cliffd(x1, x2)$sqse.d + Cliffd(x3, x4)$sqse.d) / 4

  zv <- stats::qnorm(alpha / 2)
  d.cu <-
    (d - d ^ 3 - zv * sqrt(vard) * sqrt((1 - d ^ 2) ^ 2 + zv ^ 2 * vard)) /
    (1 - d ^ 2 + zv ^ 2 * vard)
  d.cl <-
    (d - d ^ 3 + zv * sqrt(vard) * sqrt((1 - d ^ 2) ^ 2 + zv ^ 2 * vard)) /
    (1 - d ^ 2 + zv ^ 2 * vard)
  d.sig = d.cu < 0 | d.cl > 0

  # Find the average phat and the combined variances for the full experiment

  B1.BMP = calculatePhat(x2, x1)
  B2.BMP = calculatePhat(x4, x3)
  phat1 = B1.BMP$phat
  phat2 = B2.BMP$phat
  phat = (phat1 + phat2) / 2
  se1 = B1.BMP$s.e.
  se2 = B2.BMP$s.e.
  se = 4 * N * sqrt((se1 ^ 2 + se2 ^ 2) / 4)
  phat.se = sqrt((B1.BMP$s.e ^ 2 + B2.BMP$s.e. ^ 2) / 4)
  phat.var = phat.se ^ 2
  phat.test = 4 * N * (phat - 0.5) / se
  phat.df = B1.BMP$df + B2.BMP$df
  phat.pvalue = 2 * (1 - stats::pt(abs(phat.test), phat.df))
  phat.sig = phat.pvalue < 0.05


  output = tibble::tibble(
    N = N,
    phat = phat,
    phat.var = phat.var,
    phat.df = phat.df,
    phat.test = phat.test,
    phat.pvalue = phat.pvalue,
    phat.sig = phat.sig,
    d = d,
    vard = vard,
    d.sig = d.sig,
    cor = average.tau,
    sqse = combinedsqse,
    ctvar = ctvar,
    n1 = n1,
    n2 = n2,
    sigCVt = sigCVt,
    sigCVn = sigCVn
  )

  return(output)


}



#' @title LaplaceDist
#' @description Returns a sample of N observations from a Laplace distribution with specified mean and spread.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export LaplaceDist
#' @param N is the required sample size
#' @param mean is the required mean
#' @param spread is the spread of the function
#' @param min lower limit of the distribution, must be finite
#' @param max upper limits of the distribution, must be finite
#' @return N values from a Laplace distribution
#' @examples
#' set.seed(123)
#' LaplaceDist(10,0,1)
#' # [1] -0.55311564  0.85946218 -0.20094937  1.45258293  2.12808209 -2.39565480  0.05785263...
LaplaceDist = function(N,
                       mean,
                       spread,
                       max = 0.5,
                       min = -0.5) {
  y = stats::runif(N, min, max) # Get data from a uniform distribution
  x = mean - spread * sign(y) * log(1 - 2 * abs(y))

  return(x)
}


# Functions for parametric analysis

#' @title ExtractMAStatistics
#' @description This function extracts summary statistics from meta-analysis results obtained from the rma function of the metafor R package. If required the function transform back to standardized mean difference (effect size type "d" i.e. Hg) or point biserial correlations (effect size type "r").
#' Warning: the `ExtractMAStatistics` function works with `metafor` version 2.0-0, but changes to metafor's method of providing access to its individual results may introduce errors into the function.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export ExtractMAStatistics
#' @param maresults is the output from the rma function.
#' @param Nc is the number of participants in the control condition group.
#' @param Nt is the number of participants in the treatment condition group.
#' @param Transform is a boolean value indicating whether the outcome values need to be transformed back to standardized mean difference ("d" i.e. Hg or d) or point biserial correlations ("r"). It is defaulted to TRUE. If this parameter is set to FALSE, no transformation will be applied.
#' @param type this indicates the type of transformation required - it defaults to "d" which requests transformation from Zr to Hg, using "r" requests transformation from Zr to r.
#' @param sig indicates the number of significant digits requested in the output, the default is 4; it rounds the values of mean, pvalue, upper and lower bound to the specified number of significant digits.
#' @param returnse default to FALSE, if set to TRUE returns the standard error of the effect size
#' @return data frame incl. summary statistics from meta-analysis results: overall mean value for the effect sizes, the p-value of the mean, the upper and lower confidence interval bounds (UB and LB), QE which is the heterogeneity test statistic and QEp which the the p-value of the heterogeneity statistic
#' @examples
#' ExpData=reproducer::KitchenhamMadeyskiBrereton.ExpData
#' #Extract the experiment basic statics
#' S1data=subset(ExpData,ExpData=="S1")
#' #Use the descriptive data to construct effect size
#' S1EffectSizes = reproducer::PrepareForMetaAnalysisGtoR(
#'   S1data$Mc,S1data$Mt,S1data$SDc,S1data$SDt,S1data$Nc,S1data$Nt)
#' # Do a random effect meta-analysis of the transformed r_pbs effect size
#' S1MA = metafor::rma(S1EffectSizes$zr, S1EffectSizes$vi)
#' # Extract summary statistics from meta-analysis results and transform back to Hg scale
#' ExtractMAStatistics(S1MA, sum(S1data$Nc),sum(S1data$Nt), TRUE, "d", 4)
#' #   A tibble: 1 x 6
#' #      mean  pvalue    UB    LB    QE   QEp
#' #     <dbl>   <dbl> <dbl> <dbl> <dbl> <dbl>
#' #   1 0.666 0.00207  1.12 0.238     4  0.41

ExtractMAStatistics = function(maresults,
                               Nc,
                               Nt,
                               Transform = TRUE,
                               type = "d",
                               sig = 4,
                               returnse = FALSE) {
  pvalue = as.numeric(maresults$pval)

  se = as.numeric(maresults$se)

  QE = as.numeric(maresults$QE)
  QEp = as.numeric(maresults$QEp)
  mean = as.numeric(maresults$beta)
  UB = as.numeric(maresults$ci.ub)
  LB = as.numeric(maresults$ci.lb)

  if (Transform & type == "d") {
    mean = reproducer::transformZrtoHg(mean, Nc, Nt)
    se = reproducer::transformZrtoHg(se, Nc, Nt)
    UB = reproducer::transformZrtoHg(UB, Nc, Nt)
    LB = reproducer::transformZrtoHg(LB, Nc, Nt)

  }
  if (Transform & type == "r") {
    mean = reproducer::transformZrtoR(mean)
    se = reproducer::transformZrtoR(se)

    UB = reproducer::transformZrtoR(UB)
    LB = reproducer::transformZrtoR(LB)
  }
  mean = signif(mean, sig)
  pvalue = signif(pvalue, sig)
  se = signif(se, sig)

  UB = signif(UB, sig)
  LB = signif(LB, sig)
  QE = signif(QE, 2)
  QEp = signif(QEp, 2)
  if (returnse)
    metaanalysisresults = tibble::tibble(mean, pvalue, se, UB, LB, QE, QEp)
  else
    metaanalysisresults = tibble::tibble(mean, pvalue, UB, LB, QE, QEp)
  return(metaanalysisresults)

}



###############################################################################
# Simulation functions


#' @title simulateRandomizedDesignEffectSizes
#' @description This simulates one of four distributions, and finds the values of ktau, phat and Cliffs d and their variances. It assumes equal group sizes. It returns values of the effect sizes and their variance for a simulated randomized experiment with two treatments.  It returns whether to not each non-parametric effect size was significant. It also returns the parametric (unstandardized and unstandardized) Effect Size and the whether the t-test was signficiant.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export simulateRandomizedDesignEffectSizes
#' @param mean The mean used for one of the treatment groups
#' @param sd The spread used for both treatment groups. It mus be a real value greater than 0.
#' @param diff This is added to the parameter mean, to define the mean of the other treatment group. It can be a real value avd can take the value zero.
#' @param N this is the number of observations in each group. It must be an integer greater than 3.
#' @param type this specifies the underlying distribution used to generate the data. it takes the values "n" for a normal distribution, "l" for lognormal distribution,"g" for a gamma distribution, "lap" for a Laplace distribution.
#' @param StdAdj this specifes the extent of variance instability introduced by the treatment.
#' @return data frame incl. the non-parametric and parametric effect sizes and whether the effect sizes are significant at the 0.05 level.
#' @examples
#' set.seed(123)
#' simulateRandomizedDesignEffectSizes(mean=0,sd=1,diff=0.8,N=10,type="n",StdAdj=0)
#' # A tibble: 1 x 15
#' # phat varphat dfphat sigphat     d   vard sigd    cor varcor sigCVt ttestp    ES Variance StdES
#' # MedDiff
#' # <dbl>   <dbl>  <dbl> <lgl>   <dbl>  <dbl> <lgl> <dbl>  <dbl> <lgl>   <dbl> <dbl>    <dbl> <dbl>
#' # <dbl>
#' #1  0.75  0.0152   17.5 FALSE     0.5 0.0624 FALSE 0.263 0.0175 FALSE  0.0507 0.934  0.994 0.937
#' # 1.26
#' set.seed(123)

#' simulateRandomizedDesignEffectSizes(mean=0,sd=1,diff=0.8,N=10,type="l",StdAdj=0)
#' # A tibble: 1 x 19
#' #   phat varphat dfphat sigphat     d   vard sigd    cor varcor sigCVt ttestp    ES Variance
#' # StdES MedDiff transttest
#' #  <dbl>   <dbl>  <dbl> <lgl>   <dbl>  <dbl> <lgl> <dbl>  <dbl> <lgl>   <dbl> <dbl>    <dbl>
#' # <dbl>   <dbl>      <dbl>
#' #1  0.75  0.0152   17.5 FALSE     0.5 0.0624 FALSE 0.263 0.0175 FALSE  0.0958  2.41     9.01
#' # 0.802    2.32     0.0507
#' # … with 3 more variables: EStrans <dbl>, StdEStrans <dbl>, VarTrans <dbl>

simulateRandomizedDesignEffectSizes = function(mean,
                                              sd,
                                              diff,
                                              N,
                                              type = "n",
                                              StdAdj = 0) {
  if (type == "n")
  {
    x = stats::rnorm(N, mean, sd)
    y = stats::rnorm(N, mean + diff, sd + StdAdj)
  }
  if (type == "g")
  {
    rate = mean
    shape = sd
    x = stats::rgamma(N, shape, rate)
    y = stats::rgamma(N, shape + StdAdj, rate + diff)
  }
  if (type == "l")
  {
    x = stats::rlnorm(N, mean, sd)
    y = stats::rlnorm(N, mean + diff, sd + StdAdj)
    transx = log(x)
    transy = log(y)
  }

  if (type == "lap")
  {
    # Changed to use built-in defaults for max and min
    x = LaplaceDist(N, mean, sd)
    y = LaplaceDist(N, mean + diff, sd + StdAdj)
  }

  # Calculate the Kendall's tau value for the data

  dummy = c(rep(0, N), rep(1, N))
  xy = c(x, y)
  # Calculate ktau for the data set
  ktest = Kendalltaupb(dummy, xy)

  ktau = ktest$cor
  sigCVt = ktest$cit[1] > 0 | ktest$cit[2] < 0

  # Calculate Cliff's d for the data set and determine whether the value is significant
  ctest = Cliffd(y, x)

  d = ctest$d
  dvar = ctest$sqse

  d.lb = ctest$cl
  d.ub = ctest$cu

  dsig = d.lb > 0 | d.ub < 0

  # Calculate phat statistics for the data set

  ptest = calculatePhat(x, y)

  phat = ptest$phat

  phat.sig = ptest$sig.level < 0.05
  phat.df = ptest$df
  phat.var = ptest$s.e. ^ 2


  # Add the results of a t-test as a baseline

  res = stats::t.test(x, y)

  UES = base::mean(y) - base::mean(x)
  Var = (stats::var(x) + stats::var(y)) / 2
  MedDiff = stats::median(y) - stats::median(x)
  EffectSize = UES / sqrt(Var)
  pval = res$p.value

  if (type == "l")
  {
    # Check that log-normal datra gives appropriate values after transformation
    trans.t = stats::t.test(transx, transy)
    pval.trans = trans.t$p.value
    ES.trans = base::mean(transy) - base::mean(transx)
    VarTrans = (stats::var(transx) + stats::var(transy)) / 2
    StdES.trans = ES.trans / sqrt(VarTrans)
  }

  if (type == "l")
    output = tibble::tibble(
      phat = phat,
      varphat = phat.var,
      dfphat = phat.df,
      sigphat = phat.sig,
      d = d,
      vard = dvar,
      sigd = dsig,
      cor = ktau,
      varcor = ktest$consistentvar,
      sigCVt = sigCVt,
      ttestp = pval,
      ES = UES,
      Variance = Var,
      StdES = EffectSize,
      MedDiff = MedDiff,
      transttest = pval.trans,
      EStrans = ES.trans,
      StdEStrans = StdES.trans,
      VarTrans = VarTrans
    )
  else
    output = tibble::tibble(
      phat = phat,
      varphat = phat.var,
      dfphat = phat.df,
      sigphat = phat.sig,
      d = d,
      vard = dvar,
      sigd = dsig,
      cor = ktau,
      varcor = ktest$consistentvar,
      sigCVt = sigCVt,
      ttestp = pval,
      ES = UES,
      Variance = Var,
      StdES = EffectSize,
      MedDiff = MedDiff
    )
  return(output)
}

############################################################################################################################
#############################################################################################################################


#' @title RandomExperimentSimulations
#' @description This function performs multiple simulations of two-group balanced experiments for one of four distributions and a specific group size. It identifies the average value of phat, Cliff' d and Kendall's point biserial tau and their variances. It either returns the effect sizes for each non-parametric effect size or it reports the number of times the each non-parametric effect size is assessed to be significantly different from zero. We also present the values for the t-test as a comparison. For log-normal data the results of analysing the transformed data are also reported.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export RandomExperimentSimulations
#' @param mean The default mean used for both groups (one treatment group and one control group). It can be changed for the treatment group using  the parameter diff
#' @param sd This is the default spread for both groups. It must be a real value greater than 0. It can be adjusted for the treatment group using the parameter StdAdj
#' @param diff This is added to the treatment group mean. It can be a real value and can take the value zero.
#' @param N this is the number of observations in each group. It must be an integer greater than 3.
#' @param reps this identifies the number of times each experiment simulation is replicated.
#' @param type this specifies the underlying distribution used to generate the data. It takes the values "n" for a normal distribution, "l" for lognormal distribution,"g" for a gamma distribution, "lap" for a Laplace distribution.
#' @param seed This specifies the initial seed for the set of replications (default 123).
#' @param StdAdj this specifes the extent of variance instability introduced by the treatment and it must be non-negative but can be 0.
#' @param returnData If TRUE, the function returns the individual effect sizes and their variances, otherwise it returns summary statistics (default FALSE).
#' @examples
#' #RandomExperimentSimulations(mean=0,sd=1,diff=0.5,N=20, reps=500,type="n",seed=123,StdAdj=0)
#' #reps=500 may take a bit too for CRAN so let us use reps=50
#' RandomExperimentSimulations(mean=0,sd=1,diff=0.5,N=20, reps=50,type="n",seed=123,StdAdj=0)
#'# A tibble: 1 x 17
#'#   phat phatvar sigphat emp.phat.var     d   dvar  sigd emp.d.var  ktau ktauvar emp.tau.var
#'# kpowerCVt tpower    ES
#'#  <dbl>   <dbl>   <dbl>        <dbl> <dbl>  <dbl> <dbl>     <dbl> <dbl>   <dbl>       <dbl>
#'# <dbl>  <dbl> <dbl>
#'#1 0.639 0.00793    0.33      0.00790 0.278 0.0324 0.306    0.0316 0.143 0.00854     0.00832
#'# 0.328  0.338 0.500
#'# … with 3 more variables: Variance <dbl>, StdES <dbl>, MedDiff <dbl>

#' #RandomExperimentSimulations(mean=0,sd=1,diff=0.5,N=20, reps=500,type="l",seed=123,StdAdj=0)
#' #reps=500 may take a bit too for CRAN so let us use reps=50
#' RandomExperimentSimulations(mean=0,sd=1,diff=0.5,N=20, reps=50,type="l",seed=123,StdAdj=0)
#' # A tibble: 1 x 20
#' #   phat phatvar sigphat emp.phat.var     d   dvar  sigd emp.d.var  ktau ktauvar emp.tau.var
#' # kpowerCVt tpower    ES
#' #  <dbl>   <dbl>   <dbl>        <dbl> <dbl>  <dbl> <dbl>     <dbl> <dbl>   <dbl>       <dbl>
#' # <dbl>  <dbl> <dbl>
#' #1 0.639 0.00793    0.33      0.00790 0.278 0.0324 0.306    0.0316 0.143 0.00854     0.00832
#' # 0.328  0.218  1.05
#' # … with 6 more variables: Variance <dbl>, StdES <dbl>, MedDiff <dbl>, ESLog <dbl>,
#' # StdESLog <dbl>, VarLog <dbl>

#' RandomExperimentSimulations(mean=0,sd=1,diff=0.5,N=20, reps=10,type="n",seed=123,StdAdj=0,
#' returnData=TRUE)
#' #  A tibble: 10 x 3
#' #    Cliffd  PHat   StdES
#' #     <dbl> <dbl>   <dbl>
#' #  1  0.175 0.588  0.340
#' #  2  0.19  0.595  0.283
#' #  3 -0.125 0.438 -0.305
#' #  4  0.195 0.597  0.345
#' #  5  0.15  0.575  0.200
#' #  6  0.3   0.65   0.430
#' #  7  0.455 0.728  0.762
#' #  8  0.015 0.507  0.0622
#' #  9  0.16  0.58   0.340
#' # 10  0.2   0.6    0.336

RandomExperimentSimulations = function(mean,
                                       sd,
                                       diff,
                                       N,
                                       reps,
                                       type = "n",
                                       seed = 123,
                                       StdAdj = 0,
                                       returnData = FALSE) {
  phatsum = 0 # This is used to sum the value of phat across the replications
  phatvarsum = 0  # This is used to sum the value of the variance of phat across the replications
  sig.phat = 0 # This is used to sum the number of times pat is significant across the replications
  phatsquare = 0 # This  is used to sum phat^2 and construct an empirical variance of phat
  dsum = 0 # This is used to sum the value of Cliff's d across the replications
  dvarsum = 0  # This is used to sum the value of the variance of Cliff's d across the replications
  sig.d = 0 # This is used to sum the number of times d is significant across the replications
  dsquare = 0 # This  is used to sum d^2 and construct an empirical variance of d.

  ksum = 0 # This is used to sum the value of the point biserial tau across the replications
  kvarsum = 0 # This is used to sum the value of the variance of the point biserial tau across the replications
  ksquare = 0  # This is used to sum the square of tau_pb across the replications and construct an empirical variance
  ksigCVt = 0
  tsig = 0 # This is used to count the number of significant t values across the replications
  ES = 0 # This is used to sum the value of the parametric effect size (unstandardized) across the replications
  StdES = 0 # This is used to sum the value of the parametric effect size (standardized) across the replications
  Var = 0 # This is used to sum the value of the variance across replications
  ES.l.trans = 0 # This is used to sum the transformed unstandardized effect size for lognormal data sets
  StdES.l.trans = 0 # This is used to sum the transformed standardized effect size for lognormal data sets
  Var.l.trans = 0 # This is used to sum the transformed variance for lognormal data sets


  MedDiff = 0 # Used to hold the median difference

  DataTable = NULL

  base::set.seed(seed)
  for (i in 1:reps) {
    # Call the program that generates the random data sets and calculates the sample statistics.
    res = simulateRandomizedDesignEffectSizes(mean, sd, diff, N, type, StdAdj)

    if (returnData == FALSE) {
      # Aggregate data to provide counts of significance and overall effect size averages
      # Cliff's d
      dsum = dsum + res$d
      dvarsum = dvarsum + res$vard
      if (res$sigd)
        sig.d = sig.d + 1
      dsquare = dsquare + res$d ^ 2

      # Probability of superiority
      phatsum = phatsum + res$phat
      phatvarsum = phatvarsum + res$varphat
      if (res$sigphat)
        sig.phat = sig.phat + 1
      phatsquare = phatsquare + res$phat ^ 2

      # Point biserial Kendall's tau
      ksum = ksum + res$cor
      kvarsum = kvarsum + res$varcor
      ksquare = ksquare + res$cor ^ 2
      ksigCVt = ksigCVt + if (res$sigCVt)
        1
      else
        0

      # Parametric statistics
      ES = ES + res$ES
      StdES = StdES + res$StdES
      MedDiff = MedDiff + res$MedDiff
      tsig = tsig + if (res$ttestp < .05)
        1
      else
        0
      Var = Var + res$Variance

      if (type == "l") {
        ES.l.trans = ES.l.trans + res$EStrans
        StdES.l.trans = StdES.l.trans + res$StdEStrans
        Var.l.trans = Var.l.trans + res$VarTrans
      }
    }
    else {
      # Store the outcome from each replication
      # DataTable = tibble::tibble(base::rbind(
      #  DataTable = base::rbind(
      #     DataTable,
      #     base::cbind(
      #       Cliffd = res$d,
      #       PHat = res$phat,
      #       StdES = res$StdES
      # #    )
      #   ))
      DataTable = tibble::tibble(dplyr::bind_rows(
        DataTable,
        #base::cbind(
        dplyr::bind_cols(
          Cliffd = res$d,
          PHat = res$phat,
          StdES = res$StdES
        )
      ))
    }

  }

  if (returnData == FALSE) {
    #Calculate averages of the statistics across the replications
    d = dsum / reps
    dvar = dvarsum / (reps)
    sigd = sig.d / (reps)
    emp.d.var = (dsquare - reps * d ^ 2) / (reps - 1)

    phat = phatsum / reps
    phatvar = phatvarsum / (reps)
    sigphat = sig.phat / (reps)
    emp.phat.var = (phatsquare - reps * phat ^ 2) / (reps - 1)


    ktau = ksum / reps
    ktauvar = kvarsum / reps
    emp.tau.var = (ksquare - reps * ktau ^ 2) / (reps - 1)

    kpowerCVt = ksigCVt / reps

    ES = ES / reps
    StdES = StdES / reps
    MedDiff = MedDiff / reps
    Variance = Var / reps
    tpower = tsig / reps

    if (type == "l")
    {
      # This is used for validation that the algorithms are consistent. The statistics from the transformed lognormal data can be compared with the statistics from the normal data.
      ESLog = ES.l.trans / reps
      StdESLog = StdES.l.trans / reps
      VarLog = Var.l.trans / reps
    }

    if (type == "l")	{
      outcome = tibble::tibble(
        phat,
        phatvar,
        sigphat,
        emp.phat.var,
        d,
        dvar,
        sigd,
        emp.d.var,
        ktau,
        ktauvar,
        emp.tau.var,
        kpowerCVt,
        tpower,
        ES,
        Variance,
        StdES,
        MedDiff,
        ESLog = ESLog,
        StdESLog = StdESLog,
        VarLog = VarLog
      )
    }
    else {
      outcome = tibble::tibble(
        phat,
        phatvar,
        sigphat,
        emp.phat.var,
        d,
        dvar,
        sigd,
        emp.d.var,
        ktau,
        ktauvar,
        emp.tau.var,
        kpowerCVt,
        tpower,
        ES,
        Variance,
        StdES,
        MedDiff
      )
    }
  }
  else {
    #outcome = DataTable
    DataTable = as.data.frame(DataTable)
    outcome = tibble::tibble(DataTable)
  }
  return(outcome)

}



#' @title simulateRandomizedBlockDesignEffectSizes
#' @description This simulates one of four distributions, and finds the values of ktau and Cliffs d and their variances. It simulates a randomised blocks experiment with two treatment groups and two control groups each of which being divided into two blocks. By default it assumes equal group sizes but  group spread (standard deviation can be changed, see Stadj). It returns values of both parametric and non=parametric effect sizes and their variance for simulated experiments. It returns the number of times each effect size was signficiant. For the logarithmic distribution it calcates effect sizes based on the log transformed data as well as the raw data.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export simulateRandomizedBlockDesignEffectSizes
#' @param mean The default value for all groups which can be changed for the two treatment groups using the parameter diff and for the two block 2 groups using the parameter Blockmean
#' @param sd The default spread used for all four groups unless adjusted by the StdAdj. It must be a real value greater than 0.
#' @param diff This is added to the parameter mean to obtain the required mean for treatment groups. It can be a real value and can take the value zero.
#' @param N this is the number of observations in each group. It must be an integer greater than 3.
#' @param type this specifies the underlying distribution used to generate the data. it takes the values "n" for a normal distribution, "l" for lognormal distribution,"g" for a gamma distribution, "lap" for a Laplace distribution.
#' @param alpha is the significance level for statistical tests
#' @param Blockmean if >0 an adjusment made to both group means in Block 2
#' @param BlockStdAdj if >0, an adjustment that can be made to the sd of each group in block 2
#' @param StdAdj this specifies the extent of variance instability introduced by the treatment and if >0 will be used to amend the sd parameter for both treatment groups.
#' @return data frame incl. the non-parametric and parametric effect sizes and whether the effect sizes are significant at the 0.05 level.
#' @examples
#' set.seed(123)
#' simulateRandomizedBlockDesignEffectSizes(mean=0,sd=1,diff=.5,N=10,type="n",alpha=0.05,
#' Blockmean=0.5,BlockStdAdj=0,StdAdj=0)
#' # A tibble: 1 x 23
#' #      N  phat phat.var phat.df phat.test phat.pvalue phat.sig     d   vard d.sig   cor   sqse
#' # ctvar    n1    n2
#' #  <int> <dbl>    <dbl>   <dbl>     <dbl>       <dbl> <lgl>    <dbl>  <dbl> <lgl> <dbl>  <dbl>
#' # <dbl> <int> <int>
#' #1    40  0.79  0.00587    30.2      3.79    0.000681 TRUE      0.58 0.0243 TRUE  0.305 0.0132
#' # 0.00695    20    20
#' ## … with 8 more variables: sigCVt <lgl>, sigCVn <lgl>, ttest.sig <lgl>, ES <dbl>,
#' # Variance <dbl>, StdES <dbl>,
#' #   BlockEffect <dbl>, MedianDiff <dbl>
#' set.seed(123)
#' simulateRandomizedBlockDesignEffectSizes(mean=0,sd=1,diff=.5,N=10,type="l",alpha=0.05,
#' Blockmean=0.5,BlockStdAdj=0,StdAdj=0)
#' # A tibble: 1 x 27
#' #      N  phat phat.var phat.df phat.test phat.pvalue phat.sig     d   vard d.sig   cor   sqse
#' # ctvar    n1    n2
#' #  <int> <dbl>    <dbl>   <dbl>     <dbl>       <dbl> <lgl>    <dbl>  <dbl> <lgl> <dbl>  <dbl>
#' # <dbl> <int> <int>
#' #1    40  0.79  0.00587    30.2      3.79    0.000681 TRUE      0.58 0.0243 TRUE  0.305 0.0132
#' # 0.00695    20    20
#' # … with 12 more variables: sigCVt <lgl>, sigCVn <lgl>, ttest.sig <lgl>, ES <dbl>,
#' # Variance <dbl>, StdES <dbl>,
#' #   BlockEffect <dbl>, MedianDiff <dbl>, Log.sig <lgl>, ES.Trans <dbl>, StdES.Trans <dbl>,
#' # VarTrans <dbl>

simulateRandomizedBlockDesignEffectSizes = function(mean,
                                                    sd,
                                                    diff,
                                                    N,
                                                    type = "n",
                                                    alpha = 0.05,
                                                    Blockmean = 0,
                                                    BlockStdAdj = 0,
                                                    StdAdj = 0) {
  # Generate data. x and x2 hold control data, y and y2 to hold treatment data.  x and y are ib block 1 and x2 and y2 are in block 2

  if (type == "n")
  {
    x = stats::rnorm(N, mean, sd)
    y = stats::rnorm(N, mean + diff, sd + StdAdj)
    x2 = stats::rnorm(N, mean + Blockmean, sd + BlockStdAdj)
    y2 = stats::rnorm(N, mean + diff + Blockmean, sd + BlockStdAdj + StdAdj)

  }
  if (type == "g")
  {
    shape = sd
    rate = mean
    # For a Gamma distribution shapeand rate must always be greater than zero

    if (BlockStdAdj == 0)
      BlockStdAdj = 1

    x = stats::rgamma(N, shape, rate)
    y = stats::rgamma(N, shape + StdAdj, rate + diff)
    x2 = stats::rgamma(N, shape * BlockStdAdj + Blockmean, rate)
    y2 = stats::rgamma(N, shape * BlockStdAdj + StdAdj + Blockmean, rate + diff)
  }
  if (type == "l")
  {
    x = stats::rlnorm(N, mean, sd)
    y = stats::rlnorm(N, mean + diff, sd + StdAdj)
    x2 = stats::rlnorm(N, mean + Blockmean, sd + BlockStdAdj)
    y2 = stats::rlnorm(N, mean + diff + Blockmean, sd + BlockStdAdj + StdAdj)
  }

  if (type == "lap")
  {
    x = LaplaceDist(N, mean, sd)
    y = LaplaceDist(N, mean + diff, sd + StdAdj)
    x2 = LaplaceDist(N, mean + Blockmean, sd + BlockStdAdj)
    y2 = LaplaceDist(N, mean + diff + Blockmean, sd + BlockStdAdj + StdAdj)

  }

  res = Calc4GroupNPStats(y, x, y2, x2, alpha = alpha)

  res = tibble::as_tibble(res)

  # Add the results of a t-test as a baseline
  UES = (base::mean(y) + base::mean(y2) - base::mean(x) - base::mean(x2)) /
    2

  BlockEffect = (base::mean(x2) + base::mean(y2) - base::mean(x) - base::mean(y)) /
    2

  MedDiff = (stats::median(y) + stats::median(y2) - stats::median(x) - stats::median(x2)) /
    2



  Var = (stats::var(y) + stats::var(y2) + stats::var(x) + stats::var(x2)) /
    4
  # Estimate of Cohen's d
  StdES = UES / sqrt(Var)


  # Need to use the linear contrast method for the significance test of randomized blocks data that allows for variance heterogeneity.
  newlist = list()
  newlist[[1]] = x
  newlist[[2]] = y
  newlist[[3]] = x2
  newlist[[4]] = y2
  vec = c(-1, 1, -1, 1) / 2
  res.t = RandomizedBlocksAnalysis(newlist, con = vec, alpha = 0.05)
  df.ttest = base::data.frame(res.t$psihat)
  pval = df.ttest$p.value
  ttest.sig = pval < 0.05

  StandardMetrics = list(
    ttest.sig = ttest.sig,
    ES = UES,
    Variance = Var,
    StdES = StdES,
    BlockEffect = BlockEffect,
    MedianDiff = MedDiff
  )
  StandardMetrics = tibble::as_tibble(StandardMetrics)

  if (type == "l")
  {
    loglist = list()

    loglist[[1]] = log(x)
    loglist[[2]] = log(y)
    loglist[[3]] = log(x2)
    loglist[[4]] = log(y2)

    vec = c(-1, 1, -1, 1) / 2
    logres.t = RandomizedBlocksAnalysis(loglist, con = vec, alpha = 0.05)
    df.logttest = base::data.frame(logres.t$psihat)
    pval.log = df.logttest$p.value
    Log.sig = pval.log < 0.05
    ES.trans = (base::mean(log(y)) + base::mean(log(y2)) - base::mean(log(x)) -
                  base::mean(log(x2))) / 2

    VarTrans = (stats::var(log(x)) + stats::var(log(x2)) + stats::var(log(y)) +
                  stats::var(log(y2))) / 4
    StdES.trans = ES.trans / sqrt(VarTrans)

    AdditionalMetrics = list(
      Log.sig = Log.sig,
      ES.Trans = ES.trans,
      StdES.Trans = StdES.trans,
      VarTrans = VarTrans
    )
    AdditionalMetrics = tibble::as_tibble(AdditionalMetrics)
    StandardMetrics = dplyr::bind_cols(StandardMetrics, AdditionalMetrics)

  }
  res = dplyr::bind_cols(res, StandardMetrics)

  return(res)

}


#' @title RandomizedBlocksExperimentSimulations
#' @description This function performs multiple simulations of 4 group balanced randomised Block experiments with two control groups and two treatment groups where one control group and one treatment group are assigned to block 1 and the other control group and treatment group are assigned to block 2.  The simulations are based on one of four distributions and a specific group size. The function identifies the average value of the non-paramtric effect sizes P-hat, Cliff' d and Kendall's point biserial tau and their variances and whether ot not the statistics were significant at the 0.05 level. We also present the values of the t-test as a comparison.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export RandomizedBlocksExperimentSimulations
#' @param mean The default mean for all 4 groups. The default for the two treatment groups can be altered using the parameter diff and the block mean for block 2 can be altered using the parameter Blockmean.
#' @param sd The default spread for all 4 groups. It must be a real value greater than 0. If can be altered for treatment groups using the parameter StdAdj and for Block 2 groups using BlockStdAdj
#' @param diff The is is added to the parameter mean, to define the mean of the other treatment group. It can be a real value ad can take the value zero.
#' @param N this is the number of observations in each group. It must be an interger greater than 3.
#' @param reps this identifies the number of tiume the simulation is replicated.
#' @param type this specifies the underlying distribution used to generate the data. it takes the values "n" for a normal distribution, "l" for lognormal distribution,"g" for a gamma distribution, "lap" for a Laplace distribution.
#' @param alpha is the Type 1 error level used for constructing confidence intervals
#' @param Blockmean is the effect of having two different blocks
#' @param BlockStdAdj is the variance associated with the Block mean. If Blockvar is zero it means we are treat the block effect as a fixed effect. If BlockStdAdj>0, we treat the block effect as a random effect.
#' @param StdAdj The value used to introduce heterogeneity into the treatment groups variance if required.
#' @param seed this specifies the seed value for the simulations and allows the experiment to be repeated.
#' @param returnData = FALSE if TRUE the function returns the generated data otherwise it returns summary statistics.
#' @return depending on the parameter returnData it returns the generated nonparametric and parametric values or the summary statistics
#' examples
#' RandomizedBlocksExperimentSimulations(mean=0,sd=1,diff=0.5,N=10, reps=500,type="n",alpha=0.05,Blockmean=0.5, BlockStdAdj=0,StdAdj=0,seed=123)
#' # A tibble: 1 x 18
#' #    phat varphat sigphat emp.phat.var     d   vard  sigd emp.d.var  ktau kconsistentvar emp.tau.var kpowerCVt StdES
#' #   <dbl>   <dbl>   <dbl>        <dbl> <dbl>  <dbl> <dbl>     <dbl> <dbl>          <dbl>       <dbl>     <dbl> <dbl>
#' # 1 0.640 0.00832   0.326      0.00773 0.279 0.0343 0.282    0.0309 0.147        0.00956     0.00856     0.168 0.513
#' #  … with 5 more variables: ES <dbl>, Var <dbl>, emp.StdESvar <dbl>, MedDiff <dbl>, tpower <dbl>

#' RandomizedBlocksExperimentSimulations(mean=0,sd=1,diff=0.5,N=10, reps=10,type="n",alpha=0.05,Blockmean=0.5, BlockStdAdj=0,StdAdj=0,seed=123,returnData=T)
#' # A tibble: 10 x 3
#' #   Cliffd  PHat StdES
#' #    <dbl> <dbl> <dbl>
#' # 1  0.58  0.79  1.06
#' # 3  0.37  0.685 0.761
#' # 4  0.440 0.72  0.821
#' # 5  0.13  0.565 0.240
#' # 6  0.16  0.58  0.222
#' # 7  0.38  0.69  0.580
#' # 8  0.48  0.74  0.882
#' # 9  0.11  0.555 0.181
#' #10 -0.03  0.485 0.124


RandomizedBlocksExperimentSimulations = function(mean,
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
                                                 returnData = FALSE) {
  ksum = 0 # This is used to sum the value of the point biserial tau across the replications
  k.ss = 0 # This is used to sum the square of the point biserial tau across the replications. It can be used to calculate the empirical variance of the average tau.

  kconsistentvar = 0 # This is used to sum the variance BiSerial Kendall's tau  values  across the replications
  ksigCVt = 0  # This is used to sum the number of significant Point BiSerial Kendall's tau  values  across the replications


  dsum = 0 # This is used to sum the value of Cliff's d across the replications
  dvarsum = 0 # This is used to aggregate the variance of d across the replications
  d.sig = 0 # This is used to sum the number of significant d values  across the replications
  d.ss = 0 # This is used to sum the squared value of Cliff's d across the replications. It can be used to calculate the empirical variance of the average d.

  phatsum = 0 # This is used to sum the value of phat across the replications
  phatvarsum = 0 # This is used to aggregate the variance of phat across the replications
  phat.sig = 0 # This is used to sum the number of significant phat values  across the replications
  phat.ss = 0  # This is used to sum the squared value of phat across the replications. It can be used to calculate the empirical variance of the average phat.


  tsig = 0 # This is used to count the number of significant t values across the replications
  ES = 0 # This is used to sum the value of the parametric effect size (unstandardized) across the replications
  StdES = 0 # This is used to sum the value of the parametric effect size (standardized) across the replications
  StdES.ss = 0 # This is used to sum the squared value of the parametric effect size (standardized) across the replications. It can be used to calculate the empirical variance of the overall mean value

  Var = 0 # This is used to sum the estimate of the pooled variance across replications
  MedDiff = 0 # This is used to sum the median across replications

  ES.l.trans = 0 # This is used to sum the unstandardized effect size of the log-normal data after being transformed
  StdES.l.trans = 0 # This is used to sum the standardized effect size of the log-normal data after being transformed
  Var.l.trans = 0 # This is used to sum the variance of the log-normal data after being transformed

  trans.sig = 0  # This is used to sum the number of signficant t-tests of the transformed lognormal data

  base::set.seed(seed)

  DataTable = NULL

  for (i in 1:reps) {
    # Call the program than generates the randomized block experiment data sets and calculates the sample statistics
    res = simulateRandomizedBlockDesignEffectSizes(
      mean,
      sd,
      diff,
      N,
      type,
      alpha = alpha,
      Blockmean = Blockmean,
      BlockStdAdj = BlockStdAdj,
      StdAdj = StdAdj
    )
    if (!returnData) {
      # Cliff's d
      dsum = dsum + res$d
      dvarsum = dvarsum + res$vard
      if (res$d.sig)
        d.sig = d.sig + 1
      d.ss = d.ss + res$d ^ 2

      # Probability of Speriority (phat)
      phatsum = phatsum + res$phat
      if (res$phat.sig)
        phat.sig = phat.sig + 1
      phat.ss = phat.ss + res$phat ^ 2
      phatvarsum = phatvarsum + res$phat.var

      # Point Biserial Kendall's tau
      ksum = ksum + res$cor
      k.ss = k.ss + res$cor ^ 2
      kconsistentvar = kconsistentvar + res$ctvar
      ksigCVt = ksigCVt + if (res$sigCVt)
        1
      else
        0

      # Standard parametric effect sizes
      ES = ES + res$ES
      StdES = StdES + res$StdES
      StdES.ss = StdES.ss + res$StdES ^ 2
      Var = Var + res$Variance
      MedDiff = MedDiff + res$MedianDiff
      tsig = tsig + if (res$ttest.sig)
        1
      else
        0

      if (type == "l")
      {
        trans.sig = trans.sig + if (res$Log.sig)
          1
        else
          0

        ES.l.trans = ES.l.trans + res$ES.Trans
        StdES.l.trans = StdES.l.trans + res$StdES.Trans
        Var.l.trans = Var.l.trans + res$VarTrans
      }

    }
    else {
      # Store the outcome from each replication
      DataTable = tibble::tibble(dplyr::bind_rows(
        DataTable,
        #base::cbind(
        dplyr::bind_cols(
          Cliffd = res$d,
          PHat = res$phat,
          StdES = res$StdES
        )
      ))
    }

  }
  if (!returnData) {
    #Calculate averages.
    d = dsum / reps
    vard = dvarsum / reps
    sigd = d.sig / reps
    emp.d.var = (d.ss - reps * d ^ 2) / (reps - 1)

    phat = phatsum / reps
    varphat = phatvarsum / reps
    sigphat = phat.sig / reps
    emp.phat.var = (phat.ss - reps * phat ^ 2) / (reps - 1)

    ktau = ksum / reps
    emp.tau.var = (k.ss - reps * ktau ^ 2) / (reps - 1)
    kconsistentvar = kconsistentvar / reps
    kpowerCVt = ksigCVt / reps

    ES = ES / reps
    StdES = StdES / reps
    emp.StdESvar = (StdES.ss - reps * StdES ^ 2) / (reps - 1)
    Var = Var / reps
    tpower = tsig / reps
    MedDiff = MedDiff / reps

    if (type == "l")
    {
      ESLog = ES.l.trans / reps
      StdESLog = StdES.l.trans / reps
      VarLog = Var.l.trans / reps
      Log.sig = trans.sig / reps
    }
  }
  if (!returnData) {
    if (type == "l")
      outcome = tibble::tibble(
        phat,
        varphat,
        sigphat,
        emp.phat.var,
        d,
        vard,
        sigd,
        emp.d.var,
        ktau,
        kconsistentvar,
        emp.tau.var,
        kpowerCVt,
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

    else
      outcome = tibble::tibble(
        phat,
        varphat,
        sigphat,
        emp.phat.var,
        d,
        vard,
        sigd,
        emp.d.var,
        ktau,
        kconsistentvar,
        emp.tau.var,
        kpowerCVt,
        StdES,
        ES,
        Var,
        emp.StdESvar,
        MedDiff,
        tpower
      )
  }
  else {
    outcome = DataTable
  }

  return(outcome)

}

##############################################################
#' @title NP4GroupMetaAnalysisSimulation
#' @description This function simulates data from a family of experiments, where the number of experiments in a family is defined by ther parameter Exp. It simulates data from one of four distributions and uses the data to construct four of groups of equal size (GroupSize). Two groups are assigned as control groups and their distribution is based on the parameter mean and the parameter spread, however, the mean and spread for the control group in Block 2 can be adjusted using the parameters BlockEffect and BlockStdAdj respectively. The other two groups are treatment groups and their distribution is based on the mean+diff and the spread parameter, but the distributions can be adjusted using the StdAdj, BlockEffect and BlockStdAdj parameters. The data from each experiment is analysed separately to estimate the non-parametric statistics P-hat, Cliff's d and Kendall's tau and their variances. The statistics are then meta-analysed using the method specified by the MAMethod parameter. We output the average non-parametric effect statistics across the Exp experimet as if from a single large experiment and also the results of meta-analysising each non-parametric effect size. We use the standard parametric effect sizes and their meta-analysis as baselines. All tests of significance are done at the 0.05 level.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export NP4GroupMetaAnalysisSimulation
#' @param mean The default value used for the group means in the simulated data. It can be any real number including zero.
#' @param sd The default value used for the spread of the control group and the spread of the treatment group in the simulated data. The value must be a real value greater than 0.
#' @param diff mean+diff is the value used for the mean of the treatment group. It can be zero.
#' @param GroupSize is the size of each of the 4 groups comprising one experiment. GroupSize should be an integer of 4 or more
#' @param Exp is the number of experiments being simulated. Exp should be an integer of 2 or more. It defaults to 5.
#' @param type specifies the distribution being simulated. The permitted values are "n" for the normal distribution,  "l" for the lognormal distribution, "g" for the gamma distribution and "lap" for the Laplace dsitribution. The parameter defaults to "n".
#' @param alpha The Type 1 value used in all significance tests. It defaults to 0.05.
#' @param seed specifies the seed to be used to initiate the simulation, so the simulation is repeatable. It defaults to 123.
#' @param BlockEffect is the effect of having two different blocks
#' @param BlockStdAdj is the variance associated with the Block. If BlockStdAdj is zero it means we are treat the block effect as a fixed effect. If BlockStdAdj>0, we treat the block effect as a random effect and increase the variance of Block 2 data.
#' @param StdAdj The value used to introduce heterogeneity into the treatment groups variance if required.
#' @param MAMethod defines the method used for meta-analysis
#' @param StdExp The value used to introduce heterogeneity into experiments in a family required.
#' @param returnES default to FALSE
#' @return If returnES is true, the function returns the summary meta-analysis summary statistics otherwise the function returns the effect sizes for each experiment
#' @examples
#' NP4GroupMetaAnalysisSimulation(mean=0,sd=1,diff=0.5,GroupSize=10,Exp=5,type="n",alpha=0.05,
#' seed=123,StdAdj=0,BlockEffect=0.5,BlockStdAdj=0,StdExp=0,MAMethod="PM")
#'# A tibble: 1 x 30
#'#  NumExp GroupSize AveKtau AveKtauctvar tauSigCVt AveCliffd AveCliffdvar AveCliffdsig Avephat
#'# Avephatvar Avephatsig
#'#   <dbl>     <dbl>   <dbl>        <dbl> <lgl>         <dbl>        <dbl> <lgl>          <dbl>
#'# <dbl> <lgl>
#'# 1      5        10   0.182      0.00188 TRUE          0.346      0.00673 TRUE           0.673
#'# 0.00163 TRUE
#'# … with 19 more variables: MAMean <dbl>, MAvar <dbl>, MASig <lgl>, QE <dbl>, QEp <dbl>,
#'# HetSig <lgl>, P.mean <dbl>,

#' NP4GroupMetaAnalysisSimulation(mean=0,sd=1,diff=0.724,GroupSize=10,Exp=5,type="l",alpha=0.05,
#' seed=123,StdAdj=0,BlockEffect=0.5,BlockStdAdj=0,StdExp=0,MAMethod="PM")
#' # A tibble: 1 x 30
#' #  NumExp GroupSize AveKtau AveKtauctvar tauSigCVt AveCliffd AveCliffdvar AveCliffdsig Avephat
#' # Avephatvar Avephatsig
#' #   <dbl>     <dbl>   <dbl>        <dbl> <lgl>         <dbl>        <dbl> <lgl>          <dbl>
#' #  <dbl> <lgl>
#' # 1      5        10   0.244      0.00167 TRUE          0.464      0.00593 TRUE           0.732
#' # 0.00144 TRUE
#' # … with 19 more variables: MAMean <dbl>, MAvar <dbl>, MASig <lgl>, QE <dbl>, QEp <dbl>,
#' # HetSig <lgl>, P.mean <dbl>,


#' NP4GroupMetaAnalysisSimulation(mean=0,sd=1,diff=0.5,GroupSize=10,Exp=5,type="n",alpha=0.05,
#' seed=123,StdAdj=0,BlockEffect=0.5,BlockStdAdj=0,StdExp=0,MAMethod="PM",returnES=TRUE)
#' # A tibble: 5 x 16
#' #  MeanExp VarExp StdESExp    df  tval   tpval    tciL  tciU Cliffd Cliffdvar  PHat PHatvar
#' # PHatdf     g gvar.approx
#' #    <dbl>  <dbl>    <dbl> <dbl> <dbl>   <dbl>   <dbl> <dbl>  <dbl>     <dbl> <dbl>   <dbl>
#' # <dbl> <dbl>       <dbl>
#' #1   0.940  0.783    1.06   31.3 3.36  0.00206  0.370  1.51   0.58     0.0243 0.29  0.00587
#' # 30.2 1.04       0.112
#' #2   0.372  0.943    0.383  35.0 1.21  0.234   -0.251  0.996  0.21     0.0380 0.105 0.00927
#' # 31.3 0.375      0.0977
#' #3   0.598  0.619    0.761  28.6 2.40  0.0229   0.0892 1.11   0.37     0.0336 0.185 0.00813
#' # 30.8 0.740      0.104
#' #4   0.873  1.13     0.821  28.1 2.60  0.0148   0.184  1.56   0.440    0.0333 0.220 0.00813
#' # 23.8 0.799      0.106
#' #5   0.243  1.03     0.240  31.5 0.758 0.454   -0.410  0.896  0.13     0.0390 0.065 0.00946
#' # 32.8 0.234      0.0961
#' # … with 1 more variable: Cohendvar <dbl>

NP4GroupMetaAnalysisSimulation = function(mean,
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
                                          returnES = FALSE) {
  N = GroupSize

  set.seed(seed)
  ES.cor = rep(NA, Exp) # Used to hold Kendall's tau for each experiment
  ES.ctvar = rep(NA, Exp) # Used to hold the consistent variance of Kendall's tau for each experiment
  ES.d = rep(NA, Exp) # Used to hold Cliff's d for each experiment
  ES.vard = rep(NA, Exp) # Used to hold the variance Cliff's d for each experiment
  ES.phat = rep(NA, Exp)  # Used to hold phat for each experiment
  ES.phatvar = rep(NA, Exp) # Used to hold the variance phat for each experiment
  ES.phat.df = rep(NA, Exp) # Used to hold the degrees of freedom for the test statistics for each experiment

  g = rep(NA, Exp) # Used to hold the small sample size adjusted standardized effect size
  varg.approx = rep(NA, Exp) # Used to hold the approximate variance of the small sample size adjusted standardized effect size
  varg.exact = rep(NA, Exp) # Used to hold the exact variance of the small sample size adjusted standardized effect size
  P.cor = rep(NA, Exp) # This holds the estimate of the point bi-serial correlation effect size
  P.corvar = rep(NA, Exp) # This holds the variance of normal transformation of the correlation effect size
  c = rep(NA, Exp) # This holds the small sample size adjustment factor for each experiment
  df = rep(NA, Exp) # This holds the degrees of freedom of the t test of the unstandardized effect size for each for each experiment

  Cohen.d = rep(NA, Exp) # This holds standardized mean difference effect size of each experiment
  varCohend.approx = rep(NA, Exp)  # This holds the variance of the standardized mean difference effect size of each experiment
  tval = rep(NA, Exp) # This holds t test value for each experiment
  tpval = rep(NA, Exp) # This holds t test p-value for each experiment
  tciL = rep(NA, Exp) # This holds lower bound for the effect size
  tciU = rep(NA, Exp) # This holds upper bound for the effect size
  UES = rep(NA, Exp) # This holds the unstandardized effect size
  Var = rep(NA, Exp) # This holds the pooled within group variance


  # Generate 4-group experiments for a family of experiments (the number in the family defined by the "Exp" parameter), with the underlying distribution defined by the "type" parameter
  for (i in 1:Exp) {
    if (type == "n") {
      if (StdExp == 0)
        ExpAdj = 0
      else
        ExpAdj = stats::rnorm(1, 0, StdExp)
      x = stats::rnorm(N, mean + ExpAdj, sd)
      y = stats::rnorm(N, mean + ExpAdj + diff, sd + StdAdj)
      x2 = stats::rnorm(N, mean + ExpAdj + BlockEffect, sd + BlockStdAdj)
      y2 = stats::rnorm(N,
                        mean + ExpAdj + diff + BlockEffect,
                        sd + BlockStdAdj + StdAdj)
    }
    if (type == "g")
    {
      shape = sd
      rate = mean

      if (StdExp == 0)
        ExpAdj = 0
      else
        ExpAdj = abs(stats::rnorm(1, 0, StdExp))
      if (BlockStdAdj == 0)
        BlockStdAdj = 1

      x = stats::rgamma(N, shape, rate + ExpAdj)
      y = stats::rgamma(N, shape + StdAdj, rate + ExpAdj + diff, )
      x2 = stats::rgamma(N, shape * BlockStdAdj + BlockEffect, rate + ExpAdj)
      y2 = stats::rgamma(N,
                         shape * BlockStdAdj + StdAdj + BlockEffect,
                         rate + ExpAdj + diff)
    }
    if (type == "l")
    {
      if (StdExp == 0)
        ExpAdj = 0
      else
        ExpAdj = stats::rnorm(1, 0, StdExp)
      x = stats::rlnorm(N, mean + ExpAdj, sd)
      y = stats::rlnorm(N, mean + ExpAdj + diff, sd + StdAdj)
      x2 = stats::rlnorm(N, mean + ExpAdj + BlockEffect, sd + BlockStdAdj)
      y2 = stats::rlnorm(N,
                         mean + ExpAdj + diff + BlockEffect,
                         sd + BlockStdAdj + StdAdj)
    }

    if (type == "lap")
      # Changed to use built-in defaults for max and min
    {
      if (StdExp == 0)
        ExpAdj = 0
      else
        ExpAdj = stats::rnorm(1, 0, StdExp)
      x = LaplaceDist(N, mean + ExpAdj, sd + StdAdj)
      y = LaplaceDist(N, mean + diff + ExpAdj, sd)
      x2 = LaplaceDist(N, mean + ExpAdj + BlockEffect, sd + BlockStdAdj)
      y2 = LaplaceDist(N,
                       mean + ExpAdj + diff + BlockEffect,
                       sd + BlockStdAdj + StdAdj)
    }


    # The data in each experiment is analysed to obtain the values of three non-parametric effect size: Kendall's tau, its consistent variance, Cliff's d with its consistent variance and phat with its variance plus t-test results based on the Welch method applied to a randomized blocks experiment.

    NPStats = NULL

    # Calc4GroupNPStats finds the average value of various non-parametric statistics for a randomised block experimental design
    NPStats = Calc4GroupNPStats(y, x, y2, x2, alpha = alpha)
    ES.cor[i] = as.numeric(NPStats$cor)
    ES.ctvar[i] = NPStats$ctvar
    ES.d[i] = NPStats$d
    ES.vard[i] = NPStats$vard
    ES.phat[i] = NPStats$phat - 0.5 #metafor tests difference from zero but zero effect for phat is 0.5
    ES.phatvar[i] = NPStats$phat.var
    ES.phat.df[i] = NPStats$phat.df # Holds the degrees of freedom

    # Do a parametric analysis of the data
    UES[i] = (base::mean(y) + base::mean(y2) - base::mean(x) - base::mean(x2)) /
      2
    Var[i] = (stats::var(y) + stats::var(y2) + stats::var(x) + stats::var(x2)) /
      4

    newlist = list()
    newlist[[1]] = x
    newlist[[2]] = y
    newlist[[3]] = x2
    newlist[[4]] = y2
    vec = c(-1, 1, -1, 1) / 2

    # RandomizedBlocksAnalysis does a Welch-based linear contrast analysis
    res.t = RandomizedBlocksAnalysis(newlist, con = vec, alpha = 0.05)

    ttest.ESdets = base::data.frame(res.t$psihat)

    # Cohen's d is the standardized mean difference
    Cohen.d[i] = UES[i] / sqrt(Var[i])

    # g is 	small sample size adjusted mean difference

    tdets = base::as.data.frame(res.t$test)
    df[i] = as.numeric(tdets$df)
    c[i] = reproducer::calculateSmallSampleSizeAdjustment(df[i])
    g[i] = Cohen.d[i] * c[i]
    tval[i] = tdets$test
    tpval[i] = ttest.ESdets$p.value
    tciL[i] = ttest.ESdets$ci.lower
    tciU[i] = ttest.ESdets$ci.upper


    varg.approx[i] = c[i] ^ 2 / N + g[i] ^ 2 / (2 * df[i])

    varCohend.approx[i] = 1 / N + Cohen.d[i] ^ 2 / (2 * df[i])

    # Calculate r from d

    P.cor[i] = Cohen.d[i] / sqrt(Cohen.d[i] ^ 2 + 4)
    P.cor[i] = reproducer::transformRtoZr(P.cor[i])
    P.corvar[i] = 1 / (4 * N - 3)

  }


  if (!returnES)
  {
    NumExp = Exp
    # Analyse the data from the Exp experiments as a single distributed experiment
    AveKtau = base::mean(ES.cor)
    AveKtauctvar = sum(ES.ctvar) / NumExp ^ 2
    # Use t-based bounds
    vv = stats::qt(alpha / 2, (4 * N - 3) * NumExp)
    LowerBoundt = AveKtau + vv * sqrt(AveKtauctvar)
    UpperBoundt = AveKtau - vv * sqrt(AveKtauctvar)

    tauSigCVt = UpperBoundt < 0 | LowerBoundt > 0


    AveCliffd = base::mean(ES.d)
    AveCliffdvar = sum(ES.vard) / NumExp ^ 2
    # The degrees of freedom per experiment is the number of participants less 3.
    vv = stats::qt(alpha / 2, (4 * N - 3) * NumExp)
    LowerBoundd = AveCliffd + vv * sqrt(AveCliffdvar)
    UpperBoundd = AveCliffd - vv * sqrt(AveCliffdvar)
    AveCliffdsig = UpperBoundd < 0 | LowerBoundd > 0

    Avephat = base::mean(ES.phat)
    Avephatvar = sum(ES.phatvar) / NumExp ^ 2
    phatdf = sum(ES.phat.df)
    vv = stats::qt(alpha / 2, phatdf)
    LowerBoundphat = Avephat + vv * sqrt(Avephatvar)
    UpperBoundphat = Avephat - vv * sqrt(Avephatvar)
    Avephatsig = UpperBoundphat < 0 | LowerBoundphat > 0
    # Add back the 0.5
    Avephat = Avephat + 0.5



    # Perform a standard meta-analysis of the ktau values for all simulated experiments

    meth = MAMethod

    MA.res = metafor::rma(ES.cor, ES.ctvar, method = meth)
    # Extract the data from the meta-analysis

    MAMean = as.numeric(MA.res$beta)
    MAvar = as.numeric(MA.res$se ^ 2)
    QE = as.numeric(MA.res$QE)
    QEp = as.numeric(MA.res$QEp)
    UB = as.numeric(MA.res$ci.ub)
    LB = as.numeric(MA.res$ci.lb)
    MASig = UB < 0 | LB > 0
    HetSig = QEp < 0.05

    #	Meta-analyse Cliff'd d result

    MA.dres = metafor::rma(ES.d, ES.vard, method = meth)
    pvalue.d = as.numeric(MA.dres$pval)
    d.sig = pvalue.d < 0.05
    Mean.d = as.numeric(MA.dres$beta)

    #	Meta-analysis of phat

    MA.phat = metafor::rma(ES.phat, ES.phatvar, method = meth)
    Mean.phat = as.numeric(MA.phat$beta) + 0.5 # Add the 0.5 back
    phat.sig = as.numeric(MA.phat$pval) < 0.05

    # Meta-analysis of g values for exact and approx.

    # Note when N is the same for each study the variance for the study based on the average of d is the same for each study. The variance of g is affected by the degrees of freedom which based on a Welch style analysis may be less than a simple ANOVA
    for (i in 1:Exp)
      varg.exact[i] = varStandardizedEffectSize(mean(Cohen.d[i]), 4 / N, df[i], returnVarg =
                                                  TRUE)
    MA.g.exact = metafor::rma(g, varg.exact, method = meth)
    Mean.g.exact = as.numeric(MA.g.exact$beta)
    g.exact.sig = as.numeric(MA.g.exact$pval) < 0.05

    MA.g.approx = metafor::rma(g, varg.approx, method = meth)
    Mean.g.approx = as.numeric(MA.g.approx$beta)
    g.approx.sig = as.numeric(MA.g.approx$pval) < 0.05


    # Do a Pearson correlation analysisn

    MAP.res = metafor::rma(P.cor, P.corvar, method = meth)
    MAPresults = ExtractMAStatistics(MAP.res, N, N, type = "r")
    P.mean = MAPresults$mean
    P.rsig = MAPresults$pvalue <= 0.05
    P.hetsig = MAPresults$QEp <= 0.05

    # Meta-analysis of Cohen's d

    Cohend.res = metafor::rma(Cohen.d, varCohend.approx)
    Cohend.mean = as.numeric(Cohend.res$beta)
    Cohend.sig = as.numeric(Cohend.res$pval) < 0.05

    output = tibble::tibble(
      NumExp,
      GroupSize,
      AveKtau,
      AveKtauctvar,
      tauSigCVt,
      AveCliffd,
      AveCliffdvar,
      AveCliffdsig,
      Avephat,
      Avephatvar,
      Avephatsig,
      MAMean,
      MAvar,
      MASig,
      QE,
      QEp,
      HetSig,
      P.mean,
      P.rsig,
      P.hetsig,
      Mean.phat,
      phat.sig,
      Mean.d,
      d.sig,
      Mean.g.exact,
      g.exact.sig,
      Mean.g.approx,
      g.approx.sig,
      Cohend.mean,
      Cohend.sig
    )
  }
  else{
    output = tibble::tibble(
      dplyr::bind_cols(
        MeanExp = UES,
        VarExp = Var,
        StdESExp = Cohen.d,
        df = df,
        tval = tval,
        tpval = tpval,
        tciL = tciL,
        tciU = tciU,
        Cliffd = ES.d,
        Cliffdvar = ES.vard,
        PHat = ES.phat,
        PHatvar = ES.phatvar,
        PHatdf = ES.phat.df,
        g = g,
        gvar.approx = varg.approx,
        Cohendvar = varCohend.approx
      )
    )

  }
  return(output)

}

#########################################################
#' @title NP2GroupMetaAnalysisSimulation
#' @description This function simulates data from a family of experiments. The parameter Exp deteremines the number of experiments in the family. The function simulates data from one of four distributions and uses the data to construct two of groups of equal size (GroupSize). The distribution for one  of the groups corresponds to the control and is based on the given mean and spread, the distribution for the other group corresponds to the treatment group and  is based on the mean+diff and the spread plus any variance adjusemtn (determined by the parametrt StdAdj). The data from each experiment is analysed separately to estimate three non-parametric effect sizes: the point bi-serial version of Kendall's tau, Cliff's d and the probability of superiority referred to as phat and their variances. Parametric effect sizes Cohen's d (also known as the standarized means difference, SMD) and the small sample size adjusted standardized mean difference g are also calculated to gether with their variances. The effect sizes are then meta-analysed using two main methods: the simple avarge of the effect size and the variance weighted average. The function uses the metafor package for formal meta-analysis, and the specific method of formal meta-analysis used is determined by the MAAMethod. All tests of signficance are done at the 0.05 level. If the parameter returnES is TRUE, the function returns it returns the effect sizes for each experiment in the family, otherwise it returns the meta-analysis results.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export NP2GroupMetaAnalysisSimulation
#' @param mean the value used for the mean of control group in the simulated data. It can be any real number including zero.
#' @param sd the value used for the spread of the control group and the spread of the treatment group in the simulated data. The value must be a real value greater than 0.
#' @param diff mean+diff is the value used for the mean of the treatment group. It can be zero.
#' @param GroupSize is the size of each of the 2 groups comprising one experiment. GroupSize should be an integer of 4 or more
#' @param Exp is the number of experiments being simulated. Exp should be an integer of 2 or more. It defaults to 5.
#' @param type specifies the distribution being simulated. The permited values are "n" for the normal distribution,  "l" for the lognormal distribution, "g" for the gamma distribution and "lap" for the Laplace distribution. The parameter defaults to "n".
#' @param StdAdj specifies a level used to adjust the treatment variance. It allows heterogeneity to be modelled. It defaults to zero meaning no variance heterogeneity is introduced.
#' @param seed specifies the seed to be used to initiate the simulation, so the simulation is repeatable. It defauls to 123.
#' @param StdExp specifies the adjustment made to the variance of each experiment, if it is required to simulate variance heterogeneity due to individual experiments in a family. It defaults to 0.
#' @param alpha the Type 1 error rate level use for statistical tests itb defaults to 0.05.
#' @param MAMethod the meta-analysis method needed for the call to the metafor package rma algorithm
#' @param returnES Determines the format of the output. It defaults to FALSE which causes the function to output the meta-anaysis results for the family of experiments. If set to TRUE it returns the effect sizes for each experiment.
#' @return Depending on the value of the returnES parameter, the function either returnd the effect sizes for each experiment or the aggregated resutls for the family
#' @examples
#' NP2GroupMetaAnalysisSimulation(mean=0,sd=1,diff=0.5,GroupSize=100,Exp=5,type="n",StdAdj=0,
#' alpha=0.05,seed=123,StdExp=1,MAMethod="PM",returnES=FALSE)
#' # A tibble: 1 x 30
#' #  NumExp GroupSize AveKtau AveKtauctvar tauSigCVt AveCliffd AveCliffdvar AveCliffdsig Avephat
#' #  Avephatvar Avephatsig
#' #   <dbl>     <dbl>   <dbl>        <dbl> <lgl>         <dbl>        <dbl> <lgl>          <dbl>
#' #  <dbl> <lgl>
#' # 1      5       100   0.126     0.000315 TRUE          0.250      0.00125 TRUE           0.625
#' # 0.000310 TRUE
#' # … with 19 more variables: MAMean <dbl>, MAvar <dbl>, MASig <lgl>, QE <dbl>, QEp <dbl>,
#' # HetSig <lgl>, P.mean <dbl>,

#' NP2GroupMetaAnalysisSimulation(mean=0,sd=1,diff=0.5,GroupSize=10,Exp=5,type="n",StdAdj=0,
#' alpha=0.05,seed=123,StdExp=1,MAMethod="PM",returnES=TRUE)
#' # A tibble: 5 x 17
#' # MeanExp VarExp StdESExp    df    tval  tpval   tciL   tciU Cliffd Cliffdvar  PHat PHatvar
#' # PHatdf       g gvar.exact
#' #    <dbl>  <dbl>    <dbl> <dbl>   <dbl>  <dbl>  <dbl>  <dbl>  <dbl>     <dbl> <dbl>   <dbl>
#' # <dbl>   <dbl> <lgl>
#' #1  0.226   1.03    0.223   17.9 -0.498  0.624  -1.18   0.728   0.1     0.0808  0.55 0.0197
#' # 17.8  0.213  NA
#' #2  1.00    0.622   1.27    15.4 -2.84   0.0122 -1.75  -0.251   0.64    0.0400  0.82 0.00947
#' # 15.0  1.21   NA
#' #3  0.432   0.848   0.469   18.0 -1.05   0.308  -1.30   0.434   0.36    0.0725  0.68 0.0177
#' # 15.5  0.449  NA
#' #4  0.434   0.967   0.441   13.7 -0.986  0.341  -1.38   0.512   0.22    0.0759  0.61 0.0184
#' # 14.6  0.416  NA
#' #5 -0.0342  0.782  -0.0387  14.9  0.0865 0.932  -0.809  0.878  -0.2     0.0817  0.4  0.02
#' # 17.3 -0.0367 NA
#' # … with 2 more variables: gvar.approx <dbl>, ...17 <dbl>

#' NP2GroupMetaAnalysisSimulation(mean=0,sd=1,diff=0.724,GroupSize=10,Exp=5,type="l",StdAdj=0,
#' alpha=0.05,seed=123,StdExp=1,MAMethod="PM",returnES=FALSE)
#' # A tibble: 1 x 30
#' #  NumExp GroupSize AveKtau AveKtauctvar tauSigCVt AveCliffd AveCliffdvar AveCliffdsig Avephat
#' # Avephatvar Avephatsig
#' #   <dbl>     <dbl>   <dbl>        <dbl> <lgl>         <dbl>        <dbl> <lgl>          <dbl>
#' # <dbl> <lgl>
#' # 1      5        10   0.181      0.00360 TRUE          0.344       0.0129 TRUE           0.672
#' # 0.00312 TRUE
#' # … with 19 more variables: MAMean <dbl>, MAvar <dbl>, MASig <lgl>, QE <dbl>, QEp <dbl>,
#' # HetSig <lgl>, P.mean <dbl>,
#' #   P.rsig <lgl>, P.hetsig <lgl>, Mean.phat <dbl>, phat.sig <lgl>, Mean.d <dbl>, d.sig <lgl>,
#' # Mean.g.exact <dbl>,
#' #   g.exact.sig <lgl>, Mean.g.approx <dbl>, g.approx.sig <lgl>, Cohend.mean <dbl>,
#' # Cohend.sig <lgl>

NP2GroupMetaAnalysisSimulation = function(mean,
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
                                          returnES = FALSE) {
  base::set.seed(seed)
  N = GroupSize

  ES.cor = rep(NA, Exp) # This will hold the Kendall's tau for each experiment
  ES.var = rep(NA, Exp) #This will hold the estimated consistent variance of tau for each experiment
  P.cor = rep(NA, Exp) # This holds the Pearson correlation
  P.corvar = rep(NA, Exp)#This holds the variance of normal transformation of the Pearson correlation

  VarExp = rep(NA, Exp) # This holds the pooled variance for each experiment
  MeanExp = rep(NA, Exp) # This holds the meandifference for each experiment
  StdESExp = rep(NA, Exp) # This holds the standardized mean difference

  Cliffd = rep(NA, Exp) #This holds Cliff's d
  Cliffdvar = rep(NA, Exp) # This holds the variance of Clff's d

  PHat = rep(NA, Exp) #This holds the probability of superiority phat
  PHatvar = rep(NA, Exp) #This holds the variance of phat
  PHatdf = rep(NA) # This holds the degrees of freedom of the t-test for phat


  g = rep(NA, Exp) # This holds the small sample size adjusted standardized mean difference
  gvar.exact = rep(NA, Exp) # This holds the exact variance of the small sample size adjusted
  # standardized mean difference
  gvar.approx = rep(NA, Exp) # This holds the approximate variance of the small sample size adjusted
  # standardized mean difference
  df = rep(NA, Exp) # This holds the degrees of freedom of the t-test of the mean difference
  c = rep(NA, Exp) # This holds the small sample size adjustment factor for the standardized mean
  # difference
  tval = rep(NA, Exp) # This holds the t-test value for each experiment
  tpval = rep(NA, Exp)  # This holds the p-value of the t-test for each experiment
  tciL = rep(NA, Exp) # This holds the lower confidence interval bound of the standrdized mean difference for each experiment
  tciU = rep(NA, Exp)  # This holds the upper confidence interval bound of the standrdized mean difference for each experiment
  Cohendvar = rep(NA, Exp) # This holds the variance of the standardized mean difference for each
  # experiment

  # For the point-biserial tau we need a dummy variable that takes the value zero for control group
  # data points and 1 for treatment group data points
  dummy = c(rep(0, N), rep(1, N))

  for (i in 1:Exp)
  {
    if (type == "n") {
      if (StdExp == 0)
        ExpAdj = 0
      else
        ExpAdj = stats::rnorm(1, 0, StdExp)
      x = stats::rnorm(N, mean + ExpAdj, sd)
      y = stats::rnorm(N, mean + ExpAdj + diff, sd + StdAdj)
    }
    if (type == "g")
    {
      shape = sd
      rate = mean
      if (StdExp == 0)
        ExpAdj = 0
      else
        ExpAdj = abs(stats::rnorm(1, 0, StdExp))
      x = stats::rgamma(N, shape, rate + ExpAdj)
      y = stats::rgamma(N, shape + StdAdj, rate + ExpAdj + diff)

    }
    if (type == "l")
    {
      if (StdExp == 0)
        ExpAdj = 0
      else
        ExpAdj = stats::rnorm(1, 0, StdExp)
      x = stats::rlnorm(N, mean + ExpAdj, sd)
      y = stats::rlnorm(N, mean + ExpAdj + diff, sd + StdAdj)
    }

    if (type == "lap")
      # Changed to use built-in defaults for max and min
    {
      if (StdExp == 0)
        ExpAdj = 0
      else
        ExpAdj = stats::rnorm(1, 0, StdExp)
      x = LaplaceDist(N, mean + ExpAdj, sd)
      y = LaplaceDist(N, mean + diff + ExpAdj, sd + StdAdj)
    }


    # This holds the simulated data set with each observation assigned to a group in an experiment

    # The data in each experiment is analysed to obtain the value of Kendall's tau, its consistent variance and its hypothesis testing variance

    xy = c(x, y)


    expdata = base::data.frame(xy, dummy)

    ktau = Kendalltaupb(expdata$xy, expdata$dummy)
    ES.cor[i] = ktau$cor
    ES.var[i] = ktau$consistentvar

    # Analyse the generated data for each member of the family using the cid function to obtain
    # Cliff's d and its variance
    Cliff = Cliffd(y, x)
    Cliffd[i] = Cliff$d
    Cliffdvar[i] = Cliff$sqse.d

    # Analyse the generated data for each member of the family using the bmp function to obtain phat
    # and its variance

    PHat.res = calculatePhat(x, y)
    PHat[i] = PHat.res$phat
    PHatvar[i] = PHat.res$s.e. ^ 2
    PHatdf[i] = PHat.res$df

    # Prepare to do a standard analysis for comparison
    # We can do both a standardized effect size or correlation. For the paper we use the
    # standardized mean difference effect size

    VarExp[i] = (stats::var(x) + stats::var(y)) / 2
    MeanExp[i] = base::mean(y) - base::mean(x)
    StdESExp[i] = (MeanExp[i]) / sqrt(VarExp[i])



    tempttest = stats::t.test(x, y)
    df[i] = as.numeric(tempttest$parameter)
    tval[i] = tempttest$statistic
    tpval[i] = tempttest$p.value
    tciL[i] = tempttest$conf.int[1]
    tciU[i] = tempttest$conf.int[2]
    c[i] = reproducer::calculateSmallSampleSizeAdjustment(df[i])
    g[i] = StdESExp[i] * c[i]

    # N in each group

    gvar.approx[i] = 2 * c[i] ^ 2 / N + g[i] ^ 2 / (2 * df[i])

    Cohendvar[i] = 2 / N + StdESExp[i] ^ 2 / (2 * df[i])

    Pearsonr = g[i] / sqrt(g[i] ^ 2 + 4)
    TempP = as.numeric(Pearsonr)
    TempP = reproducer::transformRtoZr(TempP)
    P.cor[i] = TempP
    P.corvar[i] = 1 / (2 * N - 3)

  }

  if (!returnES)
  {
    NumExp = Exp

    # First aggregate as if the family is a planned distributed experiment, i.e. the "experiment"
    # is a blocking factor

    # Aggregate tau_pb values


    AveKtau = base::mean(ES.cor)
    AveKtauctvar = sum(ES.var) / NumExp ^ 2

    #Obtain confidence intervals using t-distribtion with 2*N-3*NumExp degrees of freedom

    vv = stats::qt(alpha / 2, (2 * N - 3) * NumExp)
    LowerBoundt = AveKtau + vv * sqrt(AveKtauctvar)
    UpperBoundt = AveKtau - vv * sqrt(AveKtauctvar)

    tauSigCVt = UpperBoundt < 0 | LowerBoundt > 0



    #Aggregate the Cliff's d values

    AveCliffd = base::mean(Cliffd)
    AveCliffdvar = sum(Cliffdvar) / NumExp ^ 2
    vvd = stats::qt(alpha / 2, (2 * N - 2) * NumExp)
    LowerBoundd = AveCliffd + vvd * sqrt(AveCliffdvar)
    UpperBoundd = AveCliffd - vvd * sqrt(AveCliffdvar)

    AveCliffdsig = UpperBoundd < 0 | LowerBoundd > 0

    #Aggregate the phat values

    Avephat = base::mean(PHat)
    Avephatadj = Avephat - 0.5 # The null effect for phat is 0.5
    Avephatvar = sum(PHatvar) / NumExp ^ 2
    pdf = sum(PHatdf)
    vvp = stats::qt(alpha / 2, pdf * NumExp)
    # Test whether Avephatadj is different from zero, tests whether phat is different from 0.5
    LowerBoundp = Avephatadj + vvp * sqrt(Avephatvar)
    UpperBoundp = Avephatadj - vv * sqrt(Avephatvar)

    Avephatsig = UpperBoundp < 0 | LowerBoundp > 0

    # Perform a standard meta-analysis of the non-parametric statistics for the family of simulated
    # experiments


    # tau_pb meta-analysis


    meth = MAMethod

    MA.res = metafor::rma(ES.cor, ES.var, method = meth)
    # Extract the data from the analysis

    MAtauResults = ExtractMAStatistics(MA.res, 2 * N, 2 * N, Transform = FALSE)
    MAMean = MAtauResults$mean
    UB = MAtauResults$UB
    LB = MAtauResults$LB
    MASig = UB < 0 | LB > 0
    MAvar = as.numeric(MA.res$se ^ 2)
    QE = MAtauResults$QE
    QEp = MAtauResults$QEp
    HetSig = QEp < 0.05

    # Cliff's d meta-analysis
    MAd.res = metafor::rma(Cliffd, Cliffdvar, method = meth)
    MAd.Results = ExtractMAStatistics(MAd.res, 2 * N, 2 * N, Transform = FALSE)
    Mean.d = MAd.Results$mean
    d.sig = MAd.Results$pvalue < 0.05

    # Probability of Superiority phat meta-analysis
    PHatAdj = PHat - 0.5
    MAphat.res = metafor::rma(PHatAdj, PHatvar, method = meth) # meta-analyse phat-0.5
    MAphat.Results = ExtractMAStatistics(MAphat.res, 2 * N, 2 * N, Transform =
                                           FALSE)
    Mean.phat = MAphat.Results$mean + 0.5 # Add back the 0.5 value
    phat.sig = MAphat.Results$pvalue < 0.05


    # Do a Pearson correlation analysis as a comparison with tau_pb

    MAP.res = metafor::rma(P.cor, P.corvar, method = meth)

    MAPresults = ExtractMAStatistics(MAP.res, 2 * N, 2 * N, type = "r")
    P.mean = MAPresults$mean
    P.rsig = MAPresults$pvalue <= 0.05
    P.hetsig = MAPresults$QEp <= 0.05

    #Do an analysis of the small sample size standardized effect size as a comparison with Cliff's d
    MAgres.approx = metafor::rma(g, gvar.approx, method = meth)
    MAgresults.approx = ExtractMAStatistics(MAgres.approx, 2 * N, 2 * N, Transform =
                                              FALSE)
    Mean.g.approx =	MAgresults.approx$mean
    g.approx.sig = MAgresults.approx$pvalue < 0.05


    # Note when N is the same for each study, the variance for the study based on the average of d
    # is the same for each study. The variance of g is affected by the degrees of freedom which
    # based on a Welch style analysis may be less than a simple ANOVA

    for (i in 1:Exp)
      gvar.exact[i] = varStandardizedEffectSize(StdESExp[i], 2 / N, df[i], returnVarg =
                                                  TRUE)
    MA.g.exact = metafor::rma(g, gvar.exact, method = meth)
    Mean.g.exact = as.numeric(MA.g.exact$beta)
    g.exact.sig = as.numeric(MA.g.exact$pval) < 0.05

    # Cohen's d metaanalysis

    Cohend.res = metafor::rma(StdESExp, Cohendvar)
    Cohend.mean = as.numeric(Cohend.res$beta)
    Cohend.sig = as.numeric(Cohend.res$pval) < 0.05

    output = tibble::tibble(
      NumExp,
      GroupSize,
      AveKtau,
      AveKtauctvar,
      tauSigCVt,
      AveCliffd,
      AveCliffdvar,
      AveCliffdsig,
      Avephat,
      Avephatvar,
      Avephatsig,
      MAMean,
      MAvar,
      MASig,
      QE,
      QEp,
      HetSig,
      P.mean,
      P.rsig,
      P.hetsig,
      Mean.phat,
      phat.sig,
      Mean.d,
      d.sig,
      Mean.g.exact,
      g.exact.sig,
      Mean.g.approx,
      g.approx.sig,
      Cohend.mean,
      Cohend.sig
    )
  }

  else  {
    output = tibble::as_tibble(
      dplyr::bind_cols(
        MeanExp = MeanExp,
        VarExp = VarExp,
        StdESExp = StdESExp,
        df = df,
        tval = tval,
        tpval = tpval,
        tciL = tciL,
        tciU = tciU,
        Cliffd = Cliffd,
        Cliffdvar = Cliffdvar,
        PHat = PHat,
        PHatvar = PHatvar,
        PHatdf = PHatdf,
        g = g,
        gvar.exact = gvar.exact,
        gvar.approx = gvar.approx,
        Cohendvar
      )
    )
  }


  return(output)
}


#######################################################
#' @title MetaAnalysisSimulations
#' @description This function simulates data from many families of experiments. The number of families simulated is defined by the Replications parameter. The parameter Exp determines the number of experiments in each family. The function simulates data from one of four distributions and uses the data to construct two of groups of equal size (GroupSize). The experimental design of individual experiments in each family is determined by the FourGroup parameter. If FourGroup=FALSE, the basic experimental design is a balanced two group randomized experiment, otherwise the experimental design is a balanced four group experiment corresponding to a randomized blocks experiment. The function calls either NP2GroupMetaAnalysisSimulation or NP4GroupMetaAnalysisSimulation to generate and analyse data for each individual family. The function either returns the meta-analysed data from each experiment or provides summary statistics.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export MetaAnalysisSimulations
#' @param mean the value used for the mean of control group in the simulated data. It can be any real number including zero.
#' @param sd the value used for the spread of the control group and the spread of the treatment group in the simulated data. The value must be a real value greater than 0.
#' @param diff mean+diff is the value used for the mean of the treatment group. It can be zero.
#' @param GroupSize is the size of each of the groups comprising one experiment. GroupSize should be an integer of 4 or more
#' @param Exp the number of experiments being simulated. Exp should be an integer of 2 or more. It defaults to 5.
#' @param Replications The number of times the set of experiments is simulated. It defaults to 500.
#' @param type specifies the distribution being simulated. The permited values are "n" for the normal distribution,  "l" for the lognormal distribution, "g" for the gamma distribution and "lap" for the Laplace distribution. The parameter defaults to "n".
#' @param seed specifies the seed to be used to initiate the simulation, so the simulation is repeatable. It defaults to 456.
#' @param alpha This is the Type 1 probability level used in significnace tests. It defaults to 0.05.
#' @param FourGroup is a Boolean variable that determines whether the experiment is a two group experiments or a 4-Group randomised block experiment. It defaults to FALSE which means a two group experiment is the default condition.
#' @param StdAdj specifies the adjustment made to the treatment variance if it is required to simulate variance heterogeneity due to the treatment. It defaults to 0.
#' @param BlockEffect specified a fixed block effect that can be applied to randomized blocks designs.
#' @param BlockStdAdj specifies the adjustment made to the variance of one block if it is required to simulate heterogeneity due to  blocks in a randomized blocks experiment. The variance adjusment is made to the same block that any Block affect is applied to. It defaults to 0.
#' @param StdExp specifies the adjustment made to the variance of each experiment  if it is required to simulate variance heterogeneity due to individual experiments in a family. It defaults to 0.
#' @param MAMethod specifies the model to be used when experimental effect sizes ar aggregated using the R metafor package.
#' @param returnES if TRUE the function outputs the summary statistics otherwise it outputs the meta-analysis results for each family. The parameter defaults to FALSE
#' @return The parameter either returns the meta-analysis values obtained from each family or the average values of the meta-analysis over all replications.
#' @examples
#' MetaAnalysisSimulations(mean=0,sd=1,diff=0.5,GroupSize=10,type="n",Replications=25,Exp=5,
#'	seed=456,alpha=0.05,FourGroup=FALSE,StdAdj=0,BlockEffect=0,BlockStdAdj=0,StdExp=0,MAMethod="PM")
#'# A tibble: 1 x 26
#'# Averagektau Averagektauctvar AveragetauSigCVt AverageCliffd AverageCliffdvar AverageCliffdsig
#'# <dbl>            <dbl>            <dbl>         <dbl>            <dbl>            <dbl>
#'#  1       0.164          0.00375             0.72         0.311           0.0134             0.72
#'#… with 20 more variables: Averagephat <dbl>, Averagephatvar <dbl>, Averagephatsig <dbl>,
#'#  MAAverage <dbl>, MAVariance <dbl>, MASignificant <dbl>, HetSignificant <dbl>, PMAAverage <dbl>,
#'#  PMASignificant <dbl>, PMAHetSignificant <dbl>, Mean.phat <dbl>, phat.sig <dbl>, Mean.d <dbl>,
#'#  d.sig <dbl>, Mean.g.exact <dbl>, g.exact.sig <dbl>, Mean.g.approx <dbl>, g.approx.sig <dbl>,
#'#  Mean.Cohend <dbl>, Cohend.sig <dbl>

#'	MetaAnalysisSimulations(mean=0,sd=1,diff=0.734,GroupSize=100,type="l",Replications=25,Exp=5,
#'	seed=456,alpha=0.05,FourGroup=FALSE,StdAdj=0,BlockEffect=0,BlockStdAdj=0,StdExp=0,MAMethod="PM")
#'#  A tibble: 1 x 26
#'#  Averagektau Averagektauctvar AveragetauSigCVt AverageCliffd AverageCliffdvar AverageCliffdsig
#'# <dbl>            <dbl>            <dbl>         <dbl>            <dbl>            <dbl>
#'#  0.200         0.000276                1         0.399          0.00109                1
#'# … with 20 more variables: Averagephat <dbl>, Averagephatvar <dbl>, Averagephatsig <dbl>,
#'#   MAAverage <dbl>, MAVariance <dbl>, MASignificant <dbl>, HetSignificant <dbl>, PMAAverage <dbl>,
#'#   PMASignificant <dbl>, PMAHetSignificant <dbl>, Mean.phat <dbl>, phat.sig <dbl>, Mean.d <dbl>,
#'#   d.sig <dbl>, Mean.g.exact <dbl>, g.exact.sig <dbl>, Mean.g.approx <dbl>, g.approx.sig <dbl>,
#'#   Mean.Cohend <dbl>, Cohend.sig <dbl>

#'# MetaAnalysisSimulations(mean=0,sd=1,diff=0.5,GroupSize=10,type="n",Replications=10,Exp=5,
#'# seed=456,alpha=0.05,FourGroup=TRUE,StdAdj=0,BlockEffect=0.5,BlockStdAdj=0,StdExp=0,
#'# MAMethod="PM",returnES=TRUE)
#'# A tibble: 10 x 31
#'# Family NumExp GroupSize AveKtau AveKtauctvar tauSigCVt AveCliffd AveCliffdvar AveCliffdsig
#'# <int>  <dbl>     <dbl>   <dbl>        <dbl> <lgl>         <dbl>        <dbl> <lgl>
#'#  1      1      5        10   0.133      0.00206 TRUE          0.252      0.00742 TRUE
#'# 2      2      5        10   0.161      0.00190 TRUE          0.306      0.00681 TRUE
#'# 3      3      5        10   0.135      0.00206 TRUE          0.256      0.00744 TRUE
#'# 4      4      5        10   0.16       0.00191 TRUE          0.304      0.00686 TRUE
#'# 5      5      5        10   0.133      0.00178 TRUE          0.252      0.00636 TRUE
#'# 6      6      5        10   0.192      0.00178 TRUE          0.364      0.00634 TRUE
#'# 7      7      5        10   0.189      0.00185 TRUE          0.36       0.00661 TRUE
#'# 8      8      5        10   0.126      0.00199 TRUE          0.24       0.00715 TRUE
#'# 9      9      5        10   0.132      0.00200 TRUE          0.25       0.00719 TRUE
#'# 10     10      5        10   0.18       0.00181 TRUE          0.342      0.00646 TRUE
#'# # … with 22 more variables: Avephat <dbl>, Avephatvar <dbl>, Avephatsig <lgl>, MAMean <dbl>,
#'# #   MAvar <dbl>, MASig <lgl>, QE <dbl>, QEp <dbl>, HetSig <lgl>, P.mean <dbl>, P.rsig <lgl>,
#'# #   P.hetsig <lgl>, Mean.phat <dbl>, phat.sig <lgl>, Mean.d <dbl>, d.sig <lgl>,
#'# #   Mean.g.exact <dbl>, g.exact.sig <lgl>, Mean.g.approx <dbl>, g.approx.sig <lgl>,
#'# #   Cohend.mean <dbl>, Cohend.sig <lgl>

MetaAnalysisSimulations = function(mean = 0,
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
                                   returnES = F) {
  # Set up variables to hold the outcome of the simulations

  Averagektau = 0 # holds the average of the estimate of ktau across all the replications
  Averagektauctvar = 0 # holds the average of the estimation of the consistent variance ktau across all the replications
  AveragetauSigCVt = 0 # Holds the average number of times the difference between the means was assessed as signifciant using the confidence limits obtained from the consistent variance with the t-distribution
  AverageCliffd = 0 # Holds the average of Cliff's d
  AverageCliffdvar = 0 # Holds the average of the variance of Cliff's d
  AverageCliffdsig = 0 # Holds the average number of times the Cliff's d was assessed as signifciant using average hypothesis testing variance
  Averagephat = 0 # Holds the average of phat
  Averagephatvar = 0 # Holds the average of the variance of phat
  Averagephatsig = 0 # Holds the average number of times phat was assessed as signifciant using average hypothesis testing variance


  MAAverage = 0 # holds the average of the estimate of ktau obtained from the meta-analysis across all the replications
  MAVariance = 0 # holds the average of the estimate of the variance of ktau obtained from the meta-analysis across all the replications
  MASignificant = 0 # Holds the average number of times the meta-analysis estimate of ktau was significant
  HetSignificant = 0 # Holds the average number of times the meta-analysis heteogenity test was significant
  MAPcor.mean = 0 # Holds the sum of the mean values of the point biserial correlation using Pearson's formula
  MAPcor.sig = 0 # Holds the number of times the mean Pr is signficiany
  MAPcor.hetsig = 0 # Holds the number of times the meta-analysis detected heterogeneity
  Mean.g.exact = 0 # Holds the standardized effect size based on meta-analysis with the exact variance
  g.exact.sig = 0 # Holds the number of times the standardized effect size was significant with a meta-analysis using the exact variance
  Mean.g.approx = 0 # Holds the standardized effect size based on meta-analysis with the approximate variance
  g.approx.sig = 0 # Holds the number of times the standardized effect size was significant with a meta-analysis using the approximate variance
  Mean.phat = 0 # Holds the mean phat effect size  based on a meta-analysis

  phat.sig = 0 # Holds the number of times the phat effect size was significant with a meta-analysis
  Mean.d = 0 # Holds the mean Cliff's d effect size based on meta-analysis
  d.sig = 0 # Holds the number of times the Cliff's d effect size was significant with a meta-analysis

  Mean.Cohend = 0  # Holds the mean SMD effect size based on meta-analysis

  Cohend.sig = 0 # Holds the number of times the SMD effect size was significant with a meta-analysis

  DataTable = NULL

  FourGroupExp = c(rep(FourGroup, Replications))

  for (i in 1:Replications) {
    if (!FourGroupExp[i]) {
      res = NP2GroupMetaAnalysisSimulation(
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
        MAMethod = MAMethod
      )
    }
    if (FourGroupExp[i]) {
      res = NP4GroupMetaAnalysisSimulation(
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
        MAMethod = MAMethod
      )
    }


    if (!returnES)  {
      #Obtain the values from the simulated family analysis and add to the results of previous simulations
      Averagektau = Averagektau + res$AveKtau
      Averagektauctvar = Averagektauctvar + res$AveKtauctvar
      if (res$tauSigCVt) {
        AveragetauSigCVt = AveragetauSigCVt + 1
      }

      AverageCliffd = AverageCliffd + res$AveCliffd
      AverageCliffdvar = AverageCliffdvar + res$AveCliffdvar
      if (res$AveCliffdsig) {
        AverageCliffdsig = AverageCliffdsig + 1
      }

      Averagephat = Averagephat + res$Avephat
      Averagephatvar = Averagephatvar + res$Avephatvar
      if (res$Avephatsig) {
        Averagephatsig = Averagephatsig + 1
      }


      MAAverage = MAAverage + res$MAMean
      MAVariance = MAVariance + res$MAvar
      if (res$MASig) {
        MASignificant = MASignificant + 1
      }
      if (res$HetSig) {
        HetSignificant = HetSignificant + 1
      }
      MAPcor.mean = MAPcor.mean + res$P.mean
      if (res$P.rsig) {
        MAPcor.sig = MAPcor.sig + 1
      }
      if (res$P.hetsig) {
        MAPcor.hetsig = MAPcor.hetsig + 1
      }


      Mean.phat = Mean.phat + res$Mean.phat
      if (res$phat.sig) {
        phat.sig = phat.sig + 1
      }

      Mean.d = Mean.d + res$Mean.d
      if (res$d.sig) {
        d.sig = d.sig + 1
      }

      Mean.g.exact = Mean.g.exact + res$Mean.g.exact
      if (res$g.exact.sig) {
        g.exact.sig = g.exact.sig + 1
      }

      Mean.g.approx = Mean.g.approx + res$Mean.g.approx
      if (res$g.approx.sig) {
        g.approx.sig = g.approx.sig + 1
      }


      Mean.Cohend = Mean.Cohend + res$Cohend.mean
      if (res$Cohend.sig) {
        Cohend.sig = Cohend.sig + 1
      }
    }
    else {
      # Store the outcome from each replication



      DataTable = tibble::as_tibble(dplyr::bind_rows(DataTable, dplyr::bind_cols(Family = i, res)))
    }

  }

  if (!returnES) {
    #Calculate averages.
    Averagektau = signif(Averagektau / Replications, 4)
    Averagektauctvar = signif(Averagektauctvar / Replications, 4)
    AveragetauSigCVt = signif(AveragetauSigCVt / Replications, 4)

    AverageCliffd = signif(AverageCliffd / Replications, 4)
    AverageCliffdvar = signif(AverageCliffdvar / Replications, 4)
    AverageCliffdsig = signif(AverageCliffdsig / Replications, 4)

    Averagephat = signif(Averagephat / Replications, 4)
    Averagephatvar = signif(Averagephatvar / Replications, 4)
    Averagephatsig = signif(Averagephatsig / Replications, 4)


    MAAverage = signif(MAAverage / Replications, 4)
    MAVariance = signif(MAVariance / Replications, 4)
    MASignificant = signif(MASignificant / Replications, 4)
    HetSignificant = signif(HetSignificant / Replications, 4)
    PMAAverage = signif(MAPcor.mean / Replications, 4)
    PMASignificant = signif(MAPcor.sig / Replications, 4)
    PMAHetSignificant = signif(MAPcor.hetsig / Replications, 4)

    Mean.phat = signif(Mean.phat / Replications, 4)
    phat.sig = signif(phat.sig / Replications, 4)

    Mean.d = signif(Mean.d / Replications, 4)
    d.sig = signif(d.sig / Replications, 4)

    Mean.g.exact = signif(Mean.g.exact / Replications, 4)
    g.exact.sig = signif(g.exact.sig / Replications, 4)

    Mean.g.approx = signif(Mean.g.approx / Replications, 4)
    g.approx.sig = signif(g.approx.sig / Replications, 4)

    Mean.Cohend = signif(Mean.Cohend / Replications, 4)
    Cohend.sig = signif(Cohend.sig / Replications, 4)

  }

  if (!returnES)
  {
    outcome = tibble::tibble(
      Averagektau,
      Averagektauctvar,
      AveragetauSigCVt,
      AverageCliffd,
      AverageCliffdvar,
      AverageCliffdsig,
      Averagephat,
      Averagephatvar,
      Averagephatsig,
      MAAverage,
      MAVariance,
      MASignificant,
      HetSignificant,
      PMAAverage,
      PMASignificant,
      PMAHetSignificant,
      Mean.phat,
      phat.sig,
      Mean.d,
      d.sig,
      Mean.g.exact,
      g.exact.sig,
      Mean.g.approx,
      g.approx.sig,
      Mean.Cohend,
      Cohend.sig
    )
  }

  else
  {
    outcome = DataTable
  }

  return(outcome)
}

####################################################################################################################################################

# Function used to find the expected means,variance and effect sizes from different distributions

#' @title CalculateTheoreticalEffectSizes
#' @description This function constructs the theoretical effect sizes and distribution statistics four (normal, lognormal, Laplace & gamma) given specific parameter values for the distributions
#' @author Barbara Kitchenham and Lech Madeyski
#' @export CalculateTheoreticalEffectSizes
#' @param mean The theoretical central location parameter for the distribution specified by the type parameter.
#' @param std The theoretical spread parameter for the distribution specified by the type parameter.
#' @param type String identifying the distribution, "n" for normal, "ln" for lognormal, "lap" for Laplace, "g" for Gamm
#' @return dataframe containing the expected standardized effect size, mean, variance,skewness and kurtosis statistics for samples from the specifie distribution
#' @examples
#' CalculateTheoreticalEffectSizes(mean=0, std=1, type="l")
#'# A tibble: 1 x 5
#'#   RawMean RawVariance RawEffectSize RawSkewness RawKurtosis
#'#     <dbl>       <dbl>         <dbl>       <dbl>       <dbl>
#'# 1    1.65        4.67         0.763        6.18        88.5
#' CalculateTheoreticalEffectSizes(mean=0, std=1, type="n")
#'#A tibble: 1 x 5
#'#  RawMean RawVariance RawEffectSize RawSkewness RawKurtosis
#'#    <dbl>       <dbl>         <dbl>       <dbl>       <dbl>
#'#1       0           1             0           0           3

CalculateTheoreticalEffectSizes = function(mean, std, type = "n") {
  if (type == "n") {
    # The expected values of a sample from the normal distribution
    RawMean = mean
    RawVariance = std ^ 2
    RawSkewness = 0
    RawKurtosis = 3
  }
  if (type == "l") {
    # The expected values of a sample from the lognormal distribution
    RawMean = exp(mean + std ^ 2 / 2)
    RawVariance = (exp(std ^ 2) - 1) * exp(2 * mean + std ^ 2)
    RawSkewness = (exp(std ^ 2) + 2) * sqrt(exp(std ^ 2) - 1)
    RawKurtosis = exp(4 * std ^ 2) + 2 * exp(2 * std ^ 2) + 3 * exp(2 *
                                                                      std ^ 2) - 3

  }
  if (type == "g") {
    # The expected values of a sample from the gamma distribution
    shape = std
    rate = mean
    RawMean = shape / rate
    RawVariance = shape / rate ^ 2
    RawSkewness = 2 / sqrt(shape)
    RawKurtosis = 6 / shape + 3
  }
  if (type == "lap") {
    # The expected values of a sample from the Laplace distribution

    location = mean
    scale = std
    RawMean = location
    RawVariance = 2 * scale  ^ 2
    RawSkewness = 0
    RawKurtosis = 6

  }
  RawEffectSize = RawMean / sqrt(RawVariance)
  output = tibble::tibble(RawMean,
                          RawVariance,
                          RawEffectSize,
                          RawSkewness,
                          RawKurtosis)
  return(output)
}

#' @title RandomizedDesignEffectSizes
#' @description This function creates the theoretical effect sizes for random samples from one of four different distributions for specified parameter values for the diftribution specified by the type parameter. It assumes there are two samples, one corresponding to a control group and the other to the treatment group. It returns the theoretical effect sizes for a fully randomized experiment.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export RandomizedDesignEffectSizes
#' @param m1 The theoretical mean for the control group
#' @param std1 The theoretical variance for the control group
#' @param m2 The theoretical mean for the treatment group
#' @param std2 The theoretical variance for the treatment group
#' @param type String identifying the distribution, "n" for normal, "ln" for lognormal, "lap" for Laplace, "g" for Gamma
#' @return dataframe containing the expected values of the unstandardized mean difference effect size, the pooled witjin group variance, the standardized mean difference effect size and the point bi-serial correlation.
#' @examples
#' RandomizedDesignEffectSizes(m1=0, std1=1, m2=1, std2=3, type = "n")
#'# A tibble: 1 x 4
#'#      ES   Var StdES  rPBS
#'#   <dbl> <dbl> <dbl> <dbl>
#'# 1     1     5 0.447 0.218

#' RandomizedDesignEffectSizes(m1=0, std1=1, m2=1, std2=3, type = "l")
#'# A tibble: 1 x 4
#'#     ES        Var  StdES    rPBS
#'# <dbl>      <dbl>  <dbl>   <dbl>
#'#1  243. 242552663. 0.0156 0.00780

#'  RandomizedDesignEffectSizes(m1=0, std1=1, m2=0.266, std2=1, type = "l")
#'# A tibble: 1 x 4
#'#     ES   Var StdES   rPBS
#'#  <dbl> <dbl> <dbl>  <dbl>
#'#1 0.502  6.31 0.200 0.0995

RandomizedDesignEffectSizes = function(m1, std1, m2, std2, type = "n") {
  G1.results = CalculateTheoreticalEffectSizes(m1, std1, type = type)
  G2.results = CalculateTheoreticalEffectSizes(m2, std2, type = type)
  ES = G2.results$RawMean - G1.results$RawMean
  Var = (G2.results$RawVariance + G1.results$RawVariance) / 2
  StdES = ES / sqrt(Var)
  rPBS = StdES / sqrt(StdES ^ 2 + 4)
  output = tibble::tibble(ES, Var, StdES, rPBS)
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
#' @param type String identifying the distribution, "n" for normal, "ln" for lognormal, "lap" for Laplace, "g" for Gamma
#' @return dataframe holing the expected unstandardized mean difference effect size, the pooled within group variance, the standardized effect size and the point bi-serial correlation.
#' @examples
#' RandomizedBlockDesignEffectSizes(m1=0,std1=1,m2=1,std2=1,m3=0,std3=1,m4=1,std4=1,BE=1,type="n")
#' #  A tibble: 1 x 4
#' #      ES   Var StdES  rPBS
#' #   <dbl> <dbl> <dbl> <dbl>
#' # 1     1     1     1 0.447

#' RandomizedBlockDesignEffectSizes(m1=0,std1=1,m2=1,std2=1,m3=0,std3=1,m4=1,std4=1,BE=1,type="l")
#'#  A tibble: 1 x 4
#'#      ES   Var StdES  rPBS
#'#   <dbl> <dbl> <dbl> <dbl>
#'# 1  5.27  82.2 0.581 0.279

#' RandomizedBlockDesignEffectSizes(m1=0,std1=1,m2=0.266,std2=1,m3=0,std3=1,m4=0.266,std4=1,BE=0,
#'  type = "l")
#'#  A tibble: 1 x 4
#'#      ES   Var StdES   rPBS
#'#   <dbl> <dbl> <dbl>  <dbl>
#'# 1 0.502  6.31 0.200 0.0995

RandomizedBlockDesignEffectSizes = function(m1,
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
  G1.results = CalculateTheoreticalEffectSizes(m1, std1, type = type)
  G2.results = CalculateTheoreticalEffectSizes(m2, std2, type = type)
  if (type != "g") {
    # Except for the gamma the block effect is applied to the central location parameter
    G3.results = CalculateTheoreticalEffectSizes(m3 + BE, std3, type = type)
    G4.results = CalculateTheoreticalEffectSizes(m4 + BE, std4, type = type)
  }

  if (type == "g") {
    # This function applies the blocking effect to the shape function for the gamma distibution
    G3.results = CalculateTheoreticalEffectSizes(m3 , std3 + BE, type = type)
    G4.results = CalculateTheoreticalEffectSizes(m4 , std4 + BE, type = type)
  }

  # Calculate the expected unstandardized effect size allowing for the experimental design
  ES = (G2.results$RawMean - G1.results$RawMean + G4.results$RawMean - G3.results$RawMean) /
    2
  # Calculate the within groups pooled variance
  Var = (
    G2.results$RawVariance + G1.results$RawVariance + G3.results$RawVariance +
      G4.results$RawVariance
  ) / 4
  # Calculate the standarized mean difference effect size
  StdES = ES / sqrt(Var)
  # Calculate the point bi-serial effect size
  rPBS = StdES / sqrt(StdES ^ 2 + 4)
  output = tibble::tibble(ES, Var, StdES, rPBS)
  return(output)
}
