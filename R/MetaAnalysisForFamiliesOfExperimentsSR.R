##    reproducer package for R
##    Copyright (C) 2018 Lech Madeyski and Barbara Kitchenham
## R functions to complement paper by Barbara Kitchenham, Lech Madeyski, Pearl Brereton, "Meta-analysis for Families of Experiments: A Systematic Review and Reproducibility Assessment".
## This file includes functions that, e.g., allow to construct and format five of the output tables used in the systematic review paper. It extracts the reported values for effect sizes, meta-analysis and descriptive statistics in the primary studies from txt files. It uses the descriptive statistics to re-calculate effect sizes and then performs a meta-analyses using the constructed effect sizes and compares the calculated values with the reported values.


#' @title calculateSmallSampleSizeAdjustment
#' @description Function calculates the small sample size adjustment for standardized mean effect sizes
#' @author Barbara Kitchenham and Lech Madeyski
#' @export calculateSmallSampleSizeAdjustment
#' @param df A vector of degrees of freedom
#' @param exact Default value=TRUE, If exact==TRUE the function returns the exact value of the adjustment(s) which is suitable for small values of df, if exact==FALSE the function returns the approximate version of the adjustment(s). See Hedges and Olkin 'Statistical methods for Meta-Analysis' Academic Press 1985.
#' @return small sample size adjustment value
#' @examples
#' df <- 2
#' c <- calculateSmallSampleSizeAdjustment(df)
#'
#' df=c(5,10,17)
#' adjexact=calculateSmallSampleSizeAdjustment(df)
#' # adjexact=0.8407487 0.9227456 0.9551115
#' # Hedges and Olkin values 0.8408, 0.9228,0.9551
#' adjapprox=calculateSmallSampleSizeAdjustment(df,FALSE)
#' # adjapprox=0.8421053 0.9230769 0.9552239
calculateSmallSampleSizeAdjustment = function(df, exact = TRUE) {
  exactvec = c(rep(exact, length(df)))
  c = ifelse(exactvec, sqrt(2 / df) * gamma(df / 2) / gamma((df - 1) / 2) , (1 - 3 / (4 * df - 1)))
  return(c)
}


#' @title constructEffectSizes
#' @description The function constructs various different d-style effect sizes for a set of different experiments given basic statistics from each experiment ( the mean value of the control group Mc, the mean value of the treatment group Mt, the standard deviation of the control group SDc, standard deviation of the the treatment group SDt, the number of observations (particpants) in the control group Nc, and the number of observations (participants) in the treatment group Nt). The input variables can be vectors or individual numbers but all input vectors must be of the same length. The function returns Glass's Delta, Cohen's D, point bi-serial r (based on Hedges'g unadjusted), Hedges'g and Hegdes' g adjusted for small sample size.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export constructEffectSizes
#' @param Mc is a vector containing the mean value of the control group for each experiment.
#' @param Mt is a vector containing the mean value of the treatment group for each experiment.
#' @param SDc is a vector of the standard deviations of the control group for each experiment.
#' @param SDt is a vector of the standard deviations of the the treatment group for each experiment.
#' @param Nc is a vector containing the the number of observations (particpants) in the control group for each experiment.
#' @param Nt is a vector of the number of observations (participants) in the treatment group for each experiment.
#' @return data frame composed of five effect sizes (Glass delta, Cohen's d, Hedges' g, r, Hedges' g adjusted)
#' @examples
#' constructEffectSizes(10, 15, 0.3, 0.2, 15, 15)
#'
#' Mt = c(0.633, 0.673, 0.423, 0.727, 0.631)
#' Mc = c(0.612, 0.526, 0.356, 0.618, 0.534)
#' SDt = c(0.198, 0.115, 0.172, 0.088, 0.122)
#' SDc = c(0.159, 0.089, 0.111, 0.166, 0.119)
#' Nt = c(12, 12, 14, 10, 8)
#' Nc= c(12, 12, 14, 10, 8)
#' EffectSizes=constructEffectSizes(Mc, Mt, SDc,SDt,Nt,Nc)
#' EffectSizes
#' # GlassDelta    Cohend   Hedgesg          r HedgesgAdjusted
#' # 1  0.1320755 0.1221516 0.1169513 0.05837591       0.1129107
#' # 2  1.6516854 1.4931812 1.4296121 0.58151846       1.3802200
#' # 3  0.6036036 0.4803405 0.4628677 0.22547423       0.4493641
#' # 4  0.6566265 0.8648343 0.8204538 0.37953300       0.7857047
#' # 5  0.8151261 0.8604924 0.8049169 0.37335594       0.7608781
constructEffectSizes = function(Mc, Mt, SDc, SDt, Nc, Nt) {
  GlassDelta = (Mt - Mc) / SDc
  VarPool = (SDc ^ 2 * (Nc - 1) + SDt ^ 2 * (Nt - 1)) / (Nc + Nt - 2)
  SPool = sqrt(VarPool)
  SigmaPool = sqrt((SDc ^ 2 * (Nc - 1) + SDt ^ 2 * (Nt - 1)) / (Nc + Nt))
  Cohend = (Mt - Mc) / (SigmaPool)
  Hedgesg = (Mt - Mc) / SPool
  r = transformHgtoR(Hedgesg, Nc, Nt)
  df = Nc + Nt - 2
  c = calculateSmallSampleSizeAdjustment(df)
  HedgesgAdjusted = Hedgesg * c
  effectsizes = data.frame(GlassDelta, Cohend, Hedgesg, r, HedgesgAdjusted)
  return(effectsizes)
}


#' @title transformRtoZr
#' @description The function transforms a vector of point biserial r values to their normal approximation. It also works for the correlation r.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export transformRtoZr
#' @param r A vector of r-values
#' @return value of normal approximation of point biserial r
#' @examples
#' transformRtoZr(0.4)
#' # [1] 0.4236489
#' Zr=transformRtoZr(c(0.4,0.2))
#' Zr
#' # [1] 0.4236489 0.2027326
transformRtoZr = function(r) {
  zr = 0.5 * log((1 + r) / (1 - r))
  return(zr)
}


#' @title transformZrtoR
#' @description The function transforms a vector of standardized normal variates to their equivalent r-values.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export transformZrtoR
#' @param zr A vector of standard normal variates.
#' @return value of point biserial r
#' @examples
#' transformZrtoR(0.4236489)
#' # [1] 0.4
#' transformZrtoR(c(0.4236489, 0.2027326))
#' # [1] 0.4 0.2
transformZrtoR = function(zr) {
  r = (exp(2 * zr) - 1) / (exp(2 * zr) + 1)
  return(r)
}



#' @title transformHgtoR
#' @description The functions transforms a vector of Hedges g values to their equivalent point bi-serial values.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export transformHgtoR
#' @param g A vector of Hegdes g values.
#' @param Nc A vector of numbers identifying the number of control condition participants in each group
#' @param Nt A vector of numbers identifying the number of treatment condition participants in each group
#' @return value of point biserial r
#' @examples
#' transformHgtoR(0.4, 20, 20)
#' # [1] 0.1961161
transformHgtoR = function(g, Nc, Nt) {
  # Transform Hedges' g to r
  a = (Nc + Nt) ^ 2 / (Nt * Nc)
  r = g / sqrt(g ^ 2 + a)
  return(r)
}


#' @title transformHgtoZr
#' @description The functions transforms a vector of Hedges g values to their normal approximation of point bi-serial values.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export transformHgtoZr
#' @param g value of Hedges' g
#' @param Nc the number of observations (particpants) in the first (control) group
#' @param Nt the number of observations (participants) in the second (treatment) group
#' @return value of normal approximation of point biserial r
#' @examples
#' transformHgtoZr(0.5, 20, 20)
#' # [1] 0.2474665
transformHgtoZr = function(g, Nc, Nt) {
  a = (Nc + Nt) ^ 2 / (Nt * Nc)
  r = g / sqrt(g ^ 2 + a)
  zr = transformRtoZr(r)
  return(zr)
}



#' @title calculateHg
#' @description This function calculates Hedges g and Hedges g adjusted given the basic experimental statistics - the mean values for participants, number of observations (participants), and standard deviation in both the control group and the treatment group. . Hence, the function assumes the data is held as summary statistics including the control group mean, standard deviation and sample size and equivalent values for treatment group
#' @author Barbara Kitchenham and Lech Madeyski
#' @export calculateHg
#' @param Mc is a vector containing the mean value of the control group for each experiment.
#' @param Mt is a vector containing the mean value of the treatment group for each experiment.
#' @param Nc is a vector containing the the number of observations (particpants) in the control group for each experiment.
#' @param Nt is a vector of the number of observations (participants) in the treatment group for each experiment.
#' @param SDc is a vector of the standard deviations of the control group for each experiment.
#' @param SDt is a vector of the standard deviations of the the treatment group for each experiment.
#' @return data frame composed of Hedges' g and Hedges' g adjusted effect sizes
#' @examples
#' calculateHg(10, 15, 20, 20, 2, 2)
#' #    Hg    HgAdjusted
#' # 1  2.5   2.450276
calculateHg = function(Mc, Mt, Nc, Nt, SDc, SDt) {
  # Calculate the pooled within group variance
  varpool = ((Nc - 1) * SDc ^ 2 + (Nt - 1) * SDt ^ 2) / (Nc + Nt - 2)
  spool = sqrt(varpool)
  # Calculate the unadjusted mean difference effect size
  Hg = (Mt - Mc) / spool
  # Calculate the adjusted mean difference effect size
  df = Nc + Nt - 2
  c = calculateSmallSampleSizeAdjustment(df)
  HgAdjusted = Hg * c
  results = data.frame(Hg, HgAdjusted)
  return(results)
}



#' @title transformRtoHg
#' @description This function coverts a vector of point bi-serial r values with associated sample size information back to the mean difference effect size Hedges g.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export transformRtoHg
#' @param r A vector of point bi-serial correlation values.
#' @param Nc A vector of the number of observations in the control condition for the related experiments.
#' @param Nt A vector of the number of observations in the treatment condition for the related experiments.
#' @return value of Hedges' g
#' @examples
#' transformRtoHg(c(0.4,0.2), c(20,20), c(20,20))
#' # [1] 0.8728716 0.4082483
transformRtoHg = function(r, Nc, Nt) {
  # Transform a correlation to the standardized effect size scale.
  a = (Nc + Nt) ^ 2 / (Nc * Nt)
  Hg = (r / sqrt(1 - r ^ 2)) * sqrt(a)
  return(Hg)
}




#' @title transformZrtoHgapprox
#' @description This function provides an approximate transformation from Zr to Hedges g when the number of observations in the treatment and control group are unknown. It is also used to allow the forest plots to display Hedge's g when they are based on r. It is necessary because the transformation function in the forest plot function does not allow any parameters other than effect size used. The function assumes that Nc=Nt and gives the same results as transformZrtoHg when Nc=Nt.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export transformZrtoHgapprox
#' @param Zr A vector of normalised point bi-serial values
#' @return approx. value of Hedges' g
#' @examples
#' transformZrtoHgapprox(c(0.4,0.2))
#' # [1] 0.8215047 0.4026720
transformZrtoHgapprox = function(Zr) {
  r = transformZrtoR(Zr)
  Hg = (r / sqrt(1 - r ^ 2)) * 2
  return(Hg)
}




#' @title transformZrtoHg
#' @description Transforms Zr to Hedge's g.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export transformZrtoHg
#' @param Zr the normal variate
#' @param Nc the number of observations (particpants) in the first (control) group
#' @param Nt the number of observations (participants) in the second (treatment) group
#' @return value of Hedges' g
#' @examples
#' transformZrtoHg(0.5, 20, 20)
#' #[1] 1.042191
transformZrtoHg = function(Zr, Nc, Nt) {
  # Transform from a normal variate back to the equivalent standardized mean effect size
  r = transformZrtoR(Zr)
  Hg = transformRtoHg(r, Nc, Nt)
  return(Hg)
}


#' @title PrepareForMetaAnalysisGtoR
#' @description This function calculates the standardized effect sizes and their confidence intervals, the equivalence point biserial effect size and the Zr and var(Zr) needed for input into the metafor rma function (meta analysis). In this function the point bi-serial effect size is based on the adjusted Hedges g value. The function uses the Hedges g to r transformation to prepare for meta-analysing the data where the mean values, the standard deviations, and the number of observations are available.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export PrepareForMetaAnalysisGtoR
#' @param Mc is a vector containing the mean value of the control group for each experiment.
#' @param Mt is a vector containing the mean value of the treatment group for each experiment.
#' @param SDc is a vector of the standard deviations of the control group for each experiment.
#' @param SDt is a vector of the standard deviations of the the treatment group for each experiment.
#' @param Nc is a vector containing the the number of observations (particpants) in the control group for each experiment.
#' @param Nt is a vector of the number of observations (participants) in the treatment group for each experiment.
#' @return data frame incl. calculated effect sizes (Hedges' g, Hedges' g adjusted), upper and lower confidence bounds on Hedges' g, zr, vi - variance of zr, r and pvalue
#' @examples
#' PrepareForMetaAnalysisGtoR(c(10,10), c(12,14), c(4,4), c(4,4), c(20,20), c(40,40))
#'#HGvalues.Hg HGvalues.HgAdjusted  Hgupper     Hglower        zr         vi         r       pvalue
#'#        0.5           0.4935018 1.082017 -0.06156572 0.2305901 0.01754386 0.2265882 0.0816981743
#'#        1.0           0.9870036 1.634701  0.40620071 0.4499419 0.01754386 0.4218513 0.0006813222
PrepareForMetaAnalysisGtoR = function(Mc, Mt, SDc, SDt, Nc, Nt) {
  # Calculate Adjusted Hedges's g statistic
  HGvalues = calculateHg(Mc, Mt, Nc, Nt, SDc, SDt)

  # Transform Hedge's g to r
  r = transformHgtoR(HGvalues$HgAdjusted, Nc, Nt)

  # Apply Fisher's transformation to normalise r
  zr = transformRtoZr(r)

  # Calculate variance of zr
  vi = 1 / (Nc + Nt - 3)

  # Calculate upper and lower confidence bounds for each study
  zrupper = zr + 1.96 * sqrt(vi)
  zrlower = zr - 1.96 * sqrt(vi)
  rupper = transformZrtoR(zrupper)
  rlower = transformZrtoR(zrlower)
  Hgupper = transformRtoHg(rupper, Nc, Nt)
  Hglower = transformRtoHg(rlower, Nc, Nt)
  Z = zr / sqrt(vi)
  pvalue = ifelse(Z < 0, 2 * stats::pnorm(Z), 2 * (1 - (stats::pnorm(Z))))

  results = data.frame(HGvalues$Hg,
                       HGvalues$HgAdjusted,
                       Hgupper,
                       Hglower,
                       zr,
                       vi,
                       r,
                       pvalue)
  return(results)
}



#' @title ExtractMAStatistics
#' @description This function extracts summary statistics from meta-analysis results obtained from the rma function of the metafor R package. If required the function transform back to standardized mean difference (effect size type "d" i.e. Hg) or point biserial correlations (effect size type "r").
#' Warning: the `ExtractMAStatistics` function works with `metafor` version 2.0-0, but changes to metafor's method of providing access to its individual results may introduce errors into the function.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export ExtractMAStatistics
#' @param maresults is the output from the rma function.
#' @param Nc is the number of participants in the control condition group.
#' @param Nt is the number of participants in the treatment condition group.
#' @param Transform is a boolean value indicating whether the outcome values need to be transformed back to standardized mean difference ("d" i.e. Hg) or point biserial correlations ("r"). It is defaulted to TRUE. If this parameter is set to FALSE, no transformation will be applied.
#' @param type this indicates the type of transformation required - it defaults to "d" which requests transformation from Zr to Hg, using "r" requests transformation from Zr to r.
#' @param sig indicates the number of significant digits requested in the output, the default is 4; it rounds the values of mean, pvalue, upper and lower bound to the specified number of significant digits.
#' @return data frame incl. summary statistics from meta-analysis results: overall mean value for the effect sizes, the p-value of the mean, the upper and lower confidence interval bounds (UB and LB), QE which is the heterogeneity test statistic and QEp which the the p-value of the heterogeneity statistic
#' @examples
#' ExpData=reproducer::KitchenhamMadeyskiBrereton.ExpData
#' #Extract the experiment basic statics
#' S1data=subset(ExpData,ExpData=="S1")
#' #Use the descriptive data to construct effect size
#' S1EffectSizes = reproducer::PrepareForMetaAnalysisGtoR(
#' S1data$Mc,S1data$Mt,S1data$SDc,S1data$SDt,S1data$Nc,S1data$Nt)
#' # Do a random effect meta-analysis of the transformed r_pbs effect size
#' S1MA = metafor::rma(S1EffectSizes$zr, S1EffectSizes$vi)
#' # Extract summary statistics from meta-analysis results and transform back to Hg scale
#' S1MAStats=reproducer::ExtractMAStatistics(S1MA, sum(S1data$Nc),sum(S1data$Nt), TRUE, "d", 4)
#' #    mean   pvalue    UB     LB QE  QEp
#' #1 0.6658 0.002069 1.122 0.2384  4 0.41
ExtractMAStatistics = function(maresults,
                               Nc,
                               Nt,
                               Transform = TRUE,
                               type = "d",
                               sig = 4) {
  pvalue = as.numeric(maresults$pval) # Bug fixed by LM 20189026 - patch: pvalue=as.numeric(maresults[4]) into pvalue=as.numeric(maresults$pval)
  QE = as.numeric(maresults$QE)
  QEp = as.numeric(maresults$QEp)
  mean = as.numeric(maresults$beta)   # Bug fixed by LM 20189026 - patch: maresults[1] into maresults$b
  UB = as.numeric(maresults$ci.ub)  # Bug fixed by LM 20189026 - patch: maresults[6] into maresults$ci.ub
  LB = as.numeric(maresults$ci.lb) # Bug fixed by LM 20189026 - patch: maresults[5] into maresults$ci.lb

  if (Transform & type == "d") {
    mean = transformZrtoHg(mean, Nc, Nt)
    UB = transformZrtoHg(UB, Nc, Nt)
    LB = transformZrtoHg(LB, Nc, Nt)
  }

  if (Transform & type == "r") {
    mean = transformZrtoR(mean)
    UB = transformZrtoR(UB) #bug 2 reported by LM 20180925 (Function transformZrtoHg has 3 parameters with no default values, while in line 154 it is called with only 1 value, which causes an error: argument "Nc" is missing, with no default). error causes subsequent errors so I suggest to start with this one and some other will be fixed as a consequence. Proposed patch: Line 154 UB=transformZrtoHg(UB) change into UB=transformZrtoR(UB) Yes
    LB = transformZrtoR(LB)
  }

  mean = signif(mean, sig)
  pvalue = signif(pvalue, sig)
  UB = signif(UB, sig)
  LB = signif(LB, sig)
  QE = signif(QE, 2)
  QEp = signif(QEp, 2)
  metaanalysisresults = data.frame(mean, pvalue, UB, LB, QE, QEp)
  return(metaanalysisresults)
}


#' @title aggregateIndividualDocumentStatistics
#' @description This function assumes an ABBA crossover experiment has reported means and variances for each technique in each time period. We calculate the weighted mean and pooled within group variance for the observations arising from the two different sets of materials for a specific technique.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export aggregateIndividualDocumentStatistics
#' @param D1.M is a vector of mean values from a set of experiments in a family reporting observations from participants using a specific document in the first time period with either the control or the treatment technique.
#' @param D1.SD is a vector of results from the set of experiment in a family reporting the standard deviations of observations from participants using the same document in the first time period with the same technique.
#' @param D1.N is a vector of the numbers of participants in each experiment in a family, using the same document for participants using either the same technique.
#' @param D2.M is a vector of mean values of observations from participants using the alternative document in the second time period, but using the same technique.
#' @param D2.SD is a vector of the standard deviations of observations from participants using the alternative document in the second time period with the same technique.
#' @param D2.N is a vector of the numbers of participants using the same document in the second time period for participants using the same technique.
#' @return data frame incl. the overall weighted mean and pooled standard deviation
#' @examples
#' aggregateIndividualDocumentStatistics(10, 2, 20, 15, 2, 20)
#' #     M SD
#' #1 12.5  2
aggregateIndividualDocumentStatistics = function(D1.M, D1.SD, D1.N, D2.M, D2.SD, D2.N) {
  #Calculate total observations
  N = D1.N + D2.N
  #Calculate overall mean
  M = (D1.M * D1.N + D2.M * D2.N) / N
  #Calculate pooled variance
  Var = (D1.SD ^ 2 * (D1.N - 1) + D2.SD ^ 2 * (D2.N - 1)) / (N - 2)
  SD = sqrt(Var)
  CombinedData = data.frame(M, SD)
  return(CombinedData)
}


#' @title reproduceTablesOfPaperMetaAnalysisForFamiliesOfExperiments
#' @description This function reproduces five of the output tables used in the systematic review paper "Meta-analysis for Families of Experiments: A Systematic Review and Reproducibility Assessment". It extracts the reported values for effect sizes, meta-analysis and descriptive statistics in the primary studies. It uses the descriptive statistics to re-calculate effect sizes and then performs a meta-analyses using the constructed effect sizes and compares the calculated values with the reported values.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export reproduceTablesOfPaperMetaAnalysisForFamiliesOfExperiments
#' @return list incl. the data presented in five of the tables presented in the paper.
#' @examples
#' rrData = reproduceTablesOfPaperMetaAnalysisForFamiliesOfExperiments()
#' # Reproduce Table "Overall Mean Values of Effect Sizes Reported and Calculated":
#' xtable::xtable(rrData$MAStats)
#' # Reproduce Table "Calculated and Reported Effect Sizes":
#' xtable::xtable(rrData$ESdata)
#' # Report values for 3 papers that reported per document
#' rrData$MAStatsTP1=data.frame(rrData$MAStatsTP1,row.names=NULL)
#' rrData$ESTP1res=data.frame(rrData$ESTP1res,row.names=NULL)
#' xtable::xtable(rrData$MAStatsTP1)
#' xtable::xtable(rrData$ESTP1res)
#' # Report extra results for Study 8
#' # Reproduce Table "Calculating r_PB Effect Size from Probabilities"
#' xtable::xtable(rrData$GH2015extra)
reproduceTablesOfPaperMetaAnalysisForFamiliesOfExperiments = function() {
  #ExperimentEffectSizes holds the data for all studies that reported expeirment level effect sizes. Note S3 only reported experiment effect sizes for one of its experiments. ABBAEffectSizes holds the effect sizes for the three papers that reported results on a time period/document bais.
  #ExperimentEffectSizes=read.table(paste0(dataSetPath, "ReportedEffectSizes.txt"),header=TRUE)
  ExperimentEffectSizes = reproducer::KitchenhamMadeyskiBrereton.ReportedEffectSizes

  ExperimentEffectSizes = data.frame(ExperimentEffectSizes)


  #ABBAEffectSizes=read.table(paste0(dataSetPath, "ABBAReportedEffectSizes.txt"),header=TRUE)
  ABBAEffectSizes = reproducer::KitchenhamMadeyskiBrereton.ABBAReportedEffectSizes

  ABBAEffectSizes = data.frame(ABBAEffectSizes)

  # Read the input meta-analysis results reported for each study.

  #MAReportedResults=read.table(paste0(dataSetPath, "MAResults.txt"),header=TRUE)
  MAReportedResults = reproducer::KitchenhamMadeyskiBrereton.MetaAnalysisReportedResults

  MARepResShort = data.frame(MAReportedResults)

  MARepResShort = reshape::rename(MARepResShort, c("Qep" = "QEp"))

  # Read the document/time period effect sizes reported by studies 3, 7, 11

  #ABBAMAResults=read.table(paste0(dataSetPath, "ABBAMAResults.txt"),header=TRUE)
  ABBAMAResults = reproducer::KitchenhamMadeyskiBrereton.ABBAMetaAnalysisReportedResults

  ABBAMARepResShort = ABBAMAResults
  ABBAMARepResShort = reshape::rename(ABBAMARepResShort, c("Qep" = "QEp"))
  ABBAMARepResShort = data.frame(ABBAMARepResShort)

  #Read files that contain the basic statistics. ExpData holds the data for all the experiments. Notes values in ExpData for studies 3, 7, 11 were gernerated from the data in DocData. DocData holds the basic statistics for the three papers that reported effect sizes for time periods not for experiments.

  #ExpData=read.table(paste0(dataSetPath, "ExpData.txt"),sep=",",header=TRUE)
  ExpData = reproducer::KitchenhamMadeyskiBrereton.ExpData

  #DocData=read.table(paste0(dataSetPath, "DocData.txt"),sep=",",header=TRUE)
  DocData = reproducer::KitchenhamMadeyskiBrereton.DocData

  #***********************************************************
  # This code calculates effect sizes and meta-analysis from the baic statistics and builds up tables to compare the calculated results with the reported results.

  #***********************************************************
  # Study1:S1 Abrahao 2013 Study reported 5 experiments all using a 4-group crossover design

  #Extract the experiment basic statics
  S1data = subset(ExpData, ExpData == "S1")

  #Use the descriptive data to construct effect size

  S1EffectSizes = PrepareForMetaAnalysisGtoR(S1data$Mc,
                                             S1data$Mt,
                                             S1data$SDc,
                                             S1data$SDt,
                                             S1data$Nc,
                                             S1data$Nt)
  # Do a random effect meta-analysis of the transformed r_pbs effect size
  S1MA = metafor::rma(S1EffectSizes$zr, S1EffectSizes$vi)

  # Extract summary statistics from meta-analysis results and transform back to Hg scale
  S1MAStats = ExtractMAStatistics(S1MA, sum(S1data$Nc), sum(S1data$Nt))

  # Restrict the calculated effect sizes to 4 decimal places
  S1.ES = signif(S1EffectSizes$HGvalues.HgAdjusted, 4)

  #Append the calculated meta-analysis statistics and the reported meta-analysis statistics to the file MAStats that will comprise the basic meta-analysis table in the SR paper.

  Study = "S1"
  Source = "Calc"
  # BAK Changed Type to "gIG"
  Type = "gIG"
  MAStats = data.frame(Study, Type, Source, S1MAStats, RR = NA)

  # Extract the entry for S1 from the data.frame holding the reported Meta-analysis results
  S1.MARep = subset(MARepResShort, MARepResShort$Study == "S1")

  S1.MARep$Type="dRM"
  # Asses whether the reported MA statistics have been reproduced
  MAReproduced = ifelse (abs(S1MAStats$mean - S1.MARep$mean) > 0.05, "No", "Yes")

  # Add the reported meta-analysis results and the calculated meta-analysis resutls to the MAStats file
  MAStats = rbind(MAStats, cbind(S1.MARep, RR = MAReproduced))

  # Append the entries for Study 1, to the data.frame used to hold the comparison of the calculated and reported effect sizes.
  # Fist construct a single row data.frame of the calculated effect sizes.
  # BAK Changed Type to "gIG"
  ESdata1 = data.frame(
    Study = Study,
    Type = "gIG",
    Source = "Calc",
    Design = "4GroupCO",
    Exp1 = S1.ES[1],
    Exp2 = S1.ES[2],
    Exp3 = S1.ES[3],
    Exp4 = S1.ES[4],
    Exp5 = S1.ES[5],
    RR = NA
  )

  # Next obtain the reported effect sizes from the data.frame
  S1.ESrep = subset(ExperimentEffectSizes, ExperimentEffectSizes$Study ==
                      "S1")

# BAK Changed Type to "dRM"

S1.ESrep$Type="dRM"

  # Next construct the vector of the reported effect sizes and append to the data.frame holding the calculated effect size
  ESdata1 = data.frame(rbind(ESdata1, cbind(S1.ESrep, RR = NA)))

  # Next check whether the effect sizes have been reproduced. Look at the difference between the rported and calculated effect sizes. An absolute difference less than or equal to 0.05 is our definition of reproduced.
  Diff1 = signif(abs(as.numeric(ESdata1[1, ]$Exp1) - as.numeric(ESdata1[2, ]$Exp1)), 4)
  Diff2 = signif(abs(as.numeric(ESdata1[1, ]$Exp2) - as.numeric(ESdata1[2, ]$Exp2)), 4)
  Diff3 = signif(abs(as.numeric(ESdata1[1, ]$Exp3) - as.numeric(ESdata1[2, ]$Exp3)), 4)
  Diff4 = signif(abs(as.numeric(ESdata1[1, ]$Exp4) - as.numeric(ESdata1[2, ]$Exp4)), 4)
  Diff5 = signif(abs(as.numeric(ESdata1[1, ]$Exp5) - as.numeric(ESdata1[2, ]$Exp5)), 4)
  # Count the number of times the effect size was reproduced
  ESReproduce = 0
  ESReproduce = ESReproduce + ifelse(Diff1 > 0.05, 0, 1)
  ESReproduce = ESReproduce + ifelse(Diff2 > 0.05, 0, 1)
  ESReproduce = ESReproduce + ifelse(Diff3 > 0.05, 0, 1)
  ESReproduce = ESReproduce + ifelse(Diff4 > 0.05, 0, 1)
  ESReproduce = ESReproduce + ifelse(Diff5 > 0.05, 0, 1)
  # Add the resulting count to the RR field in the reported  entry of the data.frame
  ESdata1[2, ]$RR = ESReproduce

  #Copy the two effect size rows to the data.frame that will hold the table of effect size results in the SR paper.
  ESdata = ESdata1

  #***********************************************************
  # Study: S2 Scanniello 2014. Study reported 4 experiments all using a 4-group crossover design

  #Extract the experiment basic statics
  S2data = subset(ExpData, ExpData == "S2")

  # Convert the reported descritive data into effect sizes
  S2EffectSizes = PrepareForMetaAnalysisGtoR(S2data$Mc,
                                             S2data$Mt,
                                             S2data$SDc,
                                             S2data$SDt,
                                             S2data$Nc,
                                             S2data$Nt)

  # Extract the Hg values and restrict to 4 significant figures
  S2.ES = signif(S2EffectSizes$HGvalues.HgAdjusted, 4)

  # Do a random effect meta-analysis of the transformed r effect sizes
  S2MA = metafor::rma(S2EffectSizes$zr, S2EffectSizes$vi)

  #Transform the summary statistics back to Hg values

  S2MAStats = ExtractMAStatistics(S2MA, sum(S2data$Nc), sum(S2data$Nt))

  #Construct the entries for the meta-analysis results file for the SR paper

  Study = "S2"
  Type = "gIG" # BAK changed from "g"
  Source = "Calc"

  #Extract the reported data for Study S2
  S2RepRes = subset(MARepResShort, MARepResShort$Study == "S2")
  S2RepRes$Type="dIG"
  # Add the calculated data for S2 to the output table
  MAStats = rbind(MAStats, cbind(Study, Type, Source, S2MAStats, RR = NA))

  # Check whether the difference between the calculated and reported overall mean is less than 0.05 to identify reproducible results
  MAReproduced = ifelse (abs(S2MAStats$mean - S2RepRes$mean) > 0.05, "No", "Yes")

  #Add the reported meta-analysis results for S2 and the reproducibily assessment to the output file

  MAStats = rbind(MAStats, cbind(S2RepRes, RR = MAReproduced))

  # Extract the reported effect size data for study 2
  S2ESrep = subset(ExperimentEffectSizes, ExperimentEffectSizes$Study =="S2")
  S2ESrep$Type="dRM"
  # Add the calculated S2 effect sizes to a temporary data.frame
  ESdata1 = data.frame(
    Study = Study,
    Type = "gIG",
    Source = "Calc",
    Design = "4GroupCO",
    Exp1 = S2.ES[1],
    Exp2 = S2.ES[2],
    Exp3 = S2.ES[3],
    Exp4 = S2.ES[4],
    Exp5 = S2.ES[5],
    RR = NA
  )
  # Add the reported S2 effects sizes to the temporary data.frame
  ESdata1 = data.frame(rbind(ESdata1, cbind(S2ESrep, RR = NA)))

  # Calculated the absolute difference betwee the calculated and reported effect sizes for study 2
  Diff1 = signif(abs(as.numeric(ESdata1[1, ]$Exp1) - as.numeric(ESdata1[2, ]$Exp1)), 4)
  Diff2 = signif(abs(as.numeric(ESdata1[1, ]$Exp2) - as.numeric(ESdata1[2, ]$Exp2)), 4)
  Diff3 = signif(abs(as.numeric(ESdata1[1, ]$Exp3) - as.numeric(ESdata1[2, ]$Exp3)), 4)
  Diff4 = signif(abs(as.numeric(ESdata1[1, ]$Exp4) - as.numeric(ESdata1[2, ]$Exp4)), 4)

  # Count the number of times the difference was less than 0.05 and copy the final count into the related column of the output file
  ESReproduce = 0
  ESReproduce = ESReproduce + ifelse(Diff1 > 0.05, 0, 1)
  ESReproduce = ESReproduce + ifelse(Diff2 > 0.05, 0, 1)
  ESReproduce = ESReproduce + ifelse(Diff3 > 0.05, 0, 1)
  ESReproduce = ESReproduce + ifelse(Diff4 > 0.05, 0, 1)

  ESdata1[2, ]$RR = ESReproduce

  # Add the temporary data.frame to the effect size output file.
  ESdata = rbind(ESdata, ESdata1)

  #***********************************************************

  # Study 3: S3 Cruz-Lemus 2009. The authors reported the results of 5 experiments. The first four experiments used an ABBA crossover design. The final experiment used an independent groups design. For the first four experiments, the authors reported the basic statistics for each document, so the descriptive data needed to be re-calculated for each technique group. They undertook their meta-analysis on the basis of the results for each domain (i.e. document). Comoarison of effect sizes will be done on the for the first document data only.

  #Extract the experiment basic statics, for the overall experiment and the first document

  S3data = subset(ExpData, ExpData == "S3")

  Doc1Data = subset(DocData, DocData$Study == "S3" &
                      DocData$Doc == "Doc1")


  if (length(S3data$Study) < 5) {
    # The data for S3 has been changed and the values for experiments 1 to 4 need to be recalculated from the data for each document. This requires the integrating the document 1 and document 2 data into an effect size for each experiment.
    Doc2Data = subset(DocData, DocData$Study == "S3" &
                        DocData$Doc == "Doc2")

    c.Data = aggregateIndividualDocumentStatistics(Doc1Data$Mc,
                                                   Doc1Data$SDc,
                                                   Doc1Data$Nc,
                                                   Doc2Data$Mc,
                                                   Doc2Data$SDc,
                                                   Doc2Data$Nc)
    t.Data = aggregateIndividualDocumentStatistics(Doc1Data$Mt,
                                                   Doc1Data$SDt,
                                                   Doc1Data$Nt,
                                                   Doc2Data$Mt,
                                                   Doc2Data$SDt,
                                                   Doc2Data$Nt)
    r = rep("NA", 4)
    Study = rep("S3", 4)
    Exp4 = c("Exp1", "Exp2", "Exp3", "Exp4")

    S3extradata = data.frame(
      Study = Study,
      Exp = Exp4,
      Mc = signif(c.Data$M, 4),
      SDc = signif(c.Data$S, 6),
      Nc = Doc1Data$Nc,
      Mt = signif(t.Data$M, 4),
      SDt = signif(t.Data$S, 6),
      Nt = Doc1Data$Nt,
      r = r
    )
    # Store the recalculated data values back into the ExpData file, so they dont need to be recalculated again unless the document values chnage again.
    ExpData = rbind(ExpData, S3extradata)
    S3data = rbind(S3extradata, S3data)
  }
  #Calculate the Hg effect sizes that will be used for comparison with the reported effect sizes and the transformed r statistics that will be used for the meta-analysis
  S3EffectSizes = PrepareForMetaAnalysisGtoR(S3data$Mc,
                                             S3data$Mt,
                                             S3data$SDc,
                                             S3data$SDt,
                                             S3data$Nc,
                                             S3data$Nt)

  # Do a random effects analysis after transforming the small sample size adjusted mean effect size to the equivalent r_pb and then to zr
  S3MA = metafor::rma(S3EffectSizes$zr, S3EffectSizes$vi)

  # Transform the meta-analysis summary metrics to Hg Values
  S3MAStats = ExtractMAStatistics(S3MA, sum(S3data$Nc), sum(S3data$Nt))

  S3.ES = signif(S3EffectSizes$HGvalues.HgAdjusted, 4)

  #Construct the entry for the calculated meta-analysis results that will be stored in the data.frame holding the meta-analysis results for all the studies in the SR
  Study = "S3"
  Type = "gIG"
  Source = "Calc"

  # Store the calculated meta-analysis statistics in the MAStats data.frame
  MAStats = rbind(MAStats, cbind(Study, Type, Source, S3MAStats, RR = NA))

  # Extract the reported meta-analysis results
  S3RepRes = subset(MARepResShort, MARepResShort$Study == "S3")
  S3RepRes$Type="dIG"
  # Check whether the results have been reproduced
  MAReproduced = ifelse (abs(S3MAStats$mean - S3RepRes$mean) > 0.05, "No", "Yes")
  # Store the reported results in the MAStats data.frame
  MAStats = rbind(MAStats, cbind(S3RepRes, RR = MAReproduced))

  # Construct the table entries for the reported and calculated effect sizes

  # Construct a temprary data.frame holding the calculated experiment effect sizes
  ESdata1 = data.frame(
    Study = Study,
    Type = "gIG",
    Source = "Calc",
    Design = "Mixed",
    Exp1 = S3.ES[1],
    Exp2 = S3.ES[2],
    Exp3 = S3.ES[3],
    Exp4 = S3.ES[4],
    Exp5 = S3.ES[5],
    RR = NA
  )

  # Extract the reported experiment effect sizes. Only experiment 5 data exists for S3
  S3ESRep = subset(ExperimentEffectSizes, ExperimentEffectSizes$Study ==
                     "S3")

  # Add The reported results to the temporary data.frame
  S3ESRep$Type="dIG"
  ESdata1 = data.frame(rbind(ESdata1, cbind(S3ESRep, RR = NA)))

  # Check whether the results for experiment 5 were reproduced
  Diff5 = signif(abs(as.numeric(ESdata1[1, ]$Exp5) - as.numeric(ESdata1[2, ]$Exp5)), 4)

  ESReproduce = ifelse((Diff5 > 0.05), 0, 1)

  ESdata1[2, ]$RR = ESReproduce

  # Append the temporary dat.frame to the data.frame holding the effect size comparisons
  ESdata = rbind(ESdata, ESdata1)

  # Single document analysis. Do a meta-analysis of document 1.

  # Remove the document identifier column from the document 1 data

  myvars <- names(Doc1Data) %in% c("Doc")
  Doc1.Data <- Doc1Data[!myvars]

  myvars = names(S3data) %in% c("r")
  S3data.red = S3data[!myvars]

  # Doc1.Data holds the basic data for the 4 experiments that used an ABBA crossover. Add the entry for the 5th experiment to the document 1 data
  Doc1.Data = rbind(Doc1.Data, S3data.red[5, ])

  # Prepare to analyse the document 1 /time period 1 data as a set of independent groups exeperiments.
  S3Doc1EffectSizes = PrepareForMetaAnalysisGtoR(
    Doc1.Data$Mc,
    Doc1.Data$Mt,
    Doc1.Data$SDc,
    Doc1.Data$SDt,
    Doc1.Data$Nc,
    Doc1.Data$Nt
  )

  # Do a random effects analysis using the r_pb transformation and its normalization transformation

  S3Doc1MA = metafor::rma(S3Doc1EffectSizes$zr, S3Doc1EffectSizes$vi)

  # Extract the document 1 calculated meta-analysis results

  S3Doc1MAStats = ExtractMAStatistics(S3Doc1MA, sum(Doc1.Data$Nc), sum(Doc1.Data$Nt))

  # Create the first entry in the MAStatsTP1 that will hold the meta-analysis results for the studies that reported document/time period based effect sizes.
  Study = "S3"
  MAStatsTP1 = cbind(Study, Type = "gIG", Source = "Calc", S3Doc1MAStats)

  # Extract the reported meta-analysis resuls for the first time period and append the data to the MAStatsTP1 file.

  S3Doc1RepRes = subset(ABBAMARepResShort, ABBAMARepResShort$Study == "S3")

  MAStatsTP1 = rbind(MAStatsTP1, S3Doc1RepRes)

  # Create a records for the calculated Doc1 effect sizes
  S3Doc1.ES = signif(S3Doc1EffectSizes$HGvalues.HgAdjusted, 4)

  ESTP1res = data.frame(
    Study = Study,
    Type = "gIG",
    Design = "Mixed",
    Source = "Calc",
    Exp1 = S3Doc1.ES[1],
    Exp2 = S3Doc1.ES[2],
    Exp3 = S3Doc1.ES[3],
    Exp4 = S3Doc1.ES[4],
    Exp5 = S3Doc1.ES[5]
  )

  # Extract the reported document/time period 1 effect sizes
  S3Doc1ESRep = subset(ABBAEffectSizes, ABBAEffectSizes$Study == "S3")

  # Add the calculated and reported effect size records to the data.frame ESTP1res that will hold the effect size data for S3, S7 and S11
  ESTP1res = rbind(ESTP1res, S3Doc1ESRep)

  #*************************************************************************
  # Study 4: S4. Fernandez 2013 reported resuls for 3 experiments al.l of which used a 4-group ABBA crossover design

  # Obtain basic statistics for all experiments.
  S4data = subset(ExpData, ExpData == "S4")

  # Calculate effect sizes for meta-analysis
  S4EffectSizes = PrepareForMetaAnalysisGtoR(S4data$Mc,
                                             S4data$Mt,
                                             S4data$SDc,
                                             S4data$SDt,
                                             S4data$Nc,
                                             S4data$Nt)

  # Perform random effect meta-analysis
  S4MA = metafor::rma(S4EffectSizes$zr, S4EffectSizes$vi)

  # Extract the standard meta-analysis results
  S4MAStats = ExtractMAStatistics(S4MA, sum(S4data$Nc), sum(S4data$Nt))

  S4.ES = signif(S4EffectSizes$HGvalues.HgAdjusted, 4)

  # Copy the meta-analysis results to MAStats
  Study = "S4"
  Type = "gIG"
  Source = "Calc"
  MAStats = rbind(MAStats, cbind(Study, Type, Source, S4MAStats, RR = NA))
  # Obtain the reported meta-analysis resutls
  S4RepResShort = subset(MARepResShort, MARepResShort$Study == "S4")
  S4RepResShort$Type="dIG"

  # Check whether the meta-analysis results were reproduced
  MAReproduced = ifelse (abs(S4MAStats$mean - MARepResShort[4, ]$mean) > 0.05, "No", "Yes")
  # Copy the reported results and reproducibility assessment into MAStats
  MAStats = rbind(MAStats, cbind(S4RepResShort, RR = MAReproduced))

  # Compare the reported and calculated effect sizes and copy the results and the reproducibility assessment into ESdata
  ESdata1 = data.frame(
    Study = Study,
    Type = "gIG",
    Source = "Calc",
    Design = "4GroupCO",
    Exp1 = S4.ES[1],
    Exp2 = S4.ES[2],
    Exp3 = S4.ES[3],
    Exp4 = NA,
    Exp5 = NA,
    RR = NA
  )

  # Get the reported effect size values
  S4ES.rep = subset(ExperimentEffectSizes, ExperimentEffectSizes == "S4")
  S4ES.rep$Type="dIG"
  ESdata1 = data.frame(rbind(ESdata1, cbind(S4ES.rep, RR = NA)))
  Diff1 = signif(abs(as.numeric(ESdata1[1, ]$Exp1) - as.numeric(ESdata1[2, ]$Exp1)), 4)
  Diff2 = signif(abs(as.numeric(ESdata1[1, ]$Exp2) - as.numeric(ESdata1[2, ]$Exp2)), 4)
  Diff3 = signif(abs(as.numeric(ESdata1[1, ]$Exp3) - as.numeric(ESdata1[2, ]$Exp3)), 4)

  ESReproduce = 0
  ESReproduce = ESReproduce + ifelse(Diff1 > 0.05, 0, 1)
  ESReproduce = ESReproduce + ifelse(Diff2 > 0.05, 0, 1)
  ESReproduce = ESReproduce + ifelse(Diff3 > 0.05, 0, 1)

  ESdata1[2, ]$RR = ESReproduce

  ESdata = rbind(ESdata, ESdata1)

  #***********************************************************
  #Fernandez-Saez 2016 (pre-publication date 2014) Reported 5 experiments all using the 4-group crossover design.

  # Get the reported basic statsitics for Study 5 experiments
  S5data = subset(ExpData, ExpData == "S5")

  # Calculate the effect sizes
  S5EffectSizes = PrepareForMetaAnalysisGtoR(S5data$Mc,
                                             S5data$Mt,
                                             S5data$SDc,
                                             S5data$SDt,
                                             S5data$Nc,
                                             S5data$Nt)

  # Perform a randomised effect meta-analysis
  S5MAres = metafor::rma(S5EffectSizes$zr, S5EffectSizes$vi)

  # Extract & Transform the meta-analysis summary statistics back to Hg values

  S5MAStats = ExtractMAStatistics(S5MAres, sum(S5data$Nc), sum(S5data$Nt))

  # Calculate the entries for MAStats and ESdata
  Study = "S5"

  S5.ES = signif(S5EffectSizes$HGvalues.HgAdjusted, 4)

  Type = "gIG"
  Source = "Calc"
  MAStats = rbind(MAStats, cbind(Study, Type, Source, S5MAStats, RR = NA))

  # Get the reported meta-analysis results

  S5MARepRes = subset(MARepResShort, MARepResShort$Study == "S5")
  S5MARepRes$Type="dIG"
  MAReproduced = ifelse (abs(S5MAStats$mean - S5MARepRes$mean) > 0.05, "No", "Yes")

  MAStats = rbind(MAStats, cbind(S5MARepRes, RR = MAReproduced))

  ESdata1 = data.frame(
    Study = Study,
    Type = "gIG",
    Source = "Calc",
    Design = "4GroupCO",
    Exp1 = S5.ES[1],
    Exp2 = S5.ES[2],
    Exp3 = S5.ES[3],
    Exp4 = S5.ES[4],
    Exp5 = S5.ES[5],
    RR = NA
  )

  # Index correct table entry
  S5EffectSizes = subset(ExperimentEffectSizes, ExperimentEffectSizes$Study ==
                           "S5")
  S5EffectSizes$Type="dIG"
  ESdata1 = data.frame(rbind(ESdata1, cbind(S5EffectSizes, RR = NA)))

  Diff1 = signif(abs(as.numeric(ESdata1[1, ]$Exp1) - as.numeric(ESdata1[2, ]$Exp1)), 4)
  Diff2 = signif(abs(as.numeric(ESdata1[1, ]$Exp2) - as.numeric(ESdata1[2, ]$Exp2)), 4)
  Diff3 = signif(abs(as.numeric(ESdata1[1, ]$Exp3) - as.numeric(ESdata1[2, ]$Exp3)), 4)
  Diff4 = signif(abs(as.numeric(ESdata1[1, ]$Exp4) - as.numeric(ESdata1[2, ]$Exp4)), 4)

  ESReproduce = 0
  ESReproduce = ESReproduce + ifelse(Diff1 > 0.05, 0, 1)
  ESReproduce = ESReproduce + ifelse(Diff2 > 0.05, 0, 1)
  ESReproduce = ESReproduce + ifelse(Diff3 > 0.05, 0, 1)
  ESReproduce = ESReproduce + ifelse(Diff4 > 0.05, 0, 1)

  ESdata1[2, ]$RR = ESReproduce

  ESdata = rbind(ESdata, ESdata1)

  #***********************************************************
  # Hadar 2013. Used Cohen's d. 3 experiments all used a 4-group crossover

  # Obtain basic statsitics reported in paper
  S6data = subset(ExpData, ExpData$Study == "S6")

  # Obtain effect sizes
  S6EffectSizes = PrepareForMetaAnalysisGtoR(S6data$Mc,
                                             S6data$Mt,
                                             S6data$SDc,
                                             S6data$SDt,
                                             S6data$Nc,
                                             S6data$Nt)

  # Effect sizes confirm that Hadar used the mean difference divided by the pooled within groups SD. However, we still use the small sample correction.
  S6MAres = metafor::rma(S6EffectSizes$zr, S6EffectSizes$vi)

  # Tranform the summary metrics to Hg Values

  S6MAStats = ExtractMAStatistics(S6MAres, sum(S6data$Nc), sum(S6data$Nt))

  Study = "S6"

  S6.ES = signif(S6EffectSizes$HGvalues.HgAdjusted, 4)

  Type = "gIG"
  Source = "Calc"
  MAStats = rbind(MAStats, cbind(Study, Type, Source, S6MAStats, RR = NA))
  # Index correct table entry
  S6MARep = subset(MARepResShort, MARepResShort$Study == "S6")
  S6MARep$Type="dIG"
  MAReproduced = ifelse (abs(S6MAStats$mean - S6MARep$mean) > 0.05, "No", "Yes")
  MAStats = rbind(MAStats, cbind(S6MARep, RR = MAReproduced))

  ESdata1 = data.frame(
    Study = Study,
    Type = "gIG",
    Source = "Calc",
    Design = "4GroupCO",
    Exp1 = S6.ES[1],
    Exp2 = S6.ES[2],
    Exp3 = S6.ES[3],
    Exp4 = NA,
    Exp5 = NA,
    RR = NA
  )

  S6Rep.ES = subset(ExperimentEffectSizes, ExperimentEffectSizes$Study ==
                      "S6")
  S6Rep.ES$Type="dIG"
  ESdata1 = data.frame(rbind(ESdata1, cbind(S6Rep.ES, RR = NA)))

  Diff1 = signif(abs(as.numeric(ESdata1[1, ]$Exp1) - as.numeric(ESdata1[2, ]$Exp1)), 4)
  Diff2 = signif(abs(as.numeric(ESdata1[1, ]$Exp2) - as.numeric(ESdata1[2, ]$Exp2)), 4)
  Diff3 = signif(abs(as.numeric(ESdata1[1, ]$Exp3) - as.numeric(ESdata1[2, ]$Exp3)), 4)


  ESReproduce = 0
  ESReproduce = ESReproduce + ifelse(Diff1 > 0.05, 0, 1)
  ESReproduce = ESReproduce + ifelse(Diff2 > 0.05, 0, 1)
  ESReproduce = ESReproduce + ifelse(Diff3 > 0.05, 0, 1)
  ESdata1[2, ]$RR = ESReproduce

  ESdata = rbind(ESdata, ESdata1)



  #*************************************************************
  # Study 7: S7. Teruel 2012. The authors reported results of 3 experiments all using an ABBA crossover design. The authors obtained effect sizes for each document believeing they wanted to aggregate effect sizes for each application. Since they used an AB/BA crossover, the results for the first document (Jigsaw) could be used as a between groups study.

  # Obtain calculated experiment level effect sizes
  S7data = subset(ExpData, ExpData$Study == "S7")

  # Obtain reported doc1/first time period effect sizes
  Doc1data = subset(DocData, DocData$Study == "S7" &
                      DocData$Doc == "Doc1")

  if (length(S7data$Study) < 3) {
    # The Doc data has been revised and the calculated basic statistics need to be recalculated
    Doc2data = subset(DocData, DocData$Study == "S7" &
                        DocData$Doc == "Doc2")


    c.Data = aggregateIndividualDocumentStatistics(Doc1data$Mc,
                                                   Doc1data$SDc,
                                                   Doc1data$Nc,
                                                   Doc2data$Mc,
                                                   Doc2data$SDc,
                                                   Doc2data$Nc)

    t.Data = aggregateIndividualDocumentStatistics(Doc1data$Mt,
                                                   Doc1data$SDt,
                                                   Doc1data$Nt,
                                                   Doc2data$Mt,
                                                   Doc2data$SDt,
                                                   Doc2data$Nt)


    Study = rep("S7", 3)

    Exp3 = c("Exp1", "Exp2", "Exp3")

    r3 = rep(NA, 3)

    S7data = data.frame(
      Study = Study,
      Exp = Exp3,
      Mc = signif(c.Data$M, 4),
      SDc = signif(c.Data$SD, 6),
      Nc = Doc1data$Nc,
      Mt = signif(t.Data$M, 4),
      SDt = signif(t.Data$SD, 6),
      Nt = Doc1data$Nt,
      r = r3
    )


    ExpData = rbind(ExpData, S7data)
  }

  S7EffectSizes = PrepareForMetaAnalysisGtoR(S7data$Mc,
                                             S7data$Mt,
                                             S7data$SDc,
                                             S7data$SDt,
                                             S7data$Nc,
                                             S7data$Nt)

  S7MAres = metafor::rma(S7EffectSizes$zr, S7EffectSizes$vi)


  # Transform the summary metrics to Hg Values
  S7MAStats = ExtractMAStatistics(S7MAres, sum(Doc1data$Nc), sum(Doc1data$Nt))

  Study = "S7"

  S7.ES = signif(S7EffectSizes$HGvalues.HgAdjusted, 4)

  Type = "gIG"
  Source = "Calc"
  MAStats = rbind(MAStats, cbind(Study, Type, Source, S7MAStats, RR = NA))
  # Index correct table entry
  S7MARep = subset(MARepResShort, MARepResShort$Study == "S7")
  S7MARep$Type="dIG"
  MAReproduced = ifelse (abs(S7MAStats$mean - S7MARep$mean) > 0.05, "No", "Yes")
  MAStats = rbind(MAStats, cbind(S7MARep, RR = MAReproduced))

  ESdata1 = data.frame(
    Study = Study,
    Type = "gIG",
    Source = "Calc",
    Design = "ABBACO",
    Exp1 = S7.ES[1],
    Exp2 = S7.ES[2],
    Exp3 = S7.ES[3],
    Exp4 = NA,
    Exp5 = NA,
    RR = NA
  )

  # Do we put in an empty row here to indicate no effect sizes were formally reported? Currently code for the null row is commented out
  #S7ES=subset(ExperimentEffectSizes,ExperimentEffectSizes$Study=="S7")
  # ESdata1=rbind(ESdata1,cbind(S7ES,RR=NA))

  ESdata = rbind(ESdata, ESdata1)

  # Analyse the first time period/doc1 results for reproducibility. The document 1 only results are analysed as an independent groups study. For this analysis we use the Hg values directly and the appropriate Hg variance


  S7TP1 = PrepareForMetaAnalysisGtoR(Doc1data$Mc,
                                     Doc1data$Mt,
                                     Doc1data$SDc,
                                     Doc1data$SDt,
                                     Doc1data$Nc,
                                     Doc1data$Nt)

  #
  #Construct normal approximation variance of the effect size
  c.Doc1 = calculateSmallSampleSizeAdjustment(Doc1data$Nc + Doc1data$Nt -
                                                2)

  # Use the adjusted Hg but the standard normal approximation equation, which is (n1+n2)/(n1*n2)+ g^2/(2c^2(n1+n2))
  Doc1.VdUnadjusted = ((Doc1data$Nc + Doc1data$Nt) / (Doc1data$Nc * Doc1data$Nt) +
                         S7TP1$HGvalues.HgAdjusted ^ 2 / (2 * c.Doc1 ^ 2 * (Doc1data$Nc + Doc1data$Nt))
  )

  # Do a meta-analysis of the adjusted Hg values using the ML method
  S7TP1res = metafor::rma(S7TP1$HGvalues.Hg, Doc1.VdUnadjusted, method =
                            "ML")

  #Copy the meta-analysis resutls & effect sizes into the data.frames MAStatsTP1 and  ESTP1res which hold the data for studies that aggregated results on  document basis

  S7TP1MAStats = ExtractMAStatistics(
    S7TP1res,
    sum(Doc1data$Nc),
    sum(Doc1data$Nt),
    Transform = F,
    sig = 5
  )

  MAStatsTP1 = rbind(MAStatsTP1,
                     cbind(Study, Type = "d", Source = "Calc", S7TP1MAStats))

  S7TP1MArep = subset(ABBAMARepResShort, ABBAMARepResShort$Study == "S7")

  MAStatsTP1 = rbind(MAStatsTP1, S7TP1MArep)

  S7TP1Hg = signif(S7TP1$HGvalues.Hg, 5)

  #ESTP1res=rbind(ESTP1res,cbind(Study,Type="d",Source="Calc",Design="ABBACO",Exp1=S7TP1Hg[1],Exp2=S7TP1Hg[2],Exp3=S7TP1Hg[3],Exp4=NA,Exp5=NA))
  ESTP1res = rbind(
    ESTP1res,
    data.frame(
      Study,
      Type = "d",
      Source = "Calc",
      Design = "ABBACO",
      Exp1 = S7TP1Hg[1],
      Exp2 = S7TP1Hg[2],
      Exp3 = S7TP1Hg[3],
      Exp4 = NA,
      Exp5 = NA
    )
  ) #issue: character instead of numeric values (see str(ESTP1res)), solution: cbind -> data.frame /LM

  S7ESRep = subset(ABBAEffectSizes, ABBAEffectSizes$Study == "S7")

  ESTP1res = rbind(ESTP1res, S7ESRep)

  #***********************************************************
  # Study 8: S8.Gonzalez-Huerta, 2015. This paper used the point bi-serial correlation effect size. The authors reported the results of 4 experiments all using a 4-group crossover design. They used p-values to calculate r_pb.


  S8data = subset(ExpData, ExpData$Study == "S8")

  S8EffectSizes = PrepareForMetaAnalysisGtoR(S8data$Mc,
                                             S8data$Mt,
                                             S8data$SDc,
                                             S8data$SDt,
                                             S8data$Nt,
                                             S8data$Nc)

  S8MAres = metafor::rma(S8EffectSizes$zr, S8EffectSizes$vi)

  #Transform the summary statistics back to r values

  S8MAstatistics = ExtractMAStatistics(S8MAres, sum(S8data$Nt), sum(S8data$Nt), type =
                                         "r")

  Study = "S8"

  S8.ES = signif(S8EffectSizes$r, 4)

  Type = "r"
  Source = "Calc"
  MAStats = rbind(MAStats,
                  cbind(
                    Study,
                    Type = "rpb",
                    Source,
                    S8MAstatistics,
                    RR = NA
                  ))
  # Extract the reported meta-analysis data
  S8MArep = subset(MARepResShort, MARepResShort$Study == "S8")
  S8MArep$Type="rpb"
  MAReproduced = ifelse (abs(S8MAstatistics$mean - S8MArep$mean) > 0.05, "No", "Yes")
  MAStats = rbind(MAStats, cbind(S8MArep, RR = MAReproduced))

  ESdata1 = data.frame(
    Study = Study,
    Type = "rpb",
    Source = "Calc",
    Design = "4GroupCO",
    Exp1 = S8.ES[1],
    Exp2 = S8.ES[2],
    Exp3 = S8.ES[3],
    Exp4 = S8.ES[4],
    Exp5 = NA,
    RR = NA
  )

  # Index correct table entry
  S8ESrep = subset(ExperimentEffectSizes,
                   ExperimentEffectSizes$Study == "S8" & Type == "r")
  S8ESrep$Type="rpb"
  ESdata1 = data.frame(rbind(ESdata1, cbind(S8ESrep, RR = NA)))

  Diff1 = signif(abs(as.numeric(ESdata1[1, ]$Exp1) - as.numeric(ESdata1[2, ]$Exp1)), 4)
  Diff2 = signif(abs(as.numeric(ESdata1[1, ]$Exp2) - as.numeric(ESdata1[2, ]$Exp2)), 4)
  Diff3 = signif(abs(as.numeric(ESdata1[1, ]$Exp3) - as.numeric(ESdata1[2, ]$Exp3)), 4)
  Diff4 = signif(abs(as.numeric(ESdata1[1, ]$Exp4) - as.numeric(ESdata1[2, ]$Exp4)), 4)

  ESReproduce = 0
  ESReproduce = ESReproduce + ifelse(Diff1 > 0.05, 0, 1)
  ESReproduce = ESReproduce + ifelse(Diff2 > 0.05, 0, 1)
  ESReproduce = ESReproduce + ifelse(Diff3 > 0.05, 0, 1)
  ESReproduce = ESReproduce + ifelse(Diff4 > 0.05, 0, 1)

  ESdata1[2, ]$RR = ESReproduce

  ESdata = rbind(ESdata, ESdata1)

  #Reported r_pbs values for Effectiveness
  GHrepr = c(S8ESrep$Exp1, S8ESrep$Exp2, S8ESrep$Exp3, S8ESrep$Exp4)

  N.Participants = S8data$Nc + S8data$Nt
  N.Obs = 2 * N.Participants

  # Reported p-values.Checking the r-values calculated by GH using the method of converting p-values to Z values and then coverting to r using r=Z/sqrt(N)

  S8ESprep = subset(ExperimentEffectSizes,
                    ExperimentEffectSizes$Study == "S8" & Type == "p")
  S8ESp = c(S8ESprep$Exp1, S8ESprep$Exp2, S8ESprep$Exp3, S8ESprep$Exp4)
  S8p.Z = stats::qnorm(S8ESp)
  S8r.NP = S8p.Z / sqrt(N.Participants)
  S8r.NO = S8p.Z / sqrt(N.Obs)


  GH2015extra = data.frame(rbind(
    p = S8ESp,
    Z = signif(S8p.Z, 4),
    r.NP = signif(S8r.NP, 4),
    r.NO = signif(S8r.NO, 4)
  ))
  #newColumnStatistic = c("p", "Z", "r_pb(NP)", "r_pb(NO)")
  newColumnStatistic = c("p", "Z", "$r_{pb}(NP)$", "$r_{pb}(NO)$")
  GH2015extra = cbind(newColumnStatistic, GH2015extra)
  GH2015extra = reshape::rename(
    GH2015extra,
    c(
      "newColumnStatistic" = "Statistic",
      "X1" = "Exp1",
      "X2" = "Exp2",
      "X3" = "Exp3",
      "X4" = "Exp4"
    )
  )


  #*************************************************************
  # Study 9: S9. Fernandez-Saez 2015. The experiments in this family are all between groups. We could use the correct Hedge's g variance from the reported data. However, for consistency we use the transformation approach. We assume FD is the control and RE is the treatment.

  S9data = subset(ExpData, ExpData$Study == "S9")

  S9EffectSizes = PrepareForMetaAnalysisGtoR(S9data$Mc,
                     S9data$Mt,
                     S9data$SDc,
                     S9data$SDt,
                     S9data$Nc,
                     S9data$Nt)

 # Extract the Hg values and restrict to 4 significant figures
  S9.ES=signif(S9EffectSizes$HGvalues.HgAdjusted,4)

# Do a random effect meta-analysis of the transformed r effect sizes
  S9MA = metafor::rma(S9EffectSizes$zr, S9EffectSizes$vi)


  S9MAStats=ExtractMAStatistics(S9MA,sum(S9data$Nc),sum(S9data$Nc))

  Study = "S9"

  Type = "gIG"
  Source = "Calc"
  MAStats = rbind(MAStats, cbind(Study, Type, Source, S9MAStats, RR = NA))  # bug fixed by LM 2018.09.22/24
  # Index correct table entry
  S9MARep = subset(MARepResShort, MARepResShort$Study == "S9")
  S9MARep$Type="dIG"
  MAReproduced = ifelse (abs(S9MAStats$mean - S9MARep$mean) > 0.05, "No", "Yes")
  MAStats = rbind(MAStats, cbind(S9MARep, RR = MAReproduced))

  ESdata1 = data.frame(
    Study = Study,
    Type=Type,
    Source = "Calc",
    Design = "IndGroups",
    Exp1 = S9.ES[1],
    Exp2 = S9.ES[2],
    Exp3 = S9.ES[3],
    Exp4 = S9.ES[4],
    Exp5 = S9.ES[5],
    RR = NA
  )

  # Extract reported effect sizes

  S9ESRep = subset(ExperimentEffectSizes, ExperimentEffectSizes$Study ==
                     "S9")
  S9ESRep$Type="dIG"
  ESdata1 = rbind(ESdata1, cbind(S9ESRep, RR = NA))

  Diff1 = signif(abs(as.numeric(ESdata1[1, ]$Exp1) - as.numeric(ESdata1[2, ]$Exp1)), 4)
  Diff2 = signif(abs(as.numeric(ESdata1[1, ]$Exp2) - as.numeric(ESdata1[2, ]$Exp2)), 4)
  Diff3 = signif(abs(as.numeric(ESdata1[1, ]$Exp3) - as.numeric(ESdata1[2, ]$Exp3)), 4)


  ESReproduce = 0
  ESReproduce = ESReproduce + ifelse(Diff1 > 0.05, 0, 1)
  ESReproduce = ESReproduce + ifelse(Diff2 > 0.05, 0, 1)
  ESReproduce = ESReproduce + ifelse(Diff3 > 0.05, 0, 1)

  ESdata1[2, ]$RR = ESReproduce

  ESdata = rbind(ESdata, ESdata1)

  #*************************************************************
  # Study S10: S10. Cruz-Lemus 2011. Reported resutls from 3 experiments. All experiments used 4 group crossover design. They did not report effect size values for individual experiments.


  S10data = subset(ExpData, ExpData$Study == "S10")

  S10.EffectSizes = PrepareForMetaAnalysisGtoR(S10data$Mc,
                                               S10data$Mt,
                                               S10data$SDc,
                                               S10data$SDt,
                                               S10data$Nc,
                                               S10data$Nt)

  S10MA.res = metafor::rma(S10.EffectSizes$zr, S10.EffectSizes$vi)

  S10MAStats = ExtractMAStatistics(S10MA.res, sum(S10data$Nc), sum(S10data$Nt))

  Study = "S10"

  S10.ES = signif(S10.EffectSizes$HGvalues.HgAdjusted, 4)

  Type = "gIG"
  Source = "Calc"
  MAStats = rbind(MAStats, cbind(Study, Type, Source, S10MAStats, RR = NA))
  # Index correct table entry
  S10MARep = subset(MARepResShort, MARepResShort$Study == "S10")
  S10MARep$Type="dIG"
  MAReproduced = ifelse (abs(S10MAStats$mean - S10MARep$mean) > 0.05, "No", "Yes")

  MAStats = rbind(MAStats, cbind(S10MARep, RR = MAReproduced))

  ESdata1 = data.frame(
    Study = Study,
    Type = "gIG",
    Source = "Calc",
    Design = "4GroupCO",
    Exp1 = S10.ES[1],
    Exp2 = S10.ES[2],
    Exp3 = S10.ES[3],
    Exp4 = NA,
    Exp5 = NA,
    RR = NA
  )

  # Crux-Lemus did not report experiment Hg values.

  ESdata = rbind(ESdata, ESdata1)

  #*************************************************************
  # Morales-2016 Morales et al. calculated effect sizes for each application. The data averaged over the difference applications are to be found in DescriptiveStatistics.xlsx sheet Morales 2016. The values below were extracted from DescriptiveStatistics.xlsx sheet BasicStatistics, row Morales-2016. The analysis is the same as it was for Tuerel.

  S11data = subset(ExpData, ExpData$Study == "S11")

  Doc1data = subset(DocData, DocData$Study == "S11" &
                      DocData$Doc == "Doc1")

  if (length(S11data$data) < 3) {
    Doc2data = subset(DocData, DocData$Study == "S11" &
                        DocData$Doc == "Doc2")

    c.Data = aggregateIndividualDocumentStatistics(Doc1data$Mc,
                                                   Doc1data$SDc,
                                                   Doc1data$Nc,
                                                   Doc2data$Mc,
                                                   Doc2data$SDc,
                                                   Doc2data$Nc)
    t.Data = aggregateIndividualDocumentStatistics(Doc1data$Mt,
                                                   Doc1data$SDt,
                                                   Doc1data$Nt,
                                                   Doc2data$Mt,
                                                   Doc2data$SDt,
                                                   Doc2data$Nt)

    Study = rep("S11", 3)
    Exp3 = c("Exp1", "Exp2", "Exp3")
    r3 = rep(NA, 3)
    S11data = data.frame(
      Study = Study,
      Exp = Exp3,
      Mc = signif(c.Data$M, 4),
      SDc = signif(c.Data$SD, 6),
      Nc = Doc1data$Nc,
      Mt = signif(t.Data$M, 4),
      SDt = signif(t.Data$SD, 6),
      Nt = Doc1data$Nt,
      r = r3
    )

    ExpData = rbind(ExpData, S11data)

  }

  S11.EffectSizes = PrepareForMetaAnalysisGtoR(S11data$Mc,
                                               S11data$Mt,
                                               S11data$SDc,
                                               S11data$SDt,
                                               S11data$Nc,
                                               S11data$Nt)

  S11MA.res = metafor::rma(S11.EffectSizes$zr, S11.EffectSizes$vi)

  S11MAStats = ExtractMAStatistics(S11MA.res, sum(S11data$Nc), sum(S11data$Nt))

  Study = "S11"

  S11.ES = signif(S11.EffectSizes$HGvalues.HgAdjusted, 4)

  Type = "gIG"
  Source = "Calc"
  MAStats = rbind(MAStats, cbind(Study, Type, Source, S11MAStats, RR = NA))
  # Index correct table entry
  S11MARep = subset(MARepResShort, MARepResShort$Study == "S11")
  S11MARep$Type="dIG"
  MAReproduced = ifelse (abs(S11MAStats$mean - S11MARep$mean) > 0.05, "No", "Yes")

  MAStats = rbind(MAStats, cbind(S11MARep, RR = MAReproduced))

  ESdata1 = data.frame(
    Study = Study,
    Type = Type,
    Source = "Calc",
    Design = "ABBACO",
    Exp1 = S11.ES[1],
    Exp2 = S11.ES[2],
    Exp3 = S11.ES[3],
    Exp4 = S11.ES[4],
    Exp5 = S11.ES[5],
    RR = NA
  )


  # No experiment level effect sizes were reported, so none are reported

  ESdata = rbind(ESdata, ESdata1)

  # Do analysis of Drone - single time period

  S11Doc1.EffectSizes = PrepareForMetaAnalysisGtoR(Doc1data$Mc,
                                                   Doc1data$Mt,
                                                   Doc1data$SDc,
                                                   Doc1data$SDt,
                                                   Doc1data$Nc,
                                                   Doc1data$Nt)

  S11Doc1.c = calculateSmallSampleSizeAdjustment(Doc1data$Nc + Doc1data$Nt -
                                                   2)
  S11Doc1.Vd = S11Doc1.c ^ 2 * ((Doc1data$Nc + Doc1data$Nt) / (Doc1data$Nc *
                                                                 Doc1data$Nt) + S11Doc1.EffectSizes$HGvalues.HgAdjusted ^ 2 / (2 * (Doc1data$Nc +
                                                                                                                                      Doc1data$Nt))
  )
  #sqrt(Drone.Vd)
  # Morales et al.report the standardized mean difference adjusted for small sample sizes in Table 19. They report an analysis based on the unadjusted standardized effect size d_{IG} and the medium sample size variance of d_{IG}

  S11Doc1.VdUnadjusted = ((Doc1data$Nc + Doc1data$Nt) / (Doc1data$Nc * Doc1data$Nt) +
                            S11Doc1.EffectSizes$HGvalues.HgAdjusted ^ 2 / (2 * S11Doc1.c ^ 2 * (Doc1data$Nc +
                                                                                                  Doc1data$Nt))
  )

  S11Doc1MAres = metafor::rma(S11Doc1.EffectSizes$HGvalues.Hg,
                              S11Doc1.VdUnadjusted,
                              method = "ML")
  S11MATP1Stats = ExtractMAStatistics(
    S11Doc1MAres,
    sum(Doc1data$Nc),
    sum(Doc1data$Nt),
    Transform = F,
    sig = 5
  )

  MAStatsTP1 = rbind(MAStatsTP1,
                     cbind(Study, Type = "d", Source = "Calc", S11MATP1Stats))

  MAStatsTP1rep = subset(ABBAMARepResShort, ABBAMARepResShort$Study == "S11")

  MAStatsTP1 = rbind(MAStatsTP1, MAStatsTP1rep)

  S11TP1.Hg = signif(S11Doc1.EffectSizes$HGvalues.Hg, 5)

  #ESTP1res=rbind(ESTP1res,cbind(Study,Type="d",Source="Calc",Design="ABBACO",Exp1=S11TP1.Hg[1],Exp2=S11TP1.Hg[2],Exp3=S11TP1.Hg[3],Exp4=NA,Exp5=NA))
  ESTP1res = rbind(
    ESTP1res,
    data.frame(
      Study,
      Type = "d",
      Source = "Calc",
      Design = "ABBACO",
      Exp1 = S11TP1.Hg[1],
      Exp2 = S11TP1.Hg[2],
      Exp3 = S11TP1.Hg[3],
      Exp4 = NA,
      Exp5 = NA
    )
  ) #issue: character instead of numeric values, e.g., in Exp2, solution: cbind -> data.frame /LM

  S11ESRep = subset(ABBAEffectSizes, ABBAEffectSizes$Study == "S11")

  ESTP1res = rbind(ESTP1res, S11ESRep)



  #*************************************************************
  # Study 13: S13. Laitenberger 2001 Reported the correct formula for the individual improvement effect size. They reported the correlation between repeated. Their meta-analysis aggregated p-values using the correct formula. I suggest we do a meta-analysis based the repeated measures effect size and the correct variance. Also since they reported r, we calculate the standardized mean difference effect size equivalent to a between groups study and do a meta-analysis of that effect size to demonstrate the difference.
  # This analysis is based on the values estimated from the boxplots

  S13data = subset(ExpData, ExpData$Study == "S13")

  # The study aggregates the p-values of the gRM estimates. We need to compare an aggregation of the gRM estimates to compare the overall p-values

  S13.EffectSizes = PrepareForMetaAnalysisGtoR(S13data$Mc,
                                               S13data$Mt,
                                               S13data$SDc,
                                               S13data$SDt,
                                               S13data$Nc,
                                               S13data$Nt)

  S13.dIG=S13.EffectSizes$HGvalues.Hg
  S13.gIG=S13.EffectSizes$HGvalues.HgAdjusted
  S13.gRM=S13.gIG/sqrt(1-S13data$r)

  #The following line was amended to use newly renamed function 29-08-2018
  S13.rPBS = transformHgtoR(S13.gRM, S13data$Nc, S13data$Nt)

  #The following line was amended to use newly renamed function 29-08-2018
  S13RM.z = transformRtoZr(S13.rPBS)

  S13RM.vi = 1 / (S13data$Nc + S13data$Nt-3)

  S13RMres = metafor::rma(S13RM.z, S13RM.vi)

  S13RMMAStats = ExtractMAStatistics(S13RMres, sum(S13data$Nc), sum(S13data$Nt))

  Study = "S13"

  S13RM.ES = signif(S13.gRM, 4)
  S13.dIG=signif(S13.dIG,4)
  S13.gIG=signif(S13.gIG,4)


  Type = "gRM"
  Source = "Calc"
  MAStats = rbind(MAStats, cbind(Study, Type, Source, S13RMMAStats, RR =
                                   NA))
  # Index correct table entry Note there is no study 12, entry 12 applies to study12
  S13RMMArep = subset(MARepResShort,
                      MARepResShort$Study == "S13" & MARepResShort$Type == "gRM")

  MAReproduced = ifelse (abs(S13RMMAStats$mean - S13RMMArep$mean) > 0.05, "No", "Yes")

  MAStats = rbind(MAStats, cbind(S13RMMArep, RR = MAReproduced))

  # Add the gIG and gRM effect sizes


  ESdata2 = data.frame(
    Study = Study,
    Type = "gRM",
    Source = "Calc",
    Design = "ABBACO",
    Exp1 = S13RM.ES[1],
    Exp2 = S13RM.ES[2],
    Exp3 = S13RM.ES[3],
    Exp4 = NA,
    Exp5 = NA,
    RR = NA
  )

  S13RMEffectSizes = subset(ExperimentEffectSizes,
                            ExperimentEffectSizes$Study == "S13" & Type == "gRM")
  # Index correct table entry
  ESdata2 = data.frame(rbind(ESdata2, cbind(S13RMEffectSizes, RR = NA)))
  Diff1 = signif(abs(as.numeric(ESdata2[1, ]$Exp1) - as.numeric(ESdata2[2, ]$Exp1)), 4)
  Diff2 = signif(abs(as.numeric(ESdata2[1, ]$Exp2) - as.numeric(ESdata2[2, ]$Exp2)), 4)
  Diff3 = signif(abs(as.numeric(ESdata2[1, ]$Exp3) - as.numeric(ESdata2[2, ]$Exp3)), 4)

  ESReproduce = 0
  ESReproduce = ESReproduce + ifelse(Diff1 > 0.05, 0, 1)
  ESReproduce = ESReproduce + ifelse(Diff2 > 0.05, 0, 1)
  ESReproduce = ESReproduce + ifelse(Diff3 > 0.05, 0, 1)

  ESdata2[2, ]$RR = ESReproduce

  ESdata = rbind(ESdata, ESdata2)

  ESdata1 = data.frame(
    Study = Study,
    Type = "gIG",
    Source = "Calc",
    Design = "ABBACO",
    Exp1 = S13.gIG[1],
    Exp2 = S13.gIG[2],
    Exp3 = S13.gIG[3],
    Exp4 = NA,
    Exp5 = NA,
    RR = NA
  )

  S13EffectSizes = subset(ExperimentEffectSizes,
                          ExperimentEffectSizes$Study == "S13" & Type == "gIG")

  ESdata1 = data.frame(rbind(ESdata1, cbind(S13EffectSizes, RR = NA)))
  Diff1 = signif(abs(as.numeric(ESdata1[1, ]$Exp1) - as.numeric(ESdata1[2, ]$Exp1)), 4)
  Diff2 = signif(abs(as.numeric(ESdata1[1, ]$Exp2) - as.numeric(ESdata1[2, ]$Exp2)), 4)
  Diff3 = signif(abs(as.numeric(ESdata1[1, ]$Exp3) - as.numeric(ESdata1[2, ]$Exp3)), 4)

  ESReproduce = 0
  ESReproduce = ESReproduce + ifelse(Diff1 > 0.05, 0, 1)
  ESReproduce = ESReproduce + ifelse(Diff2 > 0.05, 0, 1)
  ESReproduce = ESReproduce + ifelse(Diff3 > 0.05, 0, 1)

  ESdata1[2, ]$RR = ESReproduce

  ESdata = rbind(ESdata, ESdata1)


  # Laitenberger also aggregated p-values so we need to add those values to effect sizes and calculate aggregated value P

  # We expect t and gRM to be related so we construct the t-values from that linking equation
  S13.tval = S13RM.ES / sqrt(2 / (S13data$Nc + S13data$Nt))
  #S13.tval = S13RM.gRM / sqrt(2 / (S13data$Nc + S13data$Nt)) #LM: bug reported 20190427
  S13.pt = 1 - stats::pt(S13.tval, S13data$Nc + S13data$Nt - 2)
  S13.P = -2 * (sum(log(S13.pt)))
  # P is distributed as chi-squared with 2k degrees of freedom where k is the numberof independent experiments.
  S13.P.pvalue = 1 - stats::pchisq(S13.P, 6)

  #MAStats=rbind(MAStats,cbind(Study,Type="P",Source="Calc",mean=signif(S13.P,4),pvalue=signif(S13.P.pvalue,4),UB=NA,LB=NA,QE=NA,QEp=NA,RR=NA))
  MAStats = rbind(
    MAStats,
    data.frame(
      Study,
      Type = "P",
      Source = "Calc",
      mean = signif(S13.P, 4),
      pvalue = signif(S13.P.pvalue, 4),
      UB = NA,
      LB = NA,
      QE = NA,
      QEp = NA,
      RR = NA
    )
  ) #issue: character instead of numeric values, solution: cbind -> data.frame /LM
  # Index correct table entry

  # There is no standard size for P so MAReproduced is set to null
  MAReproduced = NA
  S13PRep = subset(MARepResShort,
                   MARepResShort$Study == "S13" & MARepResShort$Type == "P")
  MAStats = rbind(MAStats, cbind(S13PRep, RR = MAReproduced))

  S13.pt = signif(S13.pt, 4)

  ESdata3 = data.frame(
    Study = Study,
    Type = "p",
    Source = "Calc",
    Design = "ABBACO",
    Exp1 = S13.pt[1],
    Exp2 = S13.pt[2],
    Exp3 = S13.pt[3],
    Exp4 = NA,
    Exp5 = NA,
    RR = NA
  )

  S13pEffectSizes = subset(ExperimentEffectSizes,
                           ExperimentEffectSizes$Study == "S13" & Type == "p")
  ESdata3 = data.frame(rbind(ESdata3, cbind(S13pEffectSizes, RR = NA)))
  Diff1 = signif(abs(as.numeric(ESdata3[1, ]$Exp1) - as.numeric(ESdata3[2, ]$Exp1)), 4)
  Diff2 = signif(abs(as.numeric(ESdata3[1, ]$Exp2) - as.numeric(ESdata3[2, ]$Exp2)), 4)
  Diff3 = signif(abs(as.numeric(ESdata3[1, ]$Exp3) - as.numeric(ESdata3[2, ]$Exp3)), 4)
  ESReproduce = 0
  ESReproduce = ESReproduce + ifelse(Diff1 > 0.05, 0, 1)
  ESReproduce = ESReproduce + ifelse(Diff2 > 0.05, 0, 1)
  ESReproduce = ESReproduce + ifelse(Diff3 > 0.05, 0, 1)

  ESdata3[2, ]$RR = ESReproduce
  ESdata = rbind(ESdata, ESdata3)

  #*************************************************************
  # Study 14: S14. Pfahl et al. 20054. Pfahl aggregated p-values but also used values based on differences of differences. Revised analysis to be based on the reported differences to reduce problems with rounding errors.

  S14data = subset(ExpData, ExpData$Study == "S14")
  S14EffectSizes = PrepareForMetaAnalysisGtoR(S14data$Mc,
                                              S14data$Mt,
                                              S14data$SDc,
                                              S14data$SDt,
                                              S14data$Nc,
                                              S14data$Nt)

S14MA=metafor::rma(S14EffectSizes$zr,S14EffectSizes$vi)

S14MAStats=ExtractMAStatistics(S14MA,sum(S14data$Nc),sum(S14data$Nt))

# Do special update to tables for Pfhal
Study="S14"

#Phfal aggregated both the dRM value and the p-values, so we calculate both and include both in the meta-analysis table

Type="gRM"
MAStats=rbind(MAStats,cbind(Study,Type,Source="Calc",S14MAStats,RR=NA))

S14MArep=subset(MARepResShort,MARepResShort$Study=="S14"&MARepResShort$Type=="dRM")

MAReproduced= ifelse (abs(S14MAStats$mean-S14MArep$mean) > 0.05,"No","Yes")

MAStats=rbind(MAStats,cbind(S14MArep,RR=MAReproduced))


  # Aggregate the one-sided p-values from the t-test and put the calculated and reported meta-analysis rsults in MAStats

  S14.tval=S14EffectSizes$HGvalues.HgAdjusted/sqrt(1/S14data$Nc+1/S14data$Nt)

  S14.tp = 1 - stats::pt(S14.tval, S14data$Nc + S14data$Nt - 2)

  S14.P = -2 * sum(log(S14.tp))

  S14.P.pvalue = 1 - stats::pchisq(S14.P, 6)

  # Included the calculated P value resutls into MAStats

  #MAStats=rbind(MAStats,cbind(Study,Type="P",Source="Calc",mean=signif(S14.P,4),pvalue=signif(S14.P.pvalue,4),UB=NA,LB=NA,QE=NA,QEp=NA,RR=NA))
  MAStats = rbind(
    MAStats,
    data.frame(
      Study,
      Type = "P",
      Source = "Calc",
      mean = signif(S14.P, 4),
      pvalue = signif(S14.P.pvalue, 4),
      UB = NA,
      LB = NA,
      QE = NA,
      QEp = NA,
      RR = NA
    )
  ) #issue: character instead of numeric values, e.g., in pvalue, solution: cbind -> data.frame /LM


  # There is no standard size for P values so we do not do a reproducibility test.
  # Put the reported P value resutls into MAStats
  S14MAPrep = subset(MARepResShort,
                     MARepResShort$Study == "S14" & MARepResShort$Type == "P")
  MAStats = rbind(MAStats, cbind(S14MAPrep, RR = NA))

  # Put the calculated dRM effect sizes into a temporary data.file

  S14Hg.ES=signif(S14EffectSizes$HGvalues.HgAdjusted,4)


  ESdata1 = data.frame(
    Study = Study,
    Type = "dRM",
    Source = "Calc",
    Design = "PrePost",
    Exp1 = S14Hg.ES[1],
    Exp2 = S14Hg.ES[2],
    Exp3 = S14Hg.ES[3],
    Exp4 = NA,
    Exp5 = NA,
    RR = NA
  )

  # Obtain the reported dRM effect sizes
  S14HgRep = subset(
    ExperimentEffectSizes,
    ExperimentEffectSizes$Study == "S14" &
      ExperimentEffectSizes$Type == "dRM"
  )

  # Store the report dRM effect sizes in the temporary data.file
  ESdata1 = data.frame(rbind(ESdata1, cbind(S14HgRep, RR = NA)))

  # Calculated the reproducibility of the dRM effect sizes
  Diff1 = signif(abs(as.numeric(ESdata1[1, ]$Exp1) - as.numeric(ESdata1[2, ]$Exp1)), 4)
  Diff2 = signif(abs(as.numeric(ESdata1[1, ]$Exp2) - as.numeric(ESdata1[2, ]$Exp2)), 4)
  Diff3 = signif(abs(as.numeric(ESdata1[1, ]$Exp3) - as.numeric(ESdata1[2, ]$Exp3)), 4)
  ESReproduce = 0
  ESReproduce = ESReproduce + ifelse(Diff1 > 0.05, 0, 1)
  ESReproduce = ESReproduce + ifelse(Diff2 > 0.05, 0, 1)
  ESReproduce = ESReproduce + ifelse(Diff3 > 0.05, 0, 1)

  ESdata1[2, ]$RR = ESReproduce
  # Add the calculated and reported effect size results to ESdata

  ESdata = rbind(ESdata, ESdata1)

  # Obtain the calculated p values
  S14pRep = subset(
    ExperimentEffectSizes,
    ExperimentEffectSizes$Study == "S14" &
      ExperimentEffectSizes$Type == "p"
  )
  # Put the calculated p values into a temporary data.fromae
  ESdata1 = data.frame(
    Study = Study,
    Type = "p",
    Source = "Calc",
    Design = "PrePost",
    Exp1 = signif(S14.tp[1], 4),
    Exp2 = signif(S14.tp[2]),
    Exp3 = signif(S14.tp[3], 4),
    Exp4 = NA,
    Exp5 = NA,
    RR = NA
  )
  # Add the calculated p values to the temporary data.frame
  ESdata1 = data.frame(rbind(ESdata1, cbind(S14pRep, RR = NA)))
  # Assess the reproducibility of the calculated values
  Diff1 = signif(abs(as.numeric(ESdata1[1, ]$Exp1) - as.numeric(ESdata1[2, ]$Exp1)), 4)
  Diff2 = signif(abs(as.numeric(ESdata1[1, ]$Exp2) - as.numeric(ESdata1[2, ]$Exp2)), 4)
  Diff3 = signif(abs(as.numeric(ESdata1[1, ]$Exp3) - as.numeric(ESdata1[2, ]$Exp3)), 4)

  ESReproduce = 0
  ESReproduce = ESReproduce + ifelse(Diff1 > 0.05, 0, 1)
  ESReproduce = ESReproduce + ifelse(Diff2 > 0.05, 0, 1)
  ESReproduce = ESReproduce + ifelse(Diff3 > 0.05, 0, 1)

  #Add the reported and calculate effect size data to ESdata
  ESdata1[2, ]$RR = ESReproduce

  ESdata = rbind(ESdata, ESdata1)
  ####################################################################
  ####################################################################
  # Report main meta-analysis table & main effect size comparison table

  MAStats = data.frame(MAStats, row.names = NULL)
  ESdata = data.frame(ESdata, row.names = NULL)

  #library("xtable")
  #xtable(MAStats)
  #xtable(ESdata)

  # Report values for 3 papers that reported per document

  MAStatsTP1 = data.frame(MAStatsTP1, row.names = NULL)
  ESTP1res = data.frame(ESTP1res, row.names = NULL)

  #xtable(MAStatsTP1)
  #xtable(ESTP1res)

  # Report extra results for Study 8
  #xtable(GH2015extra)

  results = list(
    MAStats = MAStats,
    ESdata = ESdata,
    MAStatsTP1 = MAStatsTP1,
    ESTP1res = ESTP1res,
    GH2015extra = GH2015extra
  )
  return(results)

}
