#' @title reproduceTableWithSourceDataByCiolkowski
#' @description Function reproduces Table, which shows the effect sizes reported by Ciolkowski identifying the type of design used in each study.
#' @author Lech Madeyski
#' @export reproduceTableWithSourceDataByCiolkowski
#' @examples
#' reproduceTableWithSourceDataByCiolkowski()
reproduceTableWithSourceDataByCiolkowski <- function(){
  reproducer::printXTable(data = reproducer::Ciolkowski09ESEM.MetaAnalysis.PBRvsCBRorAR, selectedColumns = 'c("Study", "Ref.", "Control", "Within-subjects", "Cross-over", "d_ByCiolkowski", "d_ByOriginalAuthors")', tableType = "latex", alignCells = "c(rep('l',6), rep('d{1.4}',1), rep('l',1))", digits = c(0,0,0,0,0,0,4,2), caption = "Source Data", label = "tab:SourceDataUsedByCiolkowski", fontSize = "footnotesize", captionPlacement = "top", alignHeader = c( "\\multicolumn{1}{l}{Study}", "\\multicolumn{1}{l}{Ref.}", "\\multicolumn{1}{l}{Control}", "\\multicolumn{1}{p{1.2cm}}{Within-subjects}", "\\multicolumn{1}{p{1.0cm}}{Cross-over}", "\\multicolumn{1}{p{1.9cm}}{$d$ (by Ciol-kowski)}", "\\multicolumn{1}{p{2.1cm}}{$d$ (by origi-nal authors)}" ) )
}

#' @title reproduceTableWithEffectSizesBasedOnMeanDifferences()
#' @description Function reproduces Table, which shows the effect sizes based on mean differences.
#' @author Lech Madeyski
#' @export reproduceTableWithEffectSizesBasedOnMeanDifferences
#' @examples
#' reproduceTableWithEffectSizesBasedOnMeanDifferences()
reproduceTableWithEffectSizesBasedOnMeanDifferences <- function(){
  dataMK <- reproducer::MadeyskiKitchenham.MetaAnalysis.PBRvsCBRorAR
  rownames(dataMK)<-dataMK$Study
  printXTable(data=dataMK, selectedColumns='c("Study", "Ref.", "ExpDesign", "Teams", "n_PBR", "n_C", "M_PBR", "M_C", "Diff", "Inc", "V_C", "V_D")', tableType="latex", alignCells="c(rep('l',4), rep('d{2.0}',3), rep('d{1.3}',2), rep('d{2.3}',1), rep('d{3.0}',1), rep('d{1.3}',2))", digits=c(0,0,0,0,0,0,0,3,3,3,0,3,3), caption="Effect Sizes based on Mean Differences", label="tab:MDEffectSize", fontSize="footnotesize", captionPlacement="top", alignHeader=c("\\multicolumn{1}{l}{Study}", "\\multicolumn{1}{l}{Ref.}", "\\multicolumn{1}{l}{ExpDesign}", "\\multicolumn{1}{l}{$Teams$}", "\\multicolumn{1}{l}{$n_{PBR}$}", "\\multicolumn{1}{l}{$n_C$}", "\\multicolumn{1}{l}{$M_{PBR}$}", "\\multicolumn{1}{l}{$M_C$}", "\\multicolumn{1}{l}{$\\overline{\\mathit{diff}}$}", "\\multicolumn{1}{l}{$\\%Inc$}", "\\multicolumn{1}{l}{$V_C$}", "\\multicolumn{1}{l}{$V_D$}"))
}

#' @title reproduceForestPlotRandomEffects()
#' @description Function reproduces Forest Plot of a Random-Effects Meta-analysis of Mean Differences.
#' @author Lech Madeyski
#' @export reproduceForestPlotRandomEffects
#' @examples
#' reproduceForestPlotRandomEffects()
reproduceForestPlotRandomEffects <- function(){
  dataMK <- reproducer::MadeyskiKitchenham.MetaAnalysis.PBRvsCBRorAR
  res<-metafor::rma(yi=dataMK$Diff,vi=dataMK$V_D,data=dataMK,method="REML")
  metafor::forest(res, xlab="Overall Effect Size", slab=reproducer::MadeyskiKitchenham.MetaAnalysis.PBRvsCBRorAR$Study)
  graphics::text(-1,19," Study")
  graphics::text(0.9,19, "Effect Size   [95% CL]")
}

#' @title reproduceMixedEffectsAnalysisWithExperimentalDesignModerator()
#' @description Function reproduces Mixed-Effects Analysis with Experimental Design as a Moderator.
#' @author Lech Madeyski
#' @export reproduceMixedEffectsAnalysisWithExperimentalDesignModerator
#' @examples
#' reproduceMixedEffectsAnalysisWithExperimentalDesignModerator()
reproduceMixedEffectsAnalysisWithExperimentalDesignModerator <- function(){
  dataMK <- reproducer::MadeyskiKitchenham.MetaAnalysis.PBRvsCBRorAR
  resWithExpDesignMod<-metafor::rma(yi=dataMK$Diff, vi=dataMK$V_D, data=dataMK, method="REML", mods=~factor(ExpDesign)-1)
  resWithExpDesignMod
}

#' @title reproduceMixedEffectsForestPlotWithExperimentalDesignModerator()
#' @description Function reproduces Forest Plot of a Mixed Effects Meta-analysis of Mean Differences with Experimental Design as a Moderator Variable.
#' @author Lech Madeyski
#' @export reproduceMixedEffectsForestPlotWithExperimentalDesignModerator
#' @examples
#' reproduceMixedEffectsForestPlotWithExperimentalDesignModerator()
reproduceMixedEffectsForestPlotWithExperimentalDesignModerator <- function(){
  dataMK <- reproducer::MadeyskiKitchenham.MetaAnalysis.PBRvsCBRorAR
  resWithExpDesignMod<-metafor::rma(yi=dataMK$Diff, vi=dataMK$V_D, data=dataMK, method="REML", mods=~factor(ExpDesign)-1)
  metafor::forest(resWithExpDesignMod,xlab="Overall Effect Size", slab=reproducer::MadeyskiKitchenham.MetaAnalysis.PBRvsCBRorAR$Study)
  graphics::text(-1,19," Study")
  graphics::text(0.9,19, "Effect Size   [95% CL]")
}


#' @title reproduceTableWithPossibleModeratingFactors()
#' @description Function reproduces Table with possible moderating factors.
#' @author Lech Madeyski
#' @export reproduceTableWithPossibleModeratingFactors
#' @examples
#' reproduceTableWithPossibleModeratingFactors()
reproduceTableWithPossibleModeratingFactors <- function(){
  dataMK <- reproducer::MadeyskiKitchenham.MetaAnalysis.PBRvsCBRorAR
  rownames(dataMK) <- dataMK$Study
  printXTable(data=dataMK, selectedColumns='c("Study", "Ref.", "ControlType", "ParticipantsType", "TeamType", "ArtefactType", "AssociatedWithBasili")', tableType="latex", alignCells="c(rep('c',8))", digits=0, caption="Possible Moderating Factors", label="tab:Moderators", fontSize="footnotesize", captionPlacement="top", alignHeader=c("\\multicolumn{1}{c}{Study}", "\\multicolumn{1}{c}{Ref.}", "\\multicolumn{1}{c}{Control}","\\multicolumn{1}{c}{Participants}", "\\multicolumn{1}{c}{TeamType}", "\\multicolumn{1}{c}{ArtefactType}", "\\multicolumn{1}{c}{AssociatedWithBasili}"))
}


#' @title reproduceMixedEffectsAnalysisWithEstimatedVarianceAndExperimentalDesignModerator()
#' @description Function reproduces Mixed-Effects Analysis using Subject Specific Estimated Variance with Experimental Design as a Moderator.
#' @author Lech Madeyski
#' @export reproduceMixedEffectsAnalysisWithEstimatedVarianceAndExperimentalDesignModerator
#' @examples
#' reproduceMixedEffectsAnalysisWithEstimatedVarianceAndExperimentalDesignModerator()
reproduceMixedEffectsAnalysisWithEstimatedVarianceAndExperimentalDesignModerator <- function(){
  dataMK <- reproducer::MadeyskiKitchenham.MetaAnalysis.PBRvsCBRorAR
  resWithExpDesignUsingSubjectSpecificVarianceMod<-metafor::rma(yi=dataMK$Diff, vi=dataMK$V_Alt, data=dataMK,method="REML", mods=~factor(ExpDesign)-1)
  resWithExpDesignUsingSubjectSpecificVarianceMod
}


#' -------------------------------------------------------------------------------------------------------
#' Functions related to a joint paper with Barbara Kitchenham, "Effect sizes and their variance for AB/BA crossover design studies"
#' -------------------------------------------------------------------------------------------------------
#' @title getSimulationData
#' @description Function to generate the simulated data set used in a paper "Effect Sizes and their Variance for AB/BA Crossover Design Studies" by Lech Madeyski and Barbara Kitchenham
#' @author Lech Madeyski and Barbara Kitchenham
#' @export getSimulationData
#' @param var Variance among subjects is a sum of the between subjects variance and the within subjects variance
#' @param covar Covariance equal to the between subjects variance
#' @param meanA1 Mean for treatment sequence A1
#' @param treatmentDiff technique effect which is the difference between the effect of technique A and technique B
#' @param periodEffect Period effect which is the difference between period 1 and period 2
#' @param numOfSamples Number of samples ("rows" of data) required for each technique and period
#' @return Data frame:
#' 'data.frame':  4*numOfSamples obs. of  5 variables:
#'  $ pid      : int  1 2 3 4 5 6 7 8 9 10 ...
#'  $ technique: Factor w/ 2 levels "T1","T2":  ...
#'  $ period   : Factor w/ 2 levels "P1","P2":  ...
#'  $ sequence : Factor w/ 2 levels "S1","S2":  ...
#'  $ result   : num  ...
#' @examples
#' data<-getSimulationData(25, 18.75, 50, 10, 5, 500) # generate the simulated data set from the paper
getSimulationData <- function(var, covar, meanA1, treatmentDiff, periodEffect, numOfSamples)
{

  #Covariance matrix

  ##r=.75
  #var=25
  #covar=18.75

  sigmaA=matrix(c(var,covar,covar,var),nrow=2,ncol=2)

  #Mean for treatment sequence A
  #common mean mu=50, treatment difference = 10, period effect=5
  #meanA=c(50,65)
  meanA=c(meanA1,meanA1+treatmentDiff+periodEffect)
  #Mean for sequence B
  #meanB=c(60,55)
  meanB=c(meanA1+treatmentDiff,meanA1+periodEffect)


  #Sample data for sequence A
  set.seed(123)


  Num<-numOfSamples

  mydata1=MASS::mvrnorm(Num, meanA, sigmaA)

  #Sample data for Sequence B

  set.seed(9123)
  mydata2=MASS::mvrnorm(Num, meanB, sigmaA)

  mydata1=as.data.frame(mydata1)
  mydata2=as.data.frame(mydata2)
  names(mydata1)=c("AP1","AP2")
  names(mydata2)=c("BP1","BP2")

  diff1=mydata1$AP2-mydata1$AP1
  diff2=mydata2$BP2-mydata2$BP1

  treat=(mean(diff1)-mean(diff2))/2
  period=(mean(diff1)+mean(diff2))/2

  v1=mydata1$AP1
  v2=mydata1$AP2
  v3=mydata2$BP1
  v4=mydata2$BP2

  AP1=data.frame(pid=seq(1,Num),technique=rep("T1",Num),period=rep("P1",Num),sequence=rep("S1",Num), result=v1)
  AP2=data.frame(pid=seq(1,Num),technique=rep("T2",Num), period=rep("P2",Num),sequence=rep("S1",Num), result=v2)
  allA=rbind(AP1,AP2)
  BP1=data.frame(pid=seq(Num+1,2*Num),technique=rep("T2",Num),period=rep("P1",Num),sequence=rep("S2",Num), result=v3)
  BP2=data.frame(pid=seq(Num+1,2*Num),technique=rep("T1",Num),period=rep("P2",Num),sequence=rep("S2",Num),result=v4)
  allB=rbind(BP1,BP2)
  SimulationData=rbind(allA, allB)

  return(SimulationData)
}


#' @title plotOutcomesForIndividualsInEachSequenceGroup
#' @description Function to plot a figure on the outcomes for individuals in each sequence group used in a paper "Effect Sizes and their Variance for AB/BA Crossover Design Studies" by Lech Madeyski and Barbara Kitchenham
#' @author Lech Madeyski and Barbara Kitchenham
#' @export plotOutcomesForIndividualsInEachSequenceGroup
#' @param var Variance among subjects is a sum of the between subjects variance and the within subjects variance
#' @param covar Covariance equal to the between subjects variance
#' @param meanA1 Mean for treatment sequence A1
#' @param treatmentDiff technique effect which is the difference between the effect of technique A and technique B
#' @param periodEffect Period effect which is the difference between period 1 and period 2
#' @param numOfSamples Number of samples ("rows" of data) required for each technique and period
#' @return plot
#' @examples
#' myPlot<-plotOutcomesForIndividualsInEachSequenceGroup(25, 18.75, 50, 10, 5, 15)
plotOutcomesForIndividualsInEachSequenceGroup <- function(var, covar, meanA1, treatmentDiff, periodEffect, numOfSamples)
{

  #Covariance matrix

  ##r=.75
  #var=25
  #covar=18.75

  sigmaA=matrix(c(var,covar,covar,var),nrow=2,ncol=2)

  #Mean for treatment sequence A
  #common mean mu=50, treatment difference = 10, period effect=5
  #meanA=c(50,65)
  meanA=c(meanA1,meanA1+treatmentDiff+periodEffect)
  #Mean for sequence B
  #meanB=c(60,55)
  meanB=c(meanA1+treatmentDiff,meanA1+periodEffect)


  #Sample data for sequence A
  set.seed(123)

  Num=15
  mydata1=MASS::mvrnorm(Num, meanA, sigmaA)

  #Sample data for Sequence B

  set.seed(9123)
  mydata2=MASS::mvrnorm(Num, meanB, sigmaA)

  mydata1=as.data.frame(mydata1)
  mydata2=as.data.frame(mydata2)
  names(mydata1)=c("AP1","AP2")
  names(mydata2)=c("BP1","BP2")

  diff1=mydata1$AP2-mydata1$AP1
  diff2=mydata2$BP2-mydata2$BP1

  treat=(mean(diff1)-mean(diff2))/2
  period=(mean(diff1)+mean(diff2))/2

  v1=mydata1$AP1
  v2=mydata1$AP2
  v3=mydata2$BP1
  v4=mydata2$BP2

  p1=rep("P1",Num)
  p2=rep("P2",Num)
  t1=rep("T1",Num)
  t2=rep("T2",Num)
  seq1=seq(1,Num)
  seq2=seq(Num+1,2*Num)
  s1=rep("S1",Num)
  s2=rep("S2", Num)

  AP1=data.frame(pid=seq1,period=p1,sequence=s1, T1result=v1)
  BP2=data.frame(pid=seq2,period=p2,sequence=s2, T1result=v4)
  allT1=rbind(AP1,BP2)
  AP2=data.frame(pid2=seq1, period=p2,sequence=s1, T2result=v2)
  BP1=data.frame(pid2=seq2, period=p1,sequence=s2, T2result=v3)

  allT2=rbind(AP2,BP1)

  graphics::par(mfrow=c(2,1),cex=.75)
  graphics::plot(AP1$pid,AP1$T1result,ylim=c(30,80),xlim=c(1,30),ylab="Outcome",xlab="Subject ID",pch=3,
                 main="(a) Raw Data for Each Individual")
  graphics::points(BP2$pid,BP2$T1result,pch=3)
  graphics::points(AP2$pid2,AP2$T2result,pch=0)
  graphics::points(BP1$pid2,BP1$T2result,pch=0)
  graphics::abline(v=15.5)
  graphics::text(0,75,"Technique B before A", pos=4)
  graphics::text(16,75, "Technique A before B", pos=4)
  graphics::legend("bottomright",title="Technique Type", inset=0.05, c("B","A"),pch=c(3,0))

  allA=cbind(AP1,AP2)
  DiffA=allA$T2result-allA$T1result
  allA=cbind(allA,DiffA)
  newA=data.frame(pid=AP1$pid,sequence=AP1$sequence,difference=allA$DiffA)
  allB=cbind(BP1,BP2)
  DiffB=allB$T2result-allB$T1result
  allB=cbind(allB,DiffB)
  newB=data.frame(pid=BP1$pid,sequence=BP1$sequence,difference=allB$DiffB)
  newAll=rbind(newA,newB)
  graphics::boxplot(difference ~ sequence,data=newAll,ylab="Subject Difference for each Sequence", xlab="Sequence ID",main="(b) Crossover Differences for each Sequence")

}

#' @title getEffectSizesABBA
#' @description Function to calculate both effect sizes (dIG, dRM), i.e., independent groups and repeated measures standardized effect sizes and variances, for AB/BA crossover design studies. Function is used in a paper "Effect Sizes and their Variance for AB/BA Crossover Design Studies" by Lech Madeyski and Barbara Kitchenham.
#' @author Lech Madeyski and Barbara Kitchenham
#' @export getEffectSizesABBA
#' @param simulationData - data set in a form required to calculate effect sizes in AB/BA crossover experimental designs
#' @return data frame incl. calculated effect sizes and variances:
#' # dIG - independent groups standardized effect size
#' # var.dIG - variance of independent groups standardized effect size
#' # dRM - repeated measures (within-subjects) standardized effect size
#' # var.dRM - variance of repeated measures (within-subjects) standardized effect size
#' # dIG.Fromt - independent groups standardized effect size calculated from t: dIG.Fromt=t*sqrt(1-r)*sqrt((N1+N2)/(2*N1*N2))
#' # var.dIG.Fromt - variance of independent groups standardized effect size calculated from t: var.dIG.Fromt=var.t*(1-r)*((N1+N2)/(2*N1*N2))
#' # dRM.Fromt - dRM calculated from t: dRM.Fromt=t*sqrt((N1+N2)/(2*N1*N2))
#' # var.dRM.Fromt - var.dRM calculated from t: var.dRM.Fromt = var.t*((N1+N2)/(2*N1*N2))
#' # var.dRM.Fromt2 - var.dRM calculated from t or rather dRM.Fromt: var.dRM.Fromt2=(df/(df-2))*((N1+N2)/(2*N1*N2)+dRM.Fromt^2)- dRM.Fromt^2/c^2
#' # var.dRM.Approx - var.dRM calculated on a basis of Johnson and Welch (1940) report an approximate formulate for the variance of a t variable: var.dRM.Approx=((N1+N2)/(2*N1*N2)) + (dRM^2)/(2*(N1+N2-2)) #see paper and Equation 49
#' # var.dIG.Approx - var.dIG calculated on a basis of Johnson and Welch (1940) report an approximate formulate for the variance of a t variable: var.dIG.Approx=(((N1+N2)*(1-r))/(2*N1*N2)) + (dIG^2)/(2*(N1+N2-2)) #see paper and Equation 50
#' # unstandardizedES - estimated unstandardized technique effect size
#' # periodES - estimated period effect
#' # var.sig - sum of within-subjects variance and between-subjects variance
#' # var.within - within-subjects variance
#' # var.between - between-subjects variance
#' # t - t-value
#' # var.t - variance of t-variable
#' # gRM - Hedges and Olkin (1985) unbiased estimator of the repeated measures effect size gRM=dRM*c
#' # var.gRM - variance of gRM calculated as follows: var.gRM=(df/(df-2))*(((N1+N2)/(2*N1*N2))*c^2+gRM^2)- gRM^2/c^2 #Equation 56
#' # var.gRM2 - variance of gRM calculated as follows: var.gRM2=var.dRM*c^2
#' # gIG - Hedges and Olkin (1985) unbiased estimator of the independent groups effect size gIG=dIG*c
#' # var.gIG - variance of gIG calculated as follows: var.gIG=(df/(df-2))*(((N1+N2)/(2*N1*N2))*c^2+gIG^2)- gIG^2/c^2 #Equation 57
#' # var.gIG2 - variance of gRM calculated as follows: var.gIG2=var.dIG*c^2
#' # r - the correlation between the values observed for the same subject
#' @examples
#' simulationData<-getSimulationData(25, 18.75, 50, 10, 5, 500) #generate simulated data set
#' es<-getEffectSizesABBA(simulationData) #return effect sizes and variances
getEffectSizesABBA <- function(simulationData)
{
  # Replaced label "subjects" with "pid"
  cofit.re=lme4::lmer(result ~ technique+period+(1|pid),data=simulationData)
  re<-lme4::VarCorr(cofit.re)
  re.df<-data.frame(re)
  var.sig<-re.df$vcov[1]+re.df$vcov[2]
  var.within= re.df$vcov[2]
  var.between=re.df$vcov[1]

  fe<-lme4::fixef(cofit.re)
  fe.df<-data.frame(fe)
  unstandardizedES<-fe.df$fe[2]
  periodES <- fe.df$fe[3]

  t<-stats::coef(summary(cofit.re))[2,"t value"] #tRM
  dIG=unstandardizedES/sqrt(var.sig)
  dRM=unstandardizedES/sqrt(var.within)
  r=var.between/var.sig
  N1=N2=length(simulationData$pid)/4
  df=N1+N2-2
  c=1-3/(4*df-1)
  var.t=(df/(df-2))*(1+t^2)- t^2/c^2
  var.dRM.Fromt = var.t*((N1+N2)/(2*N1*N2)) #var.t/N1
  var.dRM=(df/(df-2))*((N1+N2)/(2*N1*N2)+dRM^2)- dRM^2/c^2 #Equation 42
  dRM.Fromt=t*sqrt((N1+N2)/(2*N1*N2))
  var.dRM.Fromt2=(df/(df-2))*((N1+N2)/(2*N1*N2)+dRM.Fromt^2)- dRM.Fromt^2/c^2
  var.dRM.Approx=((N1+N2)/(2*N1*N2)) + (dRM^2)/(2*(N1+N2-2))#see paper and Equation 49
  var.dIG.Approx=(((N1+N2)*(1-r))/(2*N1*N2)) + (dIG^2)/(2*(N1+N2-2))#see paper and Equation 50
  var.dIG=(df/(df-2))*((1-r)*(N1+N2)/(2*N1*N2)+dIG^2)- dIG^2/c^2
  dIG.Fromt=t*sqrt(1-r)*sqrt((N1+N2)/(2*N1*N2))
  var.dIG.Fromt=var.t*(1-r)*((N1+N2)/(2*N1*N2))

  gRM=dRM*c
  var.gRM=(df/(df-2))*(((N1+N2)/(2*N1*N2))*c^2+gRM^2)- gRM^2/c^2 #Equation 56
  var.gRM2=var.dRM*c^2

  gIG=dIG*c
  var.gIG=(df/(df-2))*(((N1+N2)*(1-r)/(2*N1*N2))*c^2+gIG^2)- gIG^2/c^2 #Equation 56
  var.gIG2=var.dIG*c^2

  r=var.between/var.sig

  effectSizes<-data.frame(dIG, var.dIG, dRM, var.dRM, dIG.Fromt, var.dIG.Fromt, dRM.Fromt, var.dRM.Fromt, var.dRM.Fromt2, var.dRM.Approx, var.dIG.Approx, unstandardizedES, periodES, var.sig, var.within, var.between, t, var.t, gRM, var.gRM, var.gRM2, gIG, var.gIG, var.gIG2, r)
  return(effectSizes)
}


#' @title getEffectSizesABBAIgnoringPeriodEffect
#' @description Function to calculate both effect sizes (dIG.ipe, dRM.ipe), i.e., independent groups and repeated measures standardized effect sizes and variances, for AB/BA crossover design studies ignoring period effect (thus wrong). Function is used in a paper "Effect Sizes and their Variance for AB/BA Crossover Design Studies" by Lech Madeyski and Barbara Kitchenham.
#' @author Lech Madeyski and Barbara Kitchenham
#' @export getEffectSizesABBAIgnoringPeriodEffect
#' @param simulationData - data set in a form required to calculate effect sizes in AB/BA crossover experimental designs
#' @return data frame incl. calculated effect sizes and variances:
#' # dIG.ipe - independent groups standardized effect size
#' # var.dIG.ipe - variance of independent groups standardized effect size
#' # dRM.ipe - repeated measures (within-subjects) standardized effect size
#' # var.dRM.ipe - variance of repeated measures (within-subjects) standardized effect size
#' # dIG.Fromt.ipe - independent groups standardized effect size calculated from t: dIG.Fromt=t*sqrt(1-r)*sqrt((N1+N2)/(2*N1*N2))
#' # var.dIG.Fromt.ipe - variance of independent groups standardized effect size calculated from t: var.dIG.Fromt=var.t*(1-r)*((N1+N2)/(2*N1*N2))
#' # dRM.Fromt.ipe - dRM calculated from t: dRM.Fromt=t*sqrt((N1+N2)/(2*N1*N2))
#' # var.dRM.Fromt.ipe - var.dRM calculated from t: var.dRM.Fromt = var.t*((N1+N2)/(2*N1*N2))
#' # var.dRM.Fromt2.ipe - var.dRM calculated from t or rather dRM.Fromt: var.dRM.Fromt2=(df/(df-2))*((N1+N2)/(2*N1*N2)+dRM.Fromt^2)- dRM.Fromt^2/c^2
#' # unstandardizedES.ipe - estimated unstandardized technique effect size
#' # var.sig.ipe - sum of within-subjects variance and between-subjects variance
#' # var.within.ipe - within-subjects variance
#' # var.between.ipe - between-subjects variance
#' # t.ipe - t-value
#' # var.t.ipe - variance of t-variable
#' @examples
#' simulationData<-getSimulationData(25, 18.75, 50, 10, 5, 500) #generate simulated data set
#' es.ipe<-getEffectSizesABBAIgnoringPeriodEffect(simulationData) #return effect sizes and variances
getEffectSizesABBAIgnoringPeriodEffect <- function(simulationData)
{
  # Replaced label "subjects" with "pid", removed period effect
  cofit.re.ipe=lme4::lmer(result ~ technique+(1|pid),data=simulationData)
  re.ipe<-lme4::VarCorr(cofit.re.ipe)
  re.ipe.df<-data.frame(re.ipe)
  var.sig.ipe<-re.ipe.df$vcov[1]+re.ipe.df$vcov[2]
  var.within.ipe<-re.ipe.df$vcov[2]
  var.between.ipe<-re.ipe.df$vcov[1]

  fe.ipe<-lme4::fixef(cofit.re.ipe)
  fe.ipe.df<-data.frame(fe.ipe)
  unstandardizedES.ipe<-fe.ipe.df$fe[2]

  t.ipe<-stats::coef(summary(cofit.re.ipe))[2,"t value"]

  dIG.ipe=unstandardizedES.ipe/sqrt(var.sig.ipe)
  dRM.ipe=unstandardizedES.ipe/sqrt(var.within.ipe)
  r.ipe=var.between.ipe/var.sig.ipe
  N1=N2=length(simulationData$pid)/4
  df=N1+N2-2
  c=1-3/(4*df-1)
  var.t.ipe=(df/(df-2))*(1+t.ipe^2)- t.ipe^2/c^2
  var.dRM.Fromt.ipe = var.t.ipe*((N1+N2)/(2*N1*N2)) #var.t/N1
  var.dRM.ipe=(df/(df-2))*((N1+N2)/(2*N1*N2)+dRM.ipe^2)- dRM.ipe^2/c^2
  dRM.Fromt.ipe=t.ipe*sqrt((N1+N2)/(2*N1*N2))
  var.dRM.Fromt2.ipe=(df/(df-2))*((N1+N2)/(2*N1*N2)+dRM.Fromt.ipe^2)- dRM.Fromt.ipe^2/c^2
  var.dIG.ipe=(df/(df-2))*((1-r.ipe)*(N1+N2)/(2*N1*N2)+dIG.ipe^2)- dIG.ipe^2/c^2
  dIG.Fromt.ipe=t.ipe*sqrt(1-r.ipe)*sqrt((N1+N2)/(2*N1*N2))
  var.dIG.Fromt.ipe=var.t.ipe*(1-r.ipe)*((N1+N2)/(2*N1*N2))

  #d<-cbind(dIG,dRM)
  #return(d)
  effectSizes.ipe<-data.frame(dIG.ipe, var.dIG.ipe, dRM.ipe, var.dRM.ipe, dIG.Fromt.ipe, var.dIG.Fromt.ipe, dRM.Fromt.ipe, var.dRM.Fromt.ipe, var.dRM.Fromt2.ipe, unstandardizedES.ipe, var.sig.ipe, var.within.ipe, var.between.ipe, t.ipe, var.t.ipe)
  return(effectSizes.ipe)
}


#' @title reproduceSimulationResultsBasedOn500Reps1000Obs
#' @description Function to calculate simulation results based on 500 repetitions of 1000 observation samples. Function is used in a paper "Effect Sizes and their Variance for AB/BA Crossover Design Studies" by Lech Madeyski and Barbara Kitchenham.
#' @author Lech Madeyski and Barbara Kitchenham
#' @export reproduceSimulationResultsBasedOn500Reps1000Obs
#' @return data frame including the following simulation results:
#' # treatmentEffect.Ave - Average Technique Effect
#' # dRM.Ave - Average dRM
#' # dRM.Var - Variance of dRM
#' # dRM.Var.Ave - Average of var(dRM)
#' # dRM.Var.ModerateSampleSizeApprox -
#' # dIG.Ave - Average dIG
#' # dIG.Var - Variance of dIG
#' # dIG.Var.Ave - Average of var(dIG)
#' # dIG.Var.ModerateSampleSizeApprox -
#' # dIG.Var.CurtinAve -
#' @examples
#' # return simulation results based on 500 repetitions of 1000 observation samples
#' simulationResultsTable500x1000<-reproduceSimulationResultsBasedOn500Reps1000Obs()
reproduceSimulationResultsBasedOn500Reps1000Obs <- function()
{
  set.seed(42)

  #Covariance matrix
  #r=.75
  var=25
  covar=18.75

  sigmaA=matrix(c(var,covar,covar,var),nrow=2,ncol=2)

  #Mean for treatment sequence A
  #common mean mu=50, treatment difference = 10, period effect=5
  meanA=c(50,65)
  #Mean for sequence B
  meanB=c(60,55)

  #Simulation aimed at validating the estimates of the effect size variance
  reps=500
  treat=c(1:reps)
  period=c(1:reps)
  diffvar=c(1:reps)
  effectsizevar=c(1:reps)
  effectsize1=c(1:reps) #dRM
  withinvar=c(1:reps)
  varwithinest=c(1:reps)
  effectsize2=c(1:reps) #dIG estimated from correct variance
  effectsize2a=c(1:reps) #dIG estimated from relationship with dRM
  effectsize2b=c(1:reps) #dIG estimated assuming r is known
  r=c(1:reps)
  lsvardrm=c(1:reps)
  lsvardrmrev=c(1:reps)
  tRM=c(1:reps) # Correct formulation of t-test for repeated measures
  vartrm=c(1:reps)
  vardrm=c(1:reps)
  vardig=c(1:reps)
  curtinvar=c(1:reps)
  lsvardig=c(1:reps) #Large sample approximation variance
  lsvardigrev=c(1:reps)
  Num=reps
  N1=N2=Num
  df=N1+N2-2
  c=1-(3/(4*df-1))


  for (i in 1:reps)
  {
    mydata1=MASS::mvrnorm(Num, meanA, sigmaA)

    #Sample data for Sequence B

    mydata2=MASS::mvrnorm(Num, meanB, sigmaA)

    mydata1=as.data.frame(mydata1)
    mydata2=as.data.frame(mydata2)
    names(mydata1)=c("AP1","AP2")
    names(mydata2)=c("BP1","BP2")

    diff1=mydata1$AP2-mydata1$AP1
    diff2=mydata2$BP2-mydata2$BP1

    treat[i]=(mean(diff1)-mean(diff2))/2
    period[i]=(mean(diff1)+mean(diff2))/2

    diffvar[i]=(stats::var(diff1)+stats::var(diff2))/2
    effectsizevar[i]=diffvar[i]/2

    effectsize1[i]=treat[i]/(sqrt(effectsizevar[i]))

    withinvar[i]=(stats::var(mydata1$AP1)+stats::var(mydata1$AP2)+stats::var(mydata2$BP1)+stats::var(mydata2$BP2))/4

    effectsize2[i]=treat[i]/sqrt(withinvar[i])

    r[i]=(withinvar[i]-(diffvar[i]/2))/withinvar[i]
    effectsize2a[i]=effectsize1[i]*sqrt(1-r[i])
    effectsize2b[i]=effectsize1[i]*sqrt(0.25)

    varwithinest[i]=diffvar[i]/(2*(1-r[i]))

    tRM[i]=effectsize1[i]*sqrt((2*N1*N2)/(N1+N2))
    vartrm[i]= (df/(df-2))*((1+tRM[i]^2))- tRM[i]^2/c^2
    vardrm[i]=vartrm[i]*(N1+N2)/(2*N1*N2)
    vardig[i]=vardrm[i]*(1-r[i])
    lsvardrm[i]=(N1+N2)/(2*N1*N2)+effectsize1[i]^2/(2*(N1+N2-2))

    lsvardig[i]=(1-r[i])*(N1+N2)/(2*N1*N2)+effectsize2[i]^2/(2*(N1+N2-2))


    curtinvar[i]=(df/(df-2))*(((N1+N2)/(4*N1*N2)+effectsize2[i]^2))- effectsize2[i]^2/c^2
  }
  treatmentEffect.Ave <- mean(treat) #Overall average of technique effect
  stats::var(treat)
  mean(period) #Overall average of period effect
  stats::var(period)
  dRM.Ave <- mean(effectsize1) #Overall average of dRM
  dRM.Var <- stats::var(effectsize1)   #Empirical estimate of variance of dRM
  mean(effectsize2) #Overall average of dIG based on within-subject within group variance
  stats::var(effectsize2) #Empirical estimate of variance of dRM
  dIG.Ave <- mean(effectsize2a) #Overall average of dIG based on relationship with dRM
  dIG.Var <- stats::var(effectsize2a) #Empirical estimate of variance of dIG
  mean(effectsize2b)
  stats::var(effectsize2b) #Empirical estimate of variance of dRM assuming r=0.75 known

  mean(r) #Overall average estimate of between subject correlation

  mean(tRM) #Overall estimate of t-value
  mean(vartrm) #Average of 500 estimates of variance of t-value
  dRM.Var.Ave <- mean(vardrm) #Average of 500 estimates of variance of dRM
  dIG.Var.Ave <- mean(vardig) #Average of 500 estimates of variance of dIG
  mean(withinvar) # Overall average within-subject within group variance estimates
  dIG.Var.CurtinAve <- mean(curtinvar) #Overall average of Curtin's estimates of variance of dIG
  dRM.Var.ModerateSampleSizeApprox <- mean(lsvardrm) # Average of 500 moderate sample estimates of var of dRM
  dIG.Var.ModerateSampleSizeApprox <- mean(lsvardig) # Average of 500 moderate sample estimates of var of dIG
  stats::var(effectsize1)*(1-mean(r)) #Estimate of varDIG based on relationsgip with varDRM

  simultationResults<-data.frame(treatmentEffect.Ave, dRM.Ave, dRM.Var, dRM.Var.Ave, dRM.Var.ModerateSampleSizeApprox, dIG.Ave, dIG.Var, dIG.Var.Ave, dIG.Var.ModerateSampleSizeApprox, dIG.Var.CurtinAve)
  return(simultationResults)
}




#' @title percentageInaccuracyOfLargeSampleVarianceApproximation
#' @description Plot the extent of inaccuracy using the large sample approximate effect size variance on 4 related graphs corresponding to the four different correlation values. Plot visualizes the relationship between sample size and effect size and the percentage inaccuracy of the large sample variance approximation. Function is used in a paper "Effect Sizes and their Variance for AB/BA Crossover Design Studies" by Lech Madeyski and Barbara Kitchenham.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export percentageInaccuracyOfLargeSampleVarianceApproximation
#' @param data - data behind the plot returned by getSimulatedCrossoverDataSets() or stored in reproducer::KitchenhamMadeyski.SimulatedCrossoverDataSets
#' @return plot described in description
#' @examples
#' data <- KitchenhamMadeyski.SimulatedCrossoverDataSets
#' myPlot <- percentageInaccuracyOfLargeSampleVarianceApproximation(data)
percentageInaccuracyOfLargeSampleVarianceApproximation <- function(data)
{
  #pdf("AccuracyPlotsAll.pdf")
  graphics::par(mfrow=c(2,2))

  d <- data[c(1:7),]
  graphics::plot(d$actualSampleSize,d$Accuracy,ylab="Approx Variance Inaccuracy with r=0", xlab="Sample Size",type="b",pch=0, lty=1,ylim=c(0,30),cex.sub=.8, main="(a)")

  d <- data[c(8:14),]
  graphics::lines(d$actualSampleSize, d$Accuracy,type="b",lty=2, pch=1)

  d <- data[c(15:21),]
  graphics::lines(d$actualSampleSize, d$Accuracy,type="b",lty=3, pch=2)

  graphics::legend("topright",inset=.01,title="Effect Size", c("0","2.5","5"),lty=c(1,2,3),pch=c(0,1,2),cex=.7)

  d <- data[c(22:28),]
  graphics::plot(d$actualSampleSize, d$Accuracy,ylab="Approx Variance Inaccuracy with r=0.25", xlab="Sample Size",type="b",pch=0, lty=1,ylim=c(0,40),cex.sub=.8,main="(b)")

  d <- data[c(29:35),]
  graphics::lines(d$actualSampleSize, d$Accuracy,type="b",lty=2, pch=1)

  d <- data[c(36:42),]
  graphics::lines(d$actualSampleSize, d$Accuracy,type="b",lty=3, pch=2)

  graphics::legend("topright",inset=.01,title="Effect Size", c("0","2.5","5"),lty=c(1,2,3),pch=c(0,1,2),cex=.7)

  d <- data[c(43:49),]
  graphics::plot(d$actualSampleSize, d$Accuracy,ylab="Approx Variance Inaccuracy with r=0.5", xlab="Sample Size",type="b",pch=0, lty=1,ylim=c(0,40),cex.sub=.8, main="(c)")

  d <- data[c(50:56),]
  graphics::lines(d$actualSampleSize, d$Accuracy,type="b",lty=2, pch=1)

  d <- data[c(57:63),]
  graphics::lines(d$actualSampleSize, d$Accuracy,type="b",lty=3, pch=2)

  graphics::legend("topright",inset=.01,title="Effect Size", c("0","2.5","5"),lty=c(1,2,3),pch=c(0,1,2),cex=.7)

  d <- data[c(64:70),]
  graphics::plot(d$actualSampleSize, d$Accuracy,ylab="Approx Variance Inaccuracy with r=0.75", xlab="Sample Size",type="b",pch=0, lty=1,ylim=c(0,40),cex.sub=.8 , main="(d)")

  d <- data[c(71:77),]
  graphics::lines(d$actualSampleSize, d$Accuracy,type="b",lty=2, pch=1)

  d <- data[c(78:84),]
  graphics::lines(d$actualSampleSize, d$Accuracy,type="b",lty=3, pch=2)

  graphics::legend("topright",inset=.01,title="Effect Size", c("0","2.5","5"),lty=c(1,2,3),pch=c(0,1,2),cex=.7)
  #dev.off()
}

#' @title proportionOfSignificantTValuesUsingCorrectAnalysis
#' @description Plots visualize the relationship between sample size, effect size and the proportion of significant t-values using the correct analysis. Function is used in a paper "Effect Sizes and their Variance for AB/BA Crossover Design Studies" by Lech Madeyski and Barbara Kitchenham.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export proportionOfSignificantTValuesUsingCorrectAnalysis
#' @param data - data behind the plot returned by getSimulatedCrossoverDataSets() or stored in reproducer::KitchenhamMadeyski.SimulatedCrossoverDataSets
#' @return plot described in description
#' @examples
#' data <- KitchenhamMadeyski.SimulatedCrossoverDataSets
#' myPlot <- proportionOfSignificantTValuesUsingCorrectAnalysis(data)
proportionOfSignificantTValuesUsingCorrectAnalysis <- function(data)
{
  #pdf("SignificantValuesPlots.pdf")
  graphics::par(mfrow=c(2,2))

  d <- data[c(1:7),]
  graphics::plot(d$actualSampleSize, d$PropSig, ylab="Prop Significant t-Values with r=0", xlab="Sample Size",type="b",pch=0, lty=1,ylim=c(0,100),xlim=c(0,100),cex.sub=.8, main="(a)")

  d <- data[c(8:14),]
  graphics::lines(d$actualSampleSize, d$PropSig,type="b",lty=2, pch=1)

  d <- data[c(15:21),]
  graphics::lines(d$actualSampleSize, d$PropSig,type="b",lty=3, pch=2)

  graphics::legend("right",inset=.01,title="Effect Size", c("0","2.5","5"),lty=c(1,2,3),pch=c(0,1,2),cex=.7)

  d <- data[c(22:28),]
  graphics::plot(d$actualSampleSize, d$PropSig, ylab="Prop Significant t-Values with r=0.25", xlab="Sample Size",type="b",pch=0, lty=1,ylim=c(0,100),xlim=c(0,100),cex.sub=.8, main="(b)")

  d <- data[c(29:35),]
  graphics::lines(d$actualSampleSize, d$PropSig,type="b",lty=2, pch=1)

  d <- data[c(36:42),]
  graphics::lines(d$actualSampleSize, d$PropSig,type="b",lty=3, pch=2)

  graphics::legend("right",inset=.01,title="Effect Size", c("0","2.5","5"),lty=c(1,2,3),pch=c(0,1,2),cex=.7)

  d <- data[c(43:49),]
  graphics::plot(d$actualSampleSize, d$PropSig, ylab="Prop Significant t-Values with r=0.5", xlab="Sample Size",type="b",pch=0, lty=1,ylim=c(0,100),xlim=c(0,100),cex.sub=.8, main="(c)")

  d <- data[c(50:56),]
  graphics::lines(d$actualSampleSize, d$PropSig, type="b",lty=2, pch=1)

  d <- data[c(57:63),]
  graphics::lines(d$actualSampleSize, d$PropSig, type="b",lty=3, pch=2)

  graphics::legend("right",inset=.01,title="Effect Size", c("0","2.5","5"),lty=c(1,2,3),pch=c(0,1,2),cex=.7)

  d <- data[c(64:70),]
  graphics::plot(d$actualSampleSize, d$PropSig, ylab="Prop Significant t-Values with r=0.75", xlab="Sample Size",type="b",pch=0, lty=1,ylim=c(0,100),xlim=c(0,100),cex.sub=.8, main="(d)")

  d <- data[c(71:77),]
  graphics::lines(d$actualSampleSize, d$PropSig, type="b",lty=2, pch=1)

  d <- data[c(78:84),]
  graphics::lines(d$actualSampleSize, d$PropSig, type="b",lty=3, pch=2)

  graphics::legend("right",inset=.01,title="Effect Size", c("0","2.5","5"),lty=c(1,2,3),pch=c(0,1,2),cex=.7)

  #dev.off()
}


#' @title proportionOfSignificantTValuesUsingIncorrectAnalysis
#' @description Plots visualize the relationship between sample size, effect size and the proportion of significant t-values using the incorrect analysis. Function is used in a paper "Effect Sizes and their Variance for AB/BA Crossover Design Studies" by Lech Madeyski and Barbara Kitchenham.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export proportionOfSignificantTValuesUsingIncorrectAnalysis
#' @param data - data behind the plot returned by getSimulatedCrossoverDataSets() or stored in reproducer::KitchenhamMadeyski.SimulatedCrossoverDataSets
#' @return plot described in description
#' @examples
#' data  <- KitchenhamMadeyski.SimulatedCrossoverDataSets
#' myPlot <- proportionOfSignificantTValuesUsingIncorrectAnalysis(data)
proportionOfSignificantTValuesUsingIncorrectAnalysis <- function(data)
{
  #pdf("SignificantValuesPlots.pdf")
  graphics::par(mfrow=c(2,2))

  d <- data[c(1:7),]
  graphics::plot(d$actualSampleSize, d$WrongTSig, ylab="Prop Significant Incorrect t-Values with r=0", xlab="Sample Size",type="b",pch=0, lty=1,ylim=c(0,100),xlim=c(0,100),cex.sub=.8, main="(a)")

  d <- data[c(8:14),]
  graphics::lines(d$actualSampleSize, d$WrongTSig,type="b",lty=2, pch=1)

  d <- data[c(15:21),]
  graphics::lines(d$actualSampleSize, d$WrongTSig, type="b",lty=3, pch=2)

  graphics::legend("right",inset=.01,title="Effect Size", c("0","2.5","5"),lty=c(1,2,3),pch=c(0,1,2),cex=.7)

  d <- data[c(22:28),]
  graphics::plot(d$actualSampleSize, d$WrongTSig, ylab="Prop Significant Incorrect t-Values with r=0.25", xlab="Sample Size",type="b",pch=0, lty=1,ylim=c(0,100),xlim=c(0,100),cex.sub=.8, main="(b)")

  d <- data[c(29:35),]
  graphics::lines(d$actualSampleSize, d$WrongTSig,type="b",lty=2, pch=1)

  d <- data[c(36:42),]
  graphics::lines(d$actualSampleSize, d$WrongTSig,type="b",lty=3, pch=2)

  graphics::legend("right",inset=.01,title="Effect Size", c("0","2.5","5"),lty=c(1,2,3),pch=c(0,1,2),cex=.7)

  d <- data[c(43:49),]
  graphics::plot(d$actualSampleSize, d$WrongTSig, ylab="Prop Significant Incorrect t-Values with r=0.5", xlab="Sample Size",type="b",pch=0, lty=1,ylim=c(0,100),xlim=c(0,100),cex.sub=.8, main="(c)")

  d <- data[c(50:56),]
  graphics::lines(d$actualSampleSize, d$WrongTSig, type="b",lty=2, pch=1)

  d <- data[c(57:63),]
  graphics::lines(d$actualSampleSize, d$WrongTSig, type="b",lty=3, pch=2)

  graphics::legend("right",inset=.01,title="Effect Size", c("0","2.5","5"),lty=c(1,2,3),pch=c(0,1,2),cex=.7)

  d <- data[c(64:70),]
  graphics::plot(d$actualSampleSize, d$WrongTSig, ylab="Prop Significant Incorrect t-Values with r=0.75", xlab="Sample Size",type="b",pch=0, lty=1,ylim=c(0,100),xlim=c(0,100),cex.sub=.8, main="(d)")

  d <- data[c(71:77),]
  graphics::lines(d$actualSampleSize, d$WrongTSig, type="b",lty=2, pch=1)

  d <- data[c(78:84),]
  graphics::lines(d$actualSampleSize, d$WrongTSig, type="b",lty=3, pch=2)

  graphics::legend("right",inset=.01,title="Effect Size", c("0","2.5","5"),lty=c(1,2,3),pch=c(0,1,2),cex=.7)

}
