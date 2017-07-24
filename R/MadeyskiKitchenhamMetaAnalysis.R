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
#' Functions related to a paper "Effect sizes and their variance for AB/BA crossover design studies" by Lech Madeyski and Barbara Kitchenham
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
#' data<-getSimulationData(25, 18.75, 50, 10, 5, 15)
getSimulationData <- function(var, covar, meanA1, treatmentDiff, periodEffect, numOfSamples)
{
  #Covariance matrix
  ##r=.75
  #var=25
  #covar=18.75
  sigmaA=matrix(c(var,covar,covar,var),nrow=2,ncol=2)

  #Mean for treatment sequence B
  #common mean mu=50, treatment difference = 10, period effect=5
  #meanB=c(50,65) SG2
  meanB=c(meanA1,meanA1+treatmentDiff+periodEffect)
  #Mean for sequence A
  #meanA=c(60,55) SG1
  meanA=c(meanA1+treatmentDiff,meanA1+periodEffect)


  #Sample data for sequence SG2
  set.seed(123)


  Num<-numOfSamples

  SG2data=MASS::mvrnorm(Num, meanB, sigmaA)

  #Sample data for Sequence SG1

  set.seed(9123)
  SG1data=MASS::mvrnorm(Num, meanA, sigmaA)

  SG1data=as.data.frame(SG1data)
  SG2data=as.data.frame(SG2data)
  names(SG1data)=c("T1P1","T2P2")
  names(SG2data)=c("T2P1","T1P2")


  p1=rep("P1",Num)
  p2=rep("P2",Num)
  t1=rep("T1",Num)
  t2=rep("T2",Num)
  seq1=seq(1,Num)
  seq2=seq(Num+1,2*Num)
  s1=rep("S1",Num)
  s2=rep("S2", Num)

  SG2P1T2=data.frame(pid=seq2,period=p1,sequence=s2, result=SG2data$T2P1,technique=t2)
  SG1P1T1=data.frame(pid=seq1,period=p1,sequence=s1, result=SG1data$T1P1,technique=t1)

  SG2P2T1=data.frame(pid=seq2, period=p2,sequence=s2, result=SG2data$T1P2,technique=t1)
  SG1P2T2=data.frame(pid=seq1, period=p2,sequence=s1, result=SG1data$T2P2,technique=t2)

  allT1=rbind(SG1P1T1, SG2P2T1)
  allT2=rbind(SG1P2T2, SG2P1T2)

  SimulationData=rbind(allT1, allT2)

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

  data<-getSimulationData(25, 18.75, 50, 10, 5, 15)


  newdataT1=data[which(data$technique=="T1"),]
  newdataT2=data[which(data$technique=="T2"),]

  SG1P1T1=newdataT1[which(newdataT1$period=="P1"),]
  SG2P2T1=newdataT1[which(newdataT1$period=="P2"),]
  SG1P2T2=newdataT2[which(newdataT2$period=="P2"),]
  SG2P1T2=newdataT2[which(newdataT2$period=="P1"),]

  DiffSG1=SG1P1T1$result-SG1P2T2$result
  SumSG1=SG1P1T1$result+SG1P2T2$result
  newSG1=data.frame(pid=SG1P1T1$pid,Sequence=SG1P1T1$sequence,Difference=DiffSG1,Sum=SumSG1)

  DiffSG2=SG2P2T1$result-SG2P1T2$result
  SumSG2=SG2P2T1$result+SG2P1T2$result
  newSG2=data.frame(pid=SG2P2T1$pid,Sequence=SG2P2T1$sequence,Difference=DiffSG2,Sum=SumSG2)

  newAll=rbind(newSG1,newSG2)


  graphics::par(mfrow=c(1,2),cex=.75)
  graphics::boxplot(Difference ~ Sequence,data=newAll,ylab="Subject Difference for each Sequence", xlab="Sequence ID",main="(b) Crossover Differences for each Sequence",cex=0.8)

  graphics::plot(newdataT1$pid,newdataT1$result,ylim=c(30,80),xlim=c(1,30),ylab="Outcome",xlab="Subject ID",pch=0,main="(a) Raw Data for Each Individual",typ="b")
  graphics::points(newdataT2$pid,newdataT2$result,pch=3,typ="b")
  graphics::abline(v=15.5)
  graphics::text(0,75,"Sequence S1", pos=4)
  graphics::text(16,75, "Sequence S2", pos=4)
  graphics::legend("bottomright",title="Technique", inset=0.05, c("1","2"),pch=c(0,3))

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
#' #OR
#' simulationData<-getSimulationData(25, 18.75,50,10,5,15)
#' es<-getEffectSizesABBA(simulationData) #return effect sizes and variances
getEffectSizesABBA <- function(simulationData)
{
  #This function calculates the correct variances
  # for all AB/BA standardized effect sizes for
  # small and medium size samples.
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
  N1=N2=length(simulationData$pid)/4 #this statement introduces the problem with the generality of the function getEffectSizesABBA  as it assumes equal numbers of values in each sequence group

  df=N1+N2-2
  c=1-3/(4*df-1)

  t.unbiased=c*t
  gIG=dIG*c
  gRM=dRM*c

  var.t=(df/(df-2))*(1+t^2)- t^2/c^2
  var.t.unbiased=(df/(df-2))*(1+t.unbiased^2)- t.unbiased^2/c^2
  var.dRM.Fromt = var.t*((N1+N2)/(2*N1*N2)) #var.t/N1
  var.dRM.Fromt.unbiased=var.t.unbiased*((N1+N2)/(2*N1*N2))
  var.dRM.biased=(df/(df-2))*((N1+N2)/(2*N1*N2)+dRM^2)- dRM^2/c^2 #Equation 42
  var.dRM.unbiased = (df/(df-2))*((N1+N2)/(2*N1*N2)+gRM^2)- gRM^2/c^2

  dRM.Fromt=t*sqrt((N1+N2)/(2*N1*N2))
  dRM.Fromt.unbiased=t.unbiased*sqrt((N1+N2)/(2*N1*N2))

  var.gRM.Fromt=var.t*((N1+N2)/(2*N1*N2))*c^2
  var.gRM.Fromt.unbiased=var.t.unbiased*((N1+N2)/(2*N1*N2))*c^2
  var.gRM.From.vardRM.unbiased= c^2*var.dRM.unbiased
  var.gRM  = c^2*(df/(df-2))*((N1+N2)/(2*N1*N2)+gRM^2)- gRM^2

  dIG.Fromt=t*sqrt(1-r)*sqrt((N1+N2)/(2*N1*N2))
  dIG.Fromt.unbiased=t.unbiased*sqrt(1-r)*sqrt((N1+N2)/(2*N1*N2))
  var.dIG.Fromt= var.t*(1-r)*((N1+N2)/(2*N1*N2))
  var.dIG.Fromt.unbiased= var.t.unbiased*(1-r)*((N1+N2)/(2*N1*N2))
  var.dIG.biased= (df/(df-2))*((1-r)*(N1+N2)/(2*N1*N2)+dIG^2)- dIG^2/c^2
  var.dIG.unbiased= (df/(df-2))*((1-r)*(N1+N2)/(2*N1*N2)+gIG^2)- gIG^2/c^2

  var.gIG.Fromt=var.t*(1-r)*((N1+N2)/(2*N1*N2))*c^2
  var.gIG.From.vardIG.unbiased=c^2*var.dIG.unbiased
  var.gIG.Fromt.unbiased=var.t.unbiased*(1-r)*((N1+N2)/(2*N1*N2))*c^2
  var.gIG=c^2*(df/(df-2))*((1-r)*(N1+N2)/(2*N1*N2)+gIG^2)- gIG^2




  var.t.Approx=1+t^2/(2*(N1+N2-2))
  var.t.Approx.unbiased=1+t.unbiased^2/(2*(N1+N2-2))

  var.gRM.Approx=((N1+N2)/(2*N1*N2)) + (gRM^2)/(2*(N1+N2-2))#see paper and Equation 49

  var.gIG.Approx=(((N1+N2)*(1-r))/(2*N1*N2)) + (gIG^2)/(2*(N1+N2-2))#see paper and Equation 50

  var.dRM.Approx.alternative=(N1+N2)/(2*N1*N2) + dRM^2/(2*(N1+N2-3.94))
  var.dIG.Approx.alternative=(1-r)*(N1+N2)/(2*N1*N2) + dIG^2/(2*(N1+N2-3.94))

  PRE.vardrm=100*(var.dRM.unbiased-var.dRM.Approx.alternative)/var.dRM.unbiased
  PRE.vardig=100*(var.dIG.unbiased-var.dIG.Approx.alternative)/var.dIG.unbiased
  PRE.vargrm=100*(var.gRM-var.gRM.Approx)/var.gRM
  PRE.vargig=100*(var.gIG-var.gIG.Approx)/var.gIG
  #effectSizes<-dgta.frame(unstandardizedES, r, t,t.unbiased,var.t, var.t.unbiased,dRM, dRM.Fromt, dRM.Fromt.unbiased, gRM,var.dRM.Fromt,var.dRM.Fromt.unbiased,var.dRM.biased, var.dRM.unbiased, var.gRM, var.gRM.Fromt, var.gRM.Fromt.unbiased,var.gRM.From.vardRM.unbiased, dIG, dIG.Fromt, dIG.Fromt.unbiased, var.dIG.Fromt, var.dIG.Fromt.unbiased,var.dIG.biased,var.dIG.unbiased, gIG, var.gIG.Fromt, var.gIG.Fromt.unbiased, var.gIG.From.vardIG.unbiased, var.gIG, var.t.Approx, var.t.Approx.unbiased, var.gRM.Approx,var.gIG.Approx, var.dRM.Approx.alternative, var.dIG.Approx.alternative,PRE.vardrm, PRE.vardig, PRE.vargrm,PRE.vargig,var.sig, var.within, var.between)

  effectSizes<-data.frame(unstandardizedES, r, t,t.unbiased, var.t.unbiased,dRM, var.dRM.unbiased,gRM, var.gRM, dIG, var.dIG.unbiased, gIG, var.gIG, var.t.Approx, var.t.Approx.unbiased, var.gRM.Approx,var.gIG.Approx, var.dRM.Approx.alternative, var.dIG.Approx.alternative,PRE.vardrm, PRE.vardig, PRE.vargrm,PRE.vargig,var.sig, var.within, var.between)
  return(effectSizes)
}


#' @title getTheoreticalEffectSizeVariancesABBA
#' @description Function provides the theoretical value of the t-statistic, variance of t, and variance of the effect sizes based on the parameters built into crossover model data simulated by the getSimilationData() function.
#' Function is used in a paper "Effect Sizes and their Variance for AB/BA Crossover Design Studies" by Lech Madeyski and Barbara Kitchenham.
#' @author Lech Madeyski and Barbara Kitchenham
#' @export getTheoreticalEffectSizeVariancesABBA
#' @param theoreticalvarW - The within subject variance used to construct the simulation, i.e., the built-in Variance - the built-in Covariance
#' @param theoreticalTechniqueEffect - The technique effect built into the crossover model data
#' @param theoreticalrho - The between subject correlation built into the crossover model simulation data
#' @param N1 - The number of subjects in sequence group 1 in the crossover model simulation
#' @param N2 - The number of subjects in sequence group 2 in the crossover model simulation
#' @return data frame incl. calculated:
#' theoreticalt - the theoretical value of the t-statistic
#' theoreticalvart - variance of t
#' theoreticalvardIG -  variance of the effect size dIG based on the parameters built into crossover model data simulated by the getSimilationData function
#' theoreticalvardRM -  variance of the effect size dRM based on the parameters built into crossover model data simulated by the getSimilationData function
#' @examples
#' # Generates data used in Table 15 of the paper
#' theoreticalEffectSizeVariances <- getTheoreticalEffectSizeVariancesABBA(6.25,-10,0.75,15,15)
getTheoreticalEffectSizeVariancesABBA <- function(theoreticalvarW,theoreticalTechniqueEffect,theoreticalrho,N1,N2)
{
  df=N1+N2-2
  c=1-3/(4*df-1)

  theoreticaldRM=theoreticalTechniqueEffect/sqrt(theoreticalvarW)
  theoreticalt=theoreticaldRM/sqrt((N1+N2)/(2*N2*N1))
  theoreticalvart=(df/(df-2))*(1+theoreticalt^2)-theoreticalt^2/c^2
  theoreticalvardRM=theoreticalvart*((N1+N2)/(2*N1*N2))
  theoreticalvardIG=theoreticalvardRM*(1-theoreticalrho)

  theoreticaleffectSizeVariancess<-data.frame(theoreticalt, theoreticalvart, theoreticalvardIG, theoreticalvardRM)

  return(theoreticaleffectSizeVariancess)
}



#' @title getEffectSizesABBAIgnoringPeriodEffect
#' @description Function to calculate both effect sizes (dIG.ipe, dRM.ipe), i.e., independent groups and repeated measures standardized effect sizes and variances, for AB/BA crossover design studies ignoring period effect (thus wrong). Function was removed in the revision of the paper "Effect Sizes and their Variance for AB/BA Crossover Design Studies" by Lech Madeyski and Barbara Kitchenham.
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
  #curtinvar=c(1:reps) # removed since 0.1.7
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


    # curtinvar[i]=(df/(df-2))*(((N1+N2)/(4*N1*N2)+effectsize2[i]^2))- effectsize2[i]^2/c^2  # removed since 0.1.7
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
  #dIG.Var.CurtinAve <- mean(curtinvar) #Overall average of Curtin's estimates of variance of dIG (removed since 0.1.7)
  dRM.Var.ModerateSampleSizeApprox <- mean(lsvardrm) # Average of 500 moderate sample estimates of var of dRM
  dIG.Var.ModerateSampleSizeApprox <- mean(lsvardig) # Average of 500 moderate sample estimates of var of dIG
  stats::var(effectsize1)*(1-mean(r)) #Estimate of varDIG based on relationsgip with varDRM

  simultationResults<-data.frame(treatmentEffect.Ave, dRM.Ave, dRM.Var, dRM.Var.Ave, dRM.Var.ModerateSampleSizeApprox, dIG.Ave, dIG.Var, dIG.Var.Ave, dIG.Var.ModerateSampleSizeApprox) # removed ", dIG.Var.CurtinAve)" since 0.1.7

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



#' @title effectSizeCI
#' @description 95% Confidence Intervals (CI) on Standardised Effect Sizes (d) for cross-over repeated-measures, before-after repeated-measures, and independent group experimental designs
#' The procedure is based on finding the upper and lower 0.025 bounds for the related t-variable.
#' The t-variable needds to be adjusted for bias by multiplying by c
#' The upper and lower bounds on the t-variable are then used to calculate to upper and lower bounds on the
#' repeated measures effect size (d_RM) by multiplying the upper and lower bound of the t-variable by sqrt((n1+n2)/(2*(n1*n2))).
#' Upper and lower bounds on the equivalent independent groups effect size (d_IG) are found by multiplying the upper and lower
#' bounds on d_RM by sqrt(1-r).
#' @author Lech Madeyski and Barbara Kitchenham
#' @export effectSizeCI
#' @param expDesign Experimental design: 1) crossover repeated measures ("CrossOverRM"), 2) before-after repeated measures (expDesign=="BeforeAfterRM"), 3) independent groups ("IG)
#' @param t t-statistics (t must be less than or equal to 37.62, the limit from the R function documentation)
#' @param n1 The number of observations in sequence group 1 (expDesign=="CrossOverRM"), the number of observations in group 1 (expDesign=="IG"), or the total number of observations (expDesign=="BeforeAfterRM")
#' @param n2 The number of observations in sequence group 2 (expDesign=="CrossOverRM") or the number of observations in group 2 (expDesign=="IG")
#' @param r The correlation between outcomes for individual subject (the within subject correlation)
#' @param epsilon The precision of the iterative procedure
#' @param maxsteps The maximum number of steps of the iterative procedure (the procedure terminates at maxsteps or earlier if CI with enough precision have been calculated)
#' @param stepsize The size of steps (influences the convergence of the calculations, i.e., the number of steps required to obtain the final result of precison defined by the epsilon)
#' @return A list of Confidence Intervals for: t-statistic (t_LB and t_UB), repeted-measures effect size d_RM (d_RM_LB, d_RM_UB), independent groups effect size (d_IG_LB, d_IG_UB)
#' @examples
#' effectSizeCI(expDesign="CrossOverRM", t=14.4, n1=15, n2=15, r=0.6401)
#' effectSizeCI(expDesign = "BeforeAfterRM", t=14.16536, n1=15, n2=0, r=0.6146771)
#' effectSizeCI(expDesign = "IG", t=-6.344175, n1=15, n2=15)
#' effectSizeCI(expDesign="CrossOverRM", t=0.5581, n1=6, n2=6, r=0.36135)
#' effectSizeCI(expDesign = "CrossOverRM", r=0.855,t=4.33, n1=7, n2=6)
effectSizeCI <- function(expDesign, t, n1, n2, r=0, epsilon=0.0000000001, maxsteps=1000, stepsize=3)
{
  #stepsize=3 #influences the number of steps required to obtain the final result (i.e., the result of precison defined by the epsilon), e.g., consider stepsize=2
  if(t>37.62)
    t <- 37.62
  switch(expDesign,
         CrossOverRM={
           df <- n1+n2-2
           d_RM <- t*sqrt((n1+n2)/(2*n1*n2))
         },
         IG={
           df <- n1+n2-2
           d_IG <- t*sqrt((1/n1)+(1/n2))
         },
         BeforeAfterRM={
           df <- n1-1
           d_RM <- t*sqrt(1/n1)
         },
         stop("Need to specify the expDesign parameter\n.")
  )

  c <- 1-3/(4*df-1)
  t.unbiased=t*c

  vart <- (df/(df-2))*(1+t.unbiased^2)-t.unbiased^2/c^2

  i <- 1
  t_LB <- t.unbiased-2*sqrt(vart)
  pt_LB <- stats::pt(t_LB, df=df, ncp=t.unbiased)
  t_UB <- t.unbiased+2*sqrt(vart)
  pt_UB <- stats::pt(t_UB, df=df, ncp=t.unbiased)
  while( i <= maxsteps & ( abs(pt_LB) >= abs(0.025+epsilon) |  abs(pt_LB) <= abs(0.025-epsilon) | abs(pt_UB) >= abs(0.975+epsilon) | abs(pt_UB) <= abs(0.975-epsilon) ) )
  {
    if(i==1)
    {
      # I'm not sure what this condition does.
      t_LB <- (t.unbiased+t_LB)/2
      t_UB <- (t.unbiased+t_UB)/2
    }
    else
    {
      if(abs(pt_LB) <= abs(0.025-epsilon))       {
        if(t_LB>0) # we need larger value of t_LB
        {t_LB <- (1+(stepsize*abs(pt_LB-0.025))) * t_LB}
        else  # we need smaller value of t_LB
        {
          t_LB <- (1-(stepsize*abs(pt_LB-0.025))) * t_LB
        }


      }
      if(abs(pt_LB) >= abs(0.025+epsilon))
      {
        if(t_LB>0) # we need smaller value of t_LB

        {t_LB <- (1-(stepsize*abs(pt_LB-0.025))) * t_LB}
        else # we need a larger value of T_LB
        {
          t_LB <- (1+(stepsize*abs(pt_LB-0.025))) * t_LB
        }

      }

      if(abs(pt_UB) <= abs(0.975-epsilon))
      {
        if(t_UB>0) # we need larger value of t_UB

        {t_UB <- (1+(stepsize*abs(pt_UB-0.975))) * t_UB}
        else # we need a smaller value of t_UB
        {
          t_UB <- (1-(stepsize*abs(pt_UB-0.975))) * t_UB
        }
      }
      if(abs(pt_UB) >= abs(0.975+epsilon))  # we need smaler value of t_LB
      {
        if(t_UB>0) # we need smaller value of t_UB

        {t_UB <- (1-(stepsize*abs(pt_UB-0.975))) * t_UB}
        else # we need a larger value of T_UB
        {t_UB <- (1+(stepsize*abs(pt_UB-0.975))) * t_UB}

      }

    }
    pt_LB <- stats::pt(t_LB, df=df, ncp=t.unbiased)
    pt_UB <- stats::pt(t_UB, df=df, ncp=t.unbiased)
    i <- i + 1
  }

  switch(expDesign,
         CrossOverRM={
           d_RM_LB <- t_LB*sqrt((n1+n2)/(2*n1*n2))
           d_RM_UB <- t_UB*sqrt((n1+n2)/(2*n1*n2))
           d_IG_LB <- d_RM_LB*sqrt(1-r)
           d_IG_UB <- d_RM_UB*sqrt(1-r)
           g_RM_LB <- c*t_LB*sqrt((n1+n2)/(2*n1*n2))
           g_RM_UB <- c*t_UB*sqrt((n1+n2)/(2*n1*n2))
           g_IG_LB <- g_RM_LB*sqrt(1-r)
           g_IG_UB <- g_RM_UB*sqrt(1-r)

         },
         IG={
           d_IG_LB <- t_LB*sqrt((n1+n2)/(n1*n2))
           d_IG_UB <- t_UB*sqrt((n1+n2)/(n1*n2))
           g_IG_LB <- c*t_LB*sqrt((n1+n2)/(n1*n2))
           g_IG_UB <- c*t_UB*sqrt((n1+n2)/(n1*n2))
           d_RM_LB <- "NA" #d_IG_LB/sqrt(2*(1-r)) Not really applicable
           d_RM_UB <- "NA" #d_IG_UB/sqrt(2*(1-r)) Not really applicable
           g_RM_LB <- "NA" #d_IG_LB/sqrt(2*(1-r)) Not really applicable
           g_RM_UB <- "NA" #d_IG_UB/sqrt(2*(1-r)) Not really applicable
         },
         BeforeAfterRM={
           d_RM_LB <- t_LB/sqrt(n1)
           d_RM_UB <- t_UB/sqrt(n1)
           d_IG_LB <- d_RM_LB*sqrt(2*(1-r))
           d_IG_UB <- d_RM_UB*sqrt(2*(1-r))
           g_RM_LB <- c*t_LB/sqrt(n1)
           g_RM_UB <- c*t_UB/sqrt(n1)
           g_IG_LB <- g_RM_LB*sqrt(2*(1-r))
           g_IG_UB <- g_RM_UB*sqrt(2*(1-r))
         },
         stop("Need to specify the expDesign parameter\n.")
  )

  result <- list(t_LB=t_LB, t_UB=t_UB, d_RM_LB=d_RM_LB, d_RM_UB=d_RM_UB, d_IG_LB=d_IG_LB, d_IG_UB=d_IG_UB,
                 g_RM_LB=g_RM_LB, g_RM_UB=g_RM_UB, g_IG_LB=g_IG_LB, g_IG_UB=g_IG_UB, i=i)
  return(result) #t=14.4, n1=15, n2=15, r=0.6401, t_LB=10.64733, t_UB=19.47515,d_RM_LB=2.749129,  d_RM_UB=5.028462, d_IG_LB=1.649248, d_IG_UB=3.016658, g_RM_LB=2.674828, g_RM_UB= 4.892558, g_IG_LB= 1.604674, g_IG_UB=2.935127,i=18
}

