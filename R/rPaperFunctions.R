#' script to obtain correlation coefficients

#' @title CalculateRLevel1
#' @description This function calculates the r value for a 2-group (2G) or 4-Group (4G) Crossover experiment for each sequence group and each outcome metric. The function returns both the exact r value and the r value based on pooled variances for each sequence group and outcome metric
#' @author Barbara Kitchenham and Lech Madeyski
#' @export CalculateRLevel1
#' @param Dataset This holds the data for each participant in a 2-group or 4-group crossover experiment in the "wide" format. I.e., there is only one entry per participant. The data set should have been generated from a long version of the data based on a variable labelled "Period" which is used to define which participant data was collected in the first period of the experiment - see function ExtractLevel1ExperimentRData.
#' @param StudyID This holds an identifier used to identify the origin of the experimental data in the output from this function.
#' @param Groups This is a list that defined the sequence group identifiers used in the dataset.
#' @param Metrics This is a list of metrics, e.g., ("Correctness","Time","Efficiency").
#' @param ExperimentName This an identifiers used to define the specific experiment in the output from this function.
#' @param Type this is a character string specifying whether the experiment is a two sequence group of four sequence group experiment.
#' @param Control this is a character string that defines the control treatment in the experiment.
#' @return table this is a tibble holding information identifying for each metric and sequence group the first time period and second time period variance, the pooled variance, the variance of the difference values and the exact r and pooled r.
#' # importFrom stats
#' # importFrom var
#' # importFrom tibble
#' @examples
#' ExperimentNames=c("EUBAS","R1UCLM","R2UCLM","R3UCLM")
#' ShortExperimentNames=c("E1","E2","E3","E4")
#' Metrics=c("Comprehension","Modification")
#' Type=c("4G", "4G", "4G", "4G")
#' Groups=c("A","B","C","D")
#' StudyID="S2"
#' Control="SC"
#'# Obtain experimental data from a file and put in wide format
#' dataset2= KitchenhamEtAl.CorrelationsAmongParticipants.Scanniello14TOSEM
#' ReshapedData=ExtractExperimentData(dataset2, ExperimentNames=ExperimentNames,
#'   idvar="ParticipantID",timevar="Period",ConvertToWide=TRUE)
#' # Calculate the correlations for each sequence group and each metric.
#'  CalculateRLevel1(Dataset=ReshapedData[[1]], StudyID, Groups=c("A","B","C","D"),
#'    ExperimentName=ShortExperimentNames[1],Metrics,Type=Type[1],Control)
#' # A tibble: 8 x 15
#' # # A tibble: 8 x 15
#' # Study Exp   Group Metric Id        n ControlFirst    var1   var2
#' # <chr> <chr> <chr> <chr>  <chr> <int> <lgl>          <dbl>  <dbl>
#' #   1 S2    E1    A     Compr… S2E1A     6 FALSE        0.0183  0.0163
#' # 2 S2    E1    B     Compr… S2E1B     6 TRUE         0.0201  0.0326
#' # 3 S2    E1    C     Compr… S2E1C     6 FALSE        0.00370 0.0155
#' # 4 S2    E1    D     Compr… S2E1D     6 TRUE         0.0173  0.0201
#' # 5 S2    E1    A     Modif… S2E1A     6 FALSE        0.0527  0.0383
#' # 6 S2    E1    B     Modif… S2E1B     6 TRUE         0.0185  0.0482
#' # 7 S2    E1    C     Modif… S2E1C     6 FALSE        0.00655 0.0244
#' # 8 S2    E1    D     Modif… S2E1D     6 TRUE         0.0222  0.0266
#' # # … with 6 more variables: varp <dbl>, ControlVarProp <dbl>,
#' # #   VarProp <dbl>, vardiff <dbl>, r <dbl>, r.p <dbl>



CalculateRLevel1 = function(Dataset,
                            StudyID,
                            Groups = c("A", "B", "C", "D"),
                            ExperimentName,
                            Metrics,
                            Type,
                            Control) {
  if (Type == "2G")
    Groups = Groups[1:2]
  else
    Groups = Groups[1:4]

  NumReps = length(Groups)
  table = NULL

  r.table = NULL
  NumMets = length(Metrics)
  N = length(Dataset$SequenceGroup.1)
  rvalidreps = 0
  for (j in 1:NumMets) {
    for (i in 1:NumReps) {
      # Analyse the wide data for each Sequence Group in a specific experiment for a specific metric
      #ss = base::subset(Dataset, SequenceGroup.1 == Groups[i])
      #Fix LM, BAK Agreed:
      ss = base::subset(Dataset, Dataset$SequenceGroup.1 == Groups[i])
      # Select the data for the specific metric corresponding to each time period
      name1 = base::paste(Metrics[j], "1", sep = ".")
      name2 = base::paste(Metrics[j], "2", sep = ".")
      res1 = base::subset(ss, select = name1)
      n1 = length(ss$SequenceGroup.1)
      res2 = base::subset(ss, select = name2)
      n2 = length(ss$SequenceGroup.2) # n1 should equal n2
      if (n1 != n2)
        return("Incorrect subset")
      # Analyse the difference data for each Metric
      diff = res2 - res1

      if (n1 > 1) {
        vardiff = as.numeric(stats::var(diff))
        var1 = as.numeric(stats::var(res1))
        var2 = as.numeric(stats::var(res2))
      }
      else {
        # Cannot construct a sensible variance
        vardiff = 0
        var1 = 0
        var2 = 0
      }

	 if (vardiff > 0 & var1 > 0 & var2 > 0) {

       # Ensure all variances are non-zero otherwise do not add to the row

      # Identify which variance corresponds to the control condition
      ControlFirst = FALSE
      if (ss$Treatment.1[1] == Control)
        ControlFirst = TRUE
      if (ControlFirst)
        controlvar = var1
      else
        controlvar = var2

     	ControlVarProp = controlvar / (var1 + var2)
      	TimePeriodVarProp = var1 / (var1 + var2)

        # Calculate r without assuming variance homogeneity
        r = as.numeric((var1 + var2 - vardiff) / (2 * sqrt(var1 * var2)))

        #r = signif(r, 4)
        #  Calculate r assuming variance homogeneity
        pooled.var = (var1 + var2) / 2
        r.p = as.numeric((2 * pooled.var - vardiff) / (2 * pooled.var))
        #r.p = signif(r.p, 4)
        #var1 = signif(var1, 4)
        #var2 = signif(var2, 4)
        #vardiff = signif(vardiff, 4)
        #ControlVarProp = signif(ControlVarProp, 3)
        #TimePeriodVarProp = signif(TimePeriodVarProp, 3)
    # Construct an Id that is unique for a specific study, experiment and sequence group but treats r-values for different metrics as repeated measures.
        Id = base::paste(StudyID, ExperimentName, sep = "")
        Id = base::paste(Id, Groups[i], sep = "")
        row = base::cbind(
          Study = StudyID,
          Exp = ExperimentName,
          Group = Groups[i],
          Metric = Metrics[j],
          Id = Id,
          n = n1,
          ControlFirst = ControlFirst,
          var1 = var1,
          var2 = var2,
          varp = pooled.var,
          ControlVarProp = ControlVarProp,
          VarProp = TimePeriodVarProp,
          vardiff = vardiff,
          r = r,
          r.p = r.p
        )


        table = tibble::as_tibble(rbind(table, row))
      }
    }

  }

  	#Coerce table columns to correct format
	table$n=as.integer(table$n)
	table$ControlFirst=as.logical(table$ControlFirst)
	table$var1=as.numeric(table$var1)
	table$var2=as.numeric(table$var2)
	table$varp=as.numeric(table$varp)
	table$ControlVarProp=as.numeric(table$ControlVarProp)
	table$VarProp=as.numeric(table$VarProp)
	table$vardiff=as.numeric(table$vardiff)
	table$r=as.numeric(table$r)
	table$r.p=as.numeric(table$r.p)


  return(table)

}



###########################################################################################################
#' @title ExtractGroupSizeData
#' @description This function constructs a table identifying the number of participants in each sequence group for a set of experiments each of which used a crossover design.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export ExtractGroupSizeData
#' @param ExpDataWide this is a list of tibbles each comprising data from one experiment in its wide format
#' @param StudyID an identifer for the group of related experiments (i.e., a family).
#' @param ShortExperimentNames a list of character strings identifying each experiment.
#' @param Type A list identifying the type of crossover "2G" or "4G" for each experiment in the family
#' @param Groups a list of the terms used to specify sequence groups in the experiments.
#' @return A tibble containing the number of participants in each sequence group in each experiment.
#' @examples
#' ExperimentNames=c("EUBAS","R1UCLM","R2UCLM","R3UCLM")
#' ShortExperimentNames=c("E1","E2","E3","E4")
#' Metrics=c("Comprehension","Modification")
#' Type=c("4G", "4G", "4G", "4G")
#' Groups=c("A","B","C","D")
#' StudyID="S2"
#' Control="SC"
#' # Obtain experimental data from a file and put in wide format
#' dataset2= KitchenhamEtAl.CorrelationsAmongParticipants.Scanniello14TOSEM
#' ReshapedData=ExtractExperimentData(dataset2, ExperimentNames=ExperimentNames,
#'   idvar="ParticipantID", timevar="Period", ConvertToWide=TRUE
#'   )
#' ExtractGroupSizeData(ReshapedData, StudyID, ShortExperimentNames, Type, Groups=Groups)
#' # A tibble: 16 x 4
#' #  Study Exp   Group     n
#' #  <chr> <chr> <chr> <int>
#' #1 S2    Exp1  A         6
#' #2 S2    Exp1  B         6
#' #3 S2    Exp1  C         6
#' #4 S2    Exp1  D         6
#' #5 S2    Exp2  A         6
#' #6 S2    Exp2  B         6
#' #7 S2    Exp2  C         5
#' #8 S2    Exp2  D         5
#' #9 S2    Exp3  A         5
#'#10 S2    Exp3  B         5
#'#11 S2    Exp3  C         6
#'#12 S2    Exp3  D         6
#'#13 S2    Exp4  A         5
#'#14 S2    Exp4  B         5
#'#15 S2    Exp4  C         4
#'#16 S2    Exp4  D         4



ExtractGroupSizeData = function(ExpDataWide,
                                StudyID,
                                ShortExperimentNames,
								Type,
                                Groups=c("A","B","C","D"))
  {

  NumExp = length(ShortExperimentNames)

  GroupDataTable = NULL
  for (i in 1:NumExp) {

    # Find the wide format data for a specific experiment
    ExpData=as.data.frame(ExpDataWide[[i]])

    # Identify how many sequence groups there are in the specific experiment (either 2 or 4) for a specific crossover
   if (Type[i]=="4G")  NumGroups=4 else  NumGroups=2

    # Find the names used to distinguish the sequence groups
    ExpGroups=Groups[1:NumGroups]

    # Set up a variable to hold the number of participants in each sequence group
    n=c(NumGroups,NA)
    for (j in 1:NumGroups){
    	#Find the data for each sequence group
    	ss=base::subset(ExpData,ExpData$SequenceGroup.1==ExpGroups[j])
    	# Find the number of participants in the specific sequence group
    	n[j]=length(ss$ParticipantID)
    	# Put the data together to specify the Study, Experiment, Sequence group and number of participants
    	row= cbind(Study=StudyID,Exp=ShortExperimentNames[i],Group=ExpGroups[j],n=n[j])
        GroupDataTable = as.data.frame(rbind(GroupDataTable, row))
       }
  }


 ################# Change to tibble ##########################################
  GroupDataTable=tibble::as_tibble(GroupDataTable)

  GroupDataTable$n=as.integer(GroupDataTable$n)
  GroupDataTable$Study=as.character(GroupDataTable$Study)
  GroupDataTable$Exp=as.character(GroupDataTable$Exp)
  GroupDataTable$Group=as.character(GroupDataTable$Group)
  return(GroupDataTable)
}

###########################################################################################

#' @title ConstructLevel1ExperimentRData
#' @description This function returns the r value for a 2-group (2G) or 4-Group (4G) Crossover experiment for a group of 1 or more experiments for each sequence group and each outcome metric. For sets of 2 or more experiments, the experiments are assumed to be replicates and to report the same sets of Metrics and have the same Control treatment and use the same squence Group identifiers, but are not necessarily the same Type. We return both the exact r value and the r value based on pooled variances for each sequence group and outcome metric.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export ConstructLevel1ExperimentRData
#' @param Data This is a list parameter each entry in the list holds the data for each participant in a 2-group or 4-group crossover experiment in the "wide" format. I.e., there is only one entry per participant. The data should have been generated from a long version of the data based on a variable labelled "Period" which is used to define which participant data was collected in the first period of the experiment - see function ExtractLevel1ExperimentRData.
#' @param StudyID This holds an identifer used to identify the origin of the experimental data in the output from this function.
#' @param Groups This is a list that defined the sequence group identifiers used in the dataset.
#' @param ExperimentNames This a list of identifiers used to define each experiment in the output from this function.
#' @param Metrics This is a list of of character strings identifying each outcome metric reported in each of the experiments in the set of replicated experiments.
#' @param Type this is a list of character strings specifying for each experiment whether the experiment is a 2-group or 4-group experiment
#' @param Control this is a character string that defines the control treatment in the experiment.
#' @return R.Data.Table this is a tibble holding information identifying for each metric and sequence group the first time period and second time period variance, the pooled variance, the variance of the difference values and the exact r and pooled r.
#' @examples
#' #
#' ShortExperimentNames=c("E1","E2","E3","E4")
#' FullExperimentNames=c("EUBAS","R1UCLM","R2UCLM","R3UCLM")
#' Metrics=c("Comprehension","Modification")
#' Groups=Groups=c("A","B","C","D")
#' Type=c(rep("4G",4))
#' StudyID="S2"
#' Control="SC"
#'# Obtain experimental data from each file and put in wide format
#' dataset2= KitchenhamEtAl.CorrelationsAmongParticipants.Scanniello14TOSEM
#' ReshapedData=ExtractExperimentData(dataset2, ExperimentNames=FullExperimentNames,
#'   idvar="ParticipantID",timevar="Period",ConvertToWide=TRUE)
#' # Calculate the correlations for each sequence group and each metric in each experiment
#' ConstructLevel1ExperimentRData(Data=ReshapedData, StudyID=StudyID,
#'   ExperimentNames=ShortExperimentNames, Groups=Groups, Metrics=Metrics, Type=Type,
#'   Control=Control
#' )
#' # # A tibble: 32 x 15
#' # Study Exp   Group Metric Id        n ControlFirst    var1   var2    varp
#' # <chr> <chr> <chr> <chr>  <chr> <int> <lgl>          <dbl>  <dbl>   <dbl>
#' #   1 S2    E1    A     Compr… S2E1A     6 FALSE        0.0183  0.0163 0.0173
#' # 2 S2    E1    B     Compr… S2E1B     6 TRUE         0.0201  0.0326 0.0263
#' # 3 S2    E1    C     Compr… S2E1C     6 FALSE        0.00370 0.0155 0.00962
#' # 4 S2    E1    D     Compr… S2E1D     6 TRUE         0.0173  0.0201 0.0187
#' # 5 S2    E1    A     Modif… S2E1A     6 FALSE        0.0527  0.0383 0.0455
#' # 6 S2    E1    B     Modif… S2E1B     6 TRUE         0.0185  0.0482 0.0333
#' # 7 S2    E1    C     Modif… S2E1C     6 FALSE        0.00655 0.0244 0.0155
#' # 8 S2    E1    D     Modif… S2E1D     6 TRUE         0.0222  0.0266 0.0244
#' # 9 S2    E2    A     Compr… S2E2A     6 FALSE        0.0194  0.0425 0.0309
#' # 10 S2    E2    B     Compr… S2E2B     6 TRUE         0.0198  0.0192 0.0195
#' # # … with 22 more rows, and 5 more variables: ControlVarProp <dbl>,
#' # #   VarProp <dbl>, vardiff <dbl>, r <dbl>, r.p <dbl>

ConstructLevel1ExperimentRData = function(Data,
                                          StudyID,
                                          ExperimentNames,
                                          Groups,
                                          Metrics,
                                          Type,
                                          Control) {
  # This calls the algorithm that constructs the correlation values for each metric and each independent sequence group in an experiment. It identifies the data for each sequence group in each experiment for each experiment and collates the reurned r-values for each sequence group


  NumExp = length(ExperimentNames)

  R.Data.Table = NULL

  for (i in 1:NumExp) {
    r.data = CalculateRLevel1(
      Data[[i]],
      Groups = Groups,
      StudyID = StudyID,
      ExperimentName = ExperimentNames[i],
      Metrics,
      Type = Type[i],
      Control = Control
    )

    R.Data.Table = tibble::as_tibble(rbind(R.Data.Table, r.data))
  }

  return(R.Data.Table)
}

#############################################################################################################

#' @title ExtractExperimentData
#' @description This function reads datasets from a defined directory in the reproducer package that hold the results of a family crossover experiments in the long format. It converts the data to the wide format if required.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export ExtractExperimentData
#' @param DataSet This is a tibble holding the data for each crossover experiment in a family (a family can include only one experiment).
#' @param ExperimentNames This is a list with the full names of each experiment.
#' @param idvar This is the name of the column that contains the data for specific participants. It is only assumed to be unique within an experiment (default idvar="ParticipantID").
#' @param timevar This is the name of the table column that defines which data was collected in a specific time period. This function assumes that there are only two time periods (default timevar="Period").
#' @param ConvertToWide This determine whether the function converts the data to the wide format (default ConvertToWide=TRUE).
#' @return A list with an entry for the data for each experiment. If ConvertToWide is TRUE, it returns the data in the wide format otherwise it returns the data as it was read. Within each list item the data is returned as a tibble
#' #importFrom stats
#' # importFrom tibble
#' # importFrom base
#' @examples
#' ExperimentNames=c("EUBAS","R1UCLM","R2UCLM","R3UCLM")
#' Metrics=c("Comprehension","Modification")
#' Groups=c("A","B","C","D")
#' Type=c(rep("4G",4))
#' StudyID="S2"
#' Control="SC"
#' # Obtain experimental data from each file and put in wide format
#' dataset2= KitchenhamEtAl.CorrelationsAmongParticipants.Scanniello14TOSEM
#' ReshapedData = ExtractExperimentData(dataset2, ExperimentNames=ExperimentNames,
#'   idvar="ParticipantID",timevar="Period",ConvertToWide=TRUE)
#' ReshapedData[[1]]
#'
#'# A tibble: 24 x 15
#'# ParticipantID ExperimentID.1 SequenceGroup.1 System.1 Treatment.1 Comprehension.1
#'# <fct>         <fct>          <fct>           <fct>    <fct>                 <dbl>
#'#   1 1             EUBAS          A               S1       AM                     0.77
#'# 2 5             EUBAS          A               S1       AM                     0.61
#'# 3 9             EUBAS          A               S1       AM                     0.61
#'# 4 13            EUBAS          A               S1       AM                     0.52
#'# 5 17            EUBAS          A               S1       AM                     0.43
#'# 6 21            EUBAS          A               S1       AM                     0.77
#'# 7 2             EUBAS          B               S1       SC                     0.92
#'# 8 6             EUBAS          B               S1       SC                     0.63
#'# 9 10            EUBAS          B               S1       SC                     0.51
#'# 10 14            EUBAS          B               S1       SC                     0.64
#'# … with 14 more rows, and 9 more variables: Modification.1 <dbl>, CrossOverID.1 <fct>,
#'#   ExperimentID.2 <fct>, SequenceGroup.2 <fct>, System.2 <fct>, Treatment.2 <fct>,
#'#   Comprehension.2 <dbl>, Modification.2 <dbl>, CrossOverID.2 <fct>


ExtractExperimentData = function(DataSet,
                                 ExperimentNames,
                                 idvar = "ParticipantID",
                                 timevar = "Period",
                                 ConvertToWide = TRUE) {
  # This algorithm reads the data for each experiment in a family (or a study)  and converts to the "wide" data format

  NumExp = length(ExperimentNames)
  ExpData = list()


  for (i in 1:NumExp) {
    # Read each of the data sets into a separate list entry

		#tempdata =as.data.frame(base::subset(DataSet,ExperimentID==ExperimentNames[i]))
    #Fix LM BAK agreed :
    tempdata =as.data.frame(base::subset(DataSet,DataSet$ExperimentID==ExperimentNames[i]))

# Find the wide version of the dataset based on the Period variable if required
		if (ConvertToWide) {
			tempdata=stats::reshape(tempdata,idvar=idvar,timevar=timevar,direction="wide")}
# Return the data in the tibble format
		ExpData[[i]] =tibble::as_tibble(tempdata)

	}

	return(ExpData)
}

#######################################################################################################
#' @title CalculateLevel2ExperimentRData
#' @description This function analyses data on r values obtained in the format obtained from  the ConstructLevel1ExperimentRData function and finds the r-value for each metric for each experiment.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export CalculateLevel2ExperimentRData
#' @param Level1Data a tibble in the format produced by the ConstructLevel1ExperimentRData function which has r-values for each sequence group in a crossover experiment
#' @param Groups This is a list that defines the sequence group labels used in the dataset.
#' @param StudyID This holds an identifer used to identify the origin of the experimental data in the output from this function.
#' @param ExperimentNames This a list of identifiers used to define each experiment in the output from this function.
#' @param Metrics This is a list of of character strings identifying each outcome metric reported in each of the experiments in the set of replicated experiments.
#' @param Type this is a list of character strings specifying for each experiment whether the experiment is a two sequence group "2G" or four sequence group "4G" experiment.
#' return RExp.Table This is a table containing the pooled data variance and the pooled difference variance for the experiment and the value r for the experiment for eachm metric
#' @examples
#' ShortExperimentNames=c("E1","E2","E3","E4")
#' FullExperimentNames=c("EUBAS","R1UCLM","R2UCLM","R3UCLM")
#' Metrics=c("Comprehension","Modification")
#' Groups=c("A","B","C","D")
#' Type=c(rep("4G",4))
#' StudyID="S2"
#' Control="SC"
#' # Obtain experimental data from each file and put in wide format
#' ReshapedData = ExtractExperimentData(
#'   KitchenhamEtAl.CorrelationsAmongParticipants.Scanniello14TOSEM,
#'   ExperimentNames=FullExperimentNames, idvar="ParticipantID",timevar="Period",ConvertToWide=TRUE)
#' Lev1Data= ConstructLevel1ExperimentRData(ReshapedData, StudyID, ShortExperimentNames, Groups,
#'   Metrics, Type, Control)
#' CalculateLevel2ExperimentRData(Lev1Data,Groups=Groups,StudyID=StudyID,
#'   ExperimentNames=ShortExperimentNames,Metrics=Metrics,Type=Type)
#' # A tibble: 8 x 10
#' #  StudyID ExpID     N Metric        PooledVar1 PooledVar2 VarProp PooledVar PooledDiffVar    r.Exp
#' #  <chr>   <chr> <int> <chr>              <dbl>      <dbl>   <dbl>     <dbl>         <dbl>    <dbl>
#' # 1 S2      S2E1     24 Comprehension     0.0148     0.0212   0.412    0.0180        0.0248  0.311
#' # 3 S2      S2E2     22 Comprehension     0.0487     0.0224   0.684    0.0356        0.0534  0.250
#' # 4 S2      S2E2     22 Modification      0.0445     0.0266   0.626    0.0356        0.0628  0.117
#' # 5 S2      S2E3     22 Comprehension     0.0353     0.0402   0.467    0.0377        0.105  -0.391
#' # 6 S2      S2E3     22 Modification      0.0433     0.0414   0.511    0.0424        0.0997 -0.176
#' # 7 S2      S2E4     18 Comprehension     0.0439     0.0237   0.649    0.0338        0.0355  0.475
#' # 8 S2      S2E4     18 Modification      0.0322     0.0592   0.353    0.0457        0.0894  0.0222


CalculateLevel2ExperimentRData = function(Level1Data,
                                          Groups,
                                          StudyID,
                                          ExperimentNames,
                                          Metrics,
                                          Type) {
  NumExp = length(ExperimentNames)

  RExp.Table = NULL


  NumMets = length(Metrics)


  for (i in 1:NumExp) {
    # Extract the variance data for a specific experiment
    # LM fix 1: ExpData = base::subset(Level1Data, Exp == ExperimentNames[i])
    # BAK Checked and agreed
    ExpData = base::subset(Level1Data, Level1Data$Exp == ExperimentNames[i])


    if (Type[i] == "2G")
      SequenceGroups = Groups[1:2]
    else
      SequenceGroups = Groups[1:4]

    NumGroups = length(SequenceGroups)
    ExpID = ExpData$Exp[1]
    ID = base::paste(StudyID, ExpID, sep = "")
    for (j in 1:NumMets) {
      # Find the sequence group statistics for a specific Metric
      # LM fix 2: MetData = base::subset(ExpData, Metric == Metrics[j])
      # Wrong
      # BAK Fix 3: MetData = base::subset(ExpData, Level1Data$Metric == Metrics[j])

     MetData = base::subset(ExpData, ExpData$Metric == Metrics[j])

      # Find the number of participants in the experiment
      NumParticipants = sum(MetData$n)

      # Calculate the pooled sequence group variance and difference variance
      PooledVar1 = sum(MetData$var1 * (MetData$n - 1)) / (NumParticipants -
                                                            NumGroups)
      PooledVar2 = sum(MetData$var2 * (MetData$n - 1)) / (NumParticipants -
                                                            NumGroups)
      VarProp = PooledVar1 / (PooledVar1 + PooledVar2)
      PooledVar = sum(MetData$varp * (MetData$n - 1)) / (NumParticipants -
                                                           NumGroups)
      PooledDiffVar = sum(MetData$vardiff * (MetData$n - 1)) / (NumParticipants -
                                                                  NumGroups)
      # Calculate r
      r = (2 * PooledVar - PooledDiffVar) / (2 * PooledVar)
      # Tidy up calculated values - premature restrictions may cause rounding errors
      #r = signif(r, 4)
      # PooledVar1 = signif(PooledVar1, 4)
      # PooledVar2 = signif(PooledVar2, 4)
      # PooledVar = signif(PooledVar, 4)
      #VarProp = signif(VarProp, 4)
      #PooledDiffVar = signif(PooledDiffVar, 4)
      # Put the output into a format that identifes the study, experiment, experiment size and metric and the average r for the experiment
      row =
        cbind(
          PooledVar1 = PooledVar1,
          PooledVar2 = PooledVar2,
          VarProp = VarProp,
          PooledVar = PooledVar,
          PooledDiffVar = PooledDiffVar,
          r.Exp = r
        )


      RExp.Table = rbind(
        RExp.Table,
        cbind(
          StudyID = StudyID,
          ExpID = ID,
          N = NumParticipants,
          Metric = Metrics[j],
          row
        )
      )

    }

  }

	RExp.Table=tibble::as_tibble(RExp.Table)

	# Coerce the data items to the correct formats
	RExp.Table$N=as.integer(RExp.Table$N)
	RExp.Table$PooledVar1=as.numeric(RExp.Table$PooledVar1)
	RExp.Table$PooledVar2=as.numeric(RExp.Table$PooledVar2)
	RExp.Table$PooledVar=as.numeric(RExp.Table$PooledVar)
	RExp.Table$VarProp=as.numeric(RExp.Table$VarProp)
	RExp.Table$PooledDiffVar=as.numeric(RExp.Table$PooledDiffVar)
	RExp.Table$r.Exp=as.numeric(RExp.Table$r.Exp)
	return(RExp.Table)

}

############################################################################################################

# Functions to help produce tables for the paper

#' @title ExtractSummaryStatisticsRandomizedExp
#' @description This function extracts data obtained from the lme4 package lmer function. It assumes a simple randomized experiment with each element having one or more repeated measures. It outputs the mean together with its standard error and confidence interval bounds.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export ExtractSummaryStatisticsRandomizedExp
#' @param lmeRA The output from the lmer function
#' @param N The total number of observations
#' @param alpha the probability level to be used when constructing the confidence interval bounds.
#' @return REA.Summary A dataframe holding the number of observations N, the overall mean value as
#' its standard error reported as by the lmer function, and its confidence interval bounds.
#' @examples
#' ShortExperimentNames=c("E1","E2","E3","E4")
#' FullExperimentNames=c("EUBAS","R1UCLM","R2UCLM","R3UCLM")
#' Metrics=c("Comprehension","Modification")
#' Groups=c("A","B","C","D")
#' Type=c(rep("4G",4))
#' StudyID="S2"
#' Control="SC"
#' ReshapedData= ExtractExperimentData(
#'   KitchenhamEtAl.CorrelationsAmongParticipants.Scanniello14TOSEM,
#'   ExperimentNames=FullExperimentNames, idvar="ParticipantID", timevar="Period",
#'   ConvertToWide=TRUE
#' )
#' NewTable= ConstructLevel1ExperimentRData(ReshapedData, StudyID, ShortExperimentNames, Groups,
#'   Metrics,Type,Control
#'   )
#' resRe=lme4::lmer(r~(1|Id),data=NewTable)
#' summary(resRe)
#' # Linear mixed model fit by REML ['lmerMod']
#' # Formula: r ~ (1 | Id)
#' # REML criterion at convergence: 47.8
#' # Scaled residuals:
#' #    Min      1Q  Median      3Q     Max
#' # -1.4382 -0.9691  0.2190  0.8649  1.4761
#' #
#' # Random effects:
#' #  Groups   Name        Variance Std.Dev.
#' #   Id       (Intercept) 0.03978  0.1994
#' #   Residual             0.20974  0.4580
#' #  Number of obs: 32, groups:  Id, 16
#' #
#' #  Fixed effects:
#' #             Estimate Std. Error t value
#' #  (Intercept)  0.06175    0.09508   0.649
#' #  N=length(NewTable$r)
#'  ExtractSummaryStatisticsRandomizedExp(lmeRA=resRe,N=32,alpha=0.05)
#' #      N    Mean      SE LowerBound UpperBound
#' #   1 32 0.06175 0.09508    -0.1319     0.2554


ExtractSummaryStatisticsRandomizedExp = function(lmeRA, N, alpha = 0.05) {
  REA.fe.t = stats::coef(summary(lmeRA))
  REA.Mean = REA.fe.t[1, "Estimate"]
  REA.SE = REA.fe.t[1, "Std. Error"]
  # Calculate confidence interval for mean r
  vv = stats::qt(alpha / 2, N)
  REA.LowerBound = REA.Mean + vv * REA.SE
  REA.UpperBound = REA.Mean - vv * REA.SE

  # Truncate the values for readability
  REA.Mean = signif(REA.Mean, 4)
  REA.SE = signif(REA.SE, 4)
  REA.LowerBound = signif(REA.LowerBound, 4)
  REA.UpperBound = signif(REA.UpperBound, 4)
  REA.Summary = data.frame(
    cbind(
      N = N,
      Mean = REA.Mean,
      SE = REA.SE,
      LowerBound = REA.LowerBound,
      UpperBound = REA.UpperBound
    )
  )

  return(REA.Summary)
}
########################################################################################################################
#' @title calculateBasicStatistics
#' @description This function calculates the following statistcs for a set of data: length, mean, median, variance, standard error of the mean, and confidence interval bounds. The input data must be a vector of 2 or more numerical values.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export calculateBasicStatistics
#' @param x The data to be summarized
#' @param alpha The probability level to be used when constructing the confidence interval bounds.
#' @return A dataframe comprising the length, mean, variance, standard error and confidence limit bounds of the input data x.
#' ShortExperimentNames=c("E1","E2","E3","E4")
#' FullExperimentNames=c("EUBAS","R1UCLM","R2UCLM","R3UCLM")
#' Groups=c("A","B","C","D")
#' Type=c(rep("4G",4))
#' StudyID="S2"
#' Control="SC"
#' ReshapedData= ExtractExperimentData(KitchenhamEtAl.CorrelationsAmongParticipants.Scanniello14TOSEM,ExperimentNames=FullExperimentNames, idvar="ParticipantID",timevar="Period",ConvertToWide=TRUE)
#' NewTable= ConstructLevel1ExperimentRData(ReshapedData,StudyID,ShortExperimentNames,Groups,Metrics,Type,Control)
#' calculateBasicStatistics(NewTable$r)
#'  #    N    Mean Median Variance      SE LowerBound UpperBound
#'  # 1 32 0.06175 0.1688   0.2482 0.08808    -0.1109     0.2344

calculateBasicStatistics = function(x, alpha = 0.05) {
  # Finds the basic statistics for a set of data
  N = length(x)
  mean_x = mean(x)
  var_x = stats::var(x)
  median_x = stats::median(x)
  se_x = sqrt(var_x / N)
  crit <- stats::qnorm(alpha / 2)
  lowerbound_x = mean_x + crit * se_x
  upperbound_x = mean_x - crit * se_x
  summary_x = data.frame(
    cbind(
      N = N,
      Mean = mean_x,
      Median = median_x,
      Variance = var_x,
      SE = se_x,
      LowerBound = lowerbound_x,
      UpperBound = upperbound_x
    )
  )
  summary_x = signif(summary_x, 4)
  return(summary_x)
}
############################################################################################################


#' @title calculateGroupSummaryStatistics
#' @description This function calculates the following statistics data within groups: length, mean, median, variance, standard error of the mean, and confidence interval bounds.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export calculateGroupSummaryStatistics
#' @param x The data to be summarized. This must be a vector of 2 or more numerical values
#' @param Group The categorical data data defining the groups. This must vector of the same length as x containing factors specifying the data groups
#' @return A dataframe comprising the number, mean, variance, standard error and confidence limit bounds of the data in each category
#' @examples
#' ShortExperimentNames=c("E1","E2","E3","E4")
#' FullExperimentNames=c("EUBAS","R1UCLM","R2UCLM","R3UCLM")
#' Metrics=c("Comprehension","Modification")
#' Groups=c("A","B","C","D")
#' Type=c(rep("4G",4))
#' StudyID="S2"
#' Control="SC"
#' ReshapedData= ExtractExperimentData(
#'   KitchenhamEtAl.CorrelationsAmongParticipants.Scanniello14TOSEM,
#'   ExperimentNames=FullExperimentNames, idvar="ParticipantID",timevar="Period",
#'   ConvertToWide=TRUE
#' )
#' NewTable= ConstructLevel1ExperimentRData(ReshapedData, StudyID,
#'   ShortExperimentNames, Groups, Metrics, Type, Control
#' )
#' SeqGroupLev=NULL
#' N.NT=length(NewTable$r)
#' for (i in 1:N.NT) {
#' 	if (NewTable$n[i]<=8) SeqGroupLev[i]=as.character(NewTable$n[i])
#' 	if (NewTable$n[i]>8) SeqGroupLev[i]=as.character(9)
#'  }
#'  calculateGroupSummaryStatistics(NewTable$r,Group=SeqGroupLev)
#' #     N    Mean  Median Variance  StDev     SE
#' #  1  4 -0.0833 -0.1699   0.2314 0.4810 0.2405
#' #  2 12  0.3658  0.4477   0.2109 0.4592 0.1326
#' #  3 16 -0.1300 -0.2214   0.1933 0.4397 0.1099


calculateGroupSummaryStatistics = function(x, Group) {
  # Calculates summary statistics for grouped data
  agg.mean = stats::aggregate(x, by = list(Group), FUN = mean)
  agg.mean = reshape::rename(agg.mean, c(x = "Mean"))
  agg.median = stats::aggregate(x, by = list(Group), FUN = stats::median)
  agg.median = reshape::rename(agg.median, c(x = "Median"))
  agg.length = stats::aggregate(x, by = list(Group), FUN = length)
  agg.length = reshape::rename(agg.length, c(x = "N"))
  agg.var = stats::aggregate(x, by = list(Group), FUN = stats::var)
  agg.var = reshape::rename(agg.var, c(x = "Variance"))
  agg.StDev = sqrt(agg.var$Variance)
  agg.SE = sqrt(agg.var$Variance / agg.length$N)

  summaryTable = NULL
  summaryTable = data.frame(
    cbind(
      N = agg.length$N,
      Mean = agg.mean$Mean,
      Median = agg.median$Median,
      Variance = agg.var$Variance,
      StDev = agg.StDev,
      SE = agg.SE
    )
  )
  summaryTable = signif(summaryTable, 4)

  return(summaryTable)
}

############################################################################################################################

#' @title rSimulations
#' @description This function simulates many datasest from the same bivariate distribution to
#' investigate the distribution of correlations for specific sample sizes.
#' @author Barbara Kitchenham and Lech Madeyski
#' @export rSimulations
#' @param mean The mean used for one of bivariate distributions - assumed to be the control condition in an experiment.
#' @param var The variance used for both treatment groups. It must be a real value greater than 0.
#' @param diff This value is added to the parameter mean to specify the mean for the other bivariate distribution - assumed to be the treatment condition in an experiment.
#' @param r This specifies the correlation coefficient to be used for the bivariate normal distribution it must be a value in the range [-1,1].
#' @param N The number of observations in each simulated bivariate normal data set.
#' @param VarAdj This value will be added to the variance of the treatment condition.
#' @param reps The number of bivariate data sets that will be simulated.
#' @param seed This specifies the seed value for the simulations and allows the experiment to be repeated.
#' @param returntSignificant If set to true the percentage of times the t-test delivered a value significant at the 0.05 level is reported (default returntSignificant=F).
#' @param returndata If set to FALSE, the function returns the summary information across all the replications  (default returndata=F). If set to TRUE the function outputs the r and variance ratio, and variance accuracy values generated in each replication.
#' @param plothist If set to T, the function outputs a histogram of the r-values, the varprop values and the accuracy values (default plothist=F).
#' @return output If returndata=F, the output returns summary information about the average of r and the variance properties across the replicated data sets. If returndata=T, the function returns the r-values obtained for each of the simulated data sets to gather with the variance ratio, the variance accuracy measure and a dummy variable indicating whether a test of significance between the mean values was significant (which is indicated by the dummy variable being set to 1) or not (which is indicated by the dummy variable being set to 0)
#' @examples
#' # output=rSimulations(mean=0,var=1,diff=0,r=0.25,N=4,reps=10000)
#' # reduced reps to pass CRAN time limits
#' output=rSimulations(mean=0,var=1,diff=0,r=0.25,N=4,reps=1000)
#' output=signif(output,4)
#' output
#' #  r.Mean r.Median  Var.r PercentNegative Mean.VarProp Variance.VarProp ...
#' # 1 0.2132   0.3128 0.3126           34.21       0.5036          0.06046 ...
#' #output=rSimulations(mean=0,var=1,diff=0.8,r=0.25,N=60,reps=10000,returntSignificant=TRUE)
#' # reduced reps to pass CRAN time limits
#' output=rSimulations(mean=0,var=1,diff=0.8,r=0.25,N=60,reps=1000,returntSignificant=TRUE)
#' output=signif(output,4)
#' output
#' #   r.Mean r.Median   Var.r PercentNegative Mean.VarProp Variance.VarProp ...
#' # 1 0.2492   0.2534 0.01529            2.62       0.5009         0.003897 ...
#' output=rSimulations(mean=0,var=1,diff=0,r=0.25,N=30,reps=10,returndata=TRUE)
#' output
#' #     rvalues   VarProp VarAccuracy VarDiffAccuracy tSig
#' #1  0.3981111 0.4276398   0.8630528       0.6974386    0
#' #2  0.2104742 0.4994285   0.7812448       0.8224174    0
#' #3  0.4252424 0.4933579   1.1568545       0.8866058    0
#' #4  0.3502651 0.6004373   0.8710482       0.7628923    0
#' #5  0.3845145 0.6029086   0.9618363       0.7998859    0
#' #6  0.1397217 0.4201069   1.1817022       1.3582855    0
#' #7  0.2311455 0.3894894   0.8322239       0.8594886    0
#' #8  0.3725047 0.5985897   1.1742117       0.9938662    0
#' #9  0.4881618 0.2712268   0.7585261       0.5723671    0
#' #10 0.1568071 0.3936400   0.9869924       1.1143561    0


rSimulations = function(mean,
                        var,
                        diff,
                        r,
                        N,
                        reps,
                        VarAdj = 0,
                        seed = 123,
                        returntSignificant = F,
                        returndata = F,
                        plothist = F) {
  # This function simulates bivariate normal distributions in order to investigate the distribution of given r-values for different sample sizes.


  set.seed(seed)

  # Set up variables to hold the data from each replication

  rvalues = c(rep(NA, reps))
  varratio = c(rep(NA, reps))
  varaccuracy = c(rep(NA, reps))
  vardiffaccuracy=c(rep(NA,reps))
  tSig = c(rep(NA, reps))

  rnegative = 0

  # Set up the parameters for the bivariate normal distribution
  meanvec = c(mean + diff, mean)
  covar = r * sqrt(var) * sqrt(var + VarAdj)
  sigma = matrix(c(var + VarAdj, covar, covar, var),
                 nrow = 2,
                 ncol = 2)


  for (i in 1:reps) {
    # Generate and analyse each simulation
    Mydataraw = MASS::mvrnorm(N, meanvec, sigma)
    Mydataraw = as.data.frame(Mydataraw)
    names(Mydataraw) = c("y", "x")

    # Find the variance and r values
    varx = stats::var(Mydataraw$x)
    vary = stats::var(Mydataraw$y)
    Mydatadiff = Mydataraw$y - Mydataraw$x

    vardiff=stats::var(Mydatadiff)

# Measure the heterogeneity of the between participants variance
    varratio[i] = varx / (varx + vary)

 # Measure the accuracy of the between participants variance
    varaccuracy[i] = (varx + vary) / (2 * var + VarAdj)

 # Measure the accuracy of the difference data variance
	 vardiffaccuracy[i]=vardiff/((2*var+VarAdj)*(1-r))

    rvalues[i] = (varx + vary - vardiff) / (2 * sqrt(varx * vary))

    if (rvalues[i] <= 0)
      rnegative = rnegative + 1

    # Do a t-test of the difference between the means
    ttest.results = stats::t.test(Mydataraw$x, Mydataraw$y)
    tSig[i] = if (ttest.results$p.value < 0.05)
      1
    else
      0

  }

  rnegative = 100 * rnegative / reps
  raverage = mean(rvalues)
  rmedian = stats::median(rvalues)
  Varr = stats::var(rvalues)
  meanvarrat = mean(varratio)
  varvarrat = stats::var(varratio)

  meanvaracc = mean(varaccuracy)
  varvaracc = stats::var(varaccuracy)

  meanvardiffacc=mean(vardiffaccuracy)
  varvardiffacc=stats::var(vardiffaccuracy)


  PercentSig = 100 * mean(tSig)

  if (plothist)
  {
    graphics::par(mfrow = c(4, 2))
    graphics::hist(rvalues, main = "Correlation Histogram", xlab = "Correlations")
    graphics::boxplot(rvalues, main = "Correlation Boxplot")

    graphics::hist(varratio, main = "Variance Ratio Histogram", xlab = "Variance Ratio Values")
    graphics::boxplot(varratio, main = "Variance Ratio Boxplot")


    graphics::hist(varaccuracy, main = "Variance Accuracy Histogram", xlab = "Variance Accuracy Values")
    graphics::boxplot(varaccuracy, main = "Variance Accuracy Boxplot")

    	graphics::hist(vardiffaccuracy,main="Diff Var Accuracy Histogram",xlab="Diff Var Accuracy Values")
	graphics::boxplot(vardiffaccuracy, main="Diff Var Accuracy Boxplot")

  }

  count = 0
  TheoreticalVarRatio = var / (2*var + VarAdj)
  LowerVarRatio = TheoreticalVarRatio / 2
  UpperVarRatio = TheoreticalVarRatio + TheoreticalVarRatio / 2
  for (i in 1:reps)
  {
    if (varratio[i] < LowerVarRatio |
        varratio[i] > UpperVarRatio)
      count = count + 1
  }
  acccount = 0
  for (i in 1:reps)
  {
    if (varaccuracy[i] < 0.5 | varaccuracy[i] > 1.5)
      acccount = acccount + 1
  }

  diffacccount = 0
  for (i in 1:reps)
  {
    if (vardiffaccuracy[i] < 0.5 | vardiffaccuracy[i] > 1.5)
      diffacccount = diffacccount + 1
  }

  if (!returndata)
  {
    output = data.frame(
      r.Mean = raverage,
      r.Median = rmedian,
      Var.r = Varr,
      PercentNegative = rnegative,
      Mean.VarProp = meanvarrat,
      Variance.VarProp = varvarrat,
      Percent.VarProp.Anomalies = 100 * count / reps,
      Mean.Varaccuracy = meanvaracc,
      Var.varaccuracy = varvaracc,
      Percent.VarAccuracy.Anomalies = 100 * acccount / reps,
      Mean.VarDiffAccuracy=meanvardiffacc,
      Var.VarDiffAccuracy= varvardiffacc,
      Percent.DiffVarAccuracy.Anomalies=100*diffacccount / reps
    )
    if (returntSignificant)
      output = data.frame(cbind(output, PercentSig = PercentSig))

  }
  else {
    output = data.frame(cbind(
      rvalues = rvalues,
      VarProp = varratio,
      VarAccuracy = varaccuracy,
      VarDiffAccuracy=vardiffaccuracy,
      tSig = tSig
    ))
  }
  return(output)

}
