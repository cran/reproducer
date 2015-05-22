#' Madeyski15SQJ.NDC data
#'
#' If you use this data set please cite: Lech Madeyski and Marian Jureczko, "Which Process Metrics Can Significantly Improve Defect Prediction Models? An Empirical Study," Software Quality Journal, 2015. DOI: 10.1007/s11219-014-9241-7
#'
#' "This paper presents an empirical evaluation in which several process metrics were investigated in order to identify the ones which significantly improve the defect prediction models based on product metrics. Data from a wide range of software projects (both, industrial and open source) were collected. The predictions of the models that use only product metrics (simple models) were compared with the predictions of the models which used product metrics, as well as one of the process metrics under scrutiny (advanced models). To decide whether the improvements were significant or not, statistical tests were performed and effect sizes were calculated. The advanced defect prediction models trained on a data set containing product metrics and additionally Number of Distinct Committers (NDC) were significantly better than the simple models without NDC, while the effect size was medium and the probability of superiority (PS) of the advanced models over simple ones was high (p=.016, r=-.29, PS=.76), which is a substantial finding useful in defect prediction. A similar result with slightly smaller PS was achieved by the advanced models trained on a data set containing product metrics and additionally all of the investigated process metrics (p=.038, r=-.29, PS=.68). The advanced models trained on a data set containing product metrics and additionally Number of Modified Lines (NML) were significantly better than the simple models without NML, but the effect size was small (p=.038, r=.06). Hence, it is reasonable to recommend the NDC process metric in building the defect prediction models." [http://dx.doi.org/10.1007/s11219-014-9241-7]
#'
#'
#' @format A data frame with variables:
#' \describe{
#' \item{Project}{In case of open source projects this field includes the name of the project as well as its version. In case of industrial projects this field includes the string "properietary" (we were not allowed to disclose the names of the analyzed industrial software projects developed by Capgemini Polska).}
#' \item{simple}{The percentage of classes that must be tested in order to find 80\% of defects in case of simple defect prediction models, i.e., using only software product metrics as predictors.}
#' \item{advanced}{The percentage of classes that must be tested in order to find 80\% of defects in case of advanced defect prediction models, using not only software product metrics but also the NDC (Number of distinct committers) process metric.}
#' }
#'
#' @source \url{http://madeyski.e-informatyka.pl/reproducible-research/}
#' @examples
#' Madeyski15SQJ.NDC
#'
"Madeyski15SQJ.NDC"

#' Madeyski15EISEJ.OpenProjects data
#'
#' If you use this data set please cite: Marian Jureczko and Lech Madeyski, "Cross-Project Defect Prediciton: An Empirical Study," (under review), 2015.
#'
#' This paper presents an analysis of 84 versions of industrial, open-source and academic projects. We have empirically evaluated whether those project types constitute separate classes of projects with regard to defect prediction. The predictions obtained from the models trained on the data from the open source projects were compared with the predictions from the other models (built on proprietary, i.e. industrial, student, open source, and not open source projects).
#'
#'
#' @format A data frame with variables:
#' \describe{
#' \item{PROP}{The percentage of classes of proprietary (i.e., industrial) projects that must be tested in order to find 80\% of defects in case of software defect prediction models built on open source projects.}
#' \item{NOTOPEN}{The percentage of classes of projects which are not open source projects that must be tested in order to find 80\% of defects in case of software defect prediction models built on open source projects.}
#' \item{STUD}{The percentage of classes of student (i.e., academic) projects that must be tested in order to find 80\% of defects in case of software defect prediction models built on open source projects.}
#' \item{OPEN}{The percentage of classes of open source projects that must be tested in order to find 80\% of defects in case of software defect prediction models built on open source projects.}
#' }
#'
#' @source \url{http://madeyski.e-informatyka.pl/reproducible-research/}
#' @examples
#' Madeyski15EISEJ.OpenProjects
#'
"Madeyski15EISEJ.OpenProjects"

#' Madeyski15EISEJ.PropProjects data
#'
#' If you use this data set please cite: Marian Jureczko and Lech Madeyski, "Cross-Project Defect Prediciton: An Empirical Study," (under review), 2015.
#'
#' @format A data frame with variables:
#' \describe{
#' \item{NOTPROP}{The percentage of classes of non-proprietary (i.e., non-industrial) projects that must be tested in order to find 80\% of defects in case of software defect prediction models built on proprietary (i.e., industrial) projects.}
#' \item{OPEN}{The percentage of classes of open source projects that must be tested in order to find 80\% of defects in case of software defect prediction models built on proprietary (i.e., industrial) projects.}
#' \item{STUD}{The percentage of classes of student (i.e., academic) projects that must be tested in order to find 80\% of defects in case of software defect prediction models built on proprietary (i.e., industrial) projects.}
#' \item{PROP}{The percentage of classes of proprietary (i.e., industrial) projects that must be tested in order to find 80\% of defects in case of software defect prediction models built on proprietary (i.e., industrial) projects.}
#' }
#'
#' @source \url{http://madeyski.e-informatyka.pl/reproducible-research/}
#' @examples
#' Madeyski15EISEJ.PropProjects
#'
"Madeyski15EISEJ.PropProjects"


#' Madeyski15EISEJ.StudProjects data
#'
#' If you use this data set please cite: Marian Jureczko and Lech Madeyski, "Cross-Project Defect Prediciton: An Empirical Study," (under review), 2015.
#'
#' @format A data frame with variables:
#' \describe{
#' \item{PROP}{The percentage of classes of proprietary (i.e., industrial) projects that must be tested in order to find 80\% of defects in case of software defect prediction models built on student (i.e., academic) projects.}
#' \item{NOTSTUD}{The percentage of classes of projects which are not student projects that must be tested in order to find 80\% of defects in case of software defect prediction models built on student (i.e., academic) projects.}
#' \item{STUD}{The percentage of classes of student (i.e., academic) projects that must be tested in order to find 80\% of defects in case of software defect prediction models built on student (i.e., academic) projects.}
#' \item{OPEN}{The percentage of classes of open source projects that must be tested in order to find 80\% of defects in case of software defect prediction models built on student (i.e., academic) projects.}
#' }
#'
#' @source \url{http://madeyski.e-informatyka.pl/reproducible-research/}
#' @examples
#' Madeyski15EISEJ.StudProjects
#'
"Madeyski15EISEJ.StudProjects"

#' Ciolkowski09ESEM.MetaAnalysis.PBRvsCBRorAR data form a set of primary studies on reading methods for software inspections. They were reported and analysed by M. Ciolkowski ("What do we know about perspective-based reading? an approach for quantitative aggregation in software engineering", in Proceedings of the 3rd International Symposium on Empirical Software Engineering and Measurement, ESEM'09, pp. 133-144, IEEE Computer Society, 2009), corrected and re-analysed by Madeyski and Kitchenham ("How variations in experimental designs impact the construction of comparable effect sizes for meta-analysis", 2015 (to be submitted)).
#'
#' If you use this data set please cite: Lech Madeyski and Barbara Kitchenham, "How variations in experimental designs impact the construction of comparable effect sizes for meta-analysis", 2015 (to be submitted).
#'
#' @format A data frame with 21 rows and 7 variables:
#' \describe{
#' \item{Study}{Name of empirical study}
#' \item{Ref.}{Reference to the paper reporting primary study or experimental run where data were originally reported}
#' \item{Control}{Control treatment: Check-Based Reading (CBR) or Ad-hoc Reading (AR)}
#' \item{Within-subjects}{Yes - if the primary study used the within-subjects experimental design, No - if the primary study did not use the within-subjects experimental design}
#' \item{Cross-over}{Yes - if the primary study used the cross-over experimental design, No - if the primary study did not use the cross-over experimental design}
#' \item{d_ByCiolkowski}{d effect size calculated by Ciolkowski}
#' \item{d_ByOriginalAuthors}{d effect size as reported by the original authors}
#' }
#' @source \url{http://madeyski.e-informatyka.pl/reproducible-research/}
#' @examples
#' Ciolkowski09ESEM.MetaAnalysis.PBRvsCBRorAR
#'
"Ciolkowski09ESEM.MetaAnalysis.PBRvsCBRorAR"


#' MadeyskiKitchenham.MetaAnalysis.PBRvsCBRorAR data form a set of primary studies on reading methods for software inspections. They were analysed by L. Madeyski and B. Kitchenham ("How variations in experimental designs impact the construction of comparable effect sizes for meta-analysis", 2015 (to be submitted)).
#'
#' If you use this data set please cite: Lech Madeyski and Barbara Kitchenham, "How variations in experimental designs impact the construction of comparable effect sizes for meta-analysis", 2015 (to be submitted).
#'
#' @format A data frame with 17 rows and 26 variables:
#' \describe{
#' \item{Study}{Name of empirical study}
#' \item{Ref.}{Reference to the paper reporting primary study or experimental run where data were originally reported}
#' \item{Teams}{The number of teams including both, PBR and Control teams}
#' \item{DesignDesc}{Experimental design description: Before-after, Between-groups, Cross-over}
#' \item{ExpDesign}{Experimental design: between-groups (BG), within-subjects cross-over (WSCO), within-subjects before-after (WSBA)}
#' \item{M_PBR}{The average proportion of defects found by teams using PBR}
#' \item{M_C}{The average proportion of defects found by teams using Control treatment: Check-Based Reading (CBR) or Ad-Hoc Reading (AR)}
#' \item{Diff}{The difference between M_PBR and M_C, i.e. Diff = M_PBR - M_C}
#' \item{Inc}{The percentage increase in defect rate detection, i.e. Inc=100*[(M_PBR-M_C)/M_C]}
#' \item{SD_C_ByAuthors}{The standard deviation of the control group values reproted by the original Authors, i.e., obtained from the papers/raw data}
#' \item{SD_C}{The standard deviation of the control group values equals SD_C_ByAuthors for studies for which the data was available OR the weighted average of SD_C_ByAuthors (i.e., 0.169) for studies where SD_C_ByAuthors is missing.}
#' \item{V_C}{The variance of the Control group observations, i.e., the variance obtained from the teams using the Control method V_C=SD_C^2}
#' \item{V_D}{The variance of the unstandardized mean difference D (between the mean value for the treatment group and the mean value for the Control group)}
#' \item{SD_C_Alt}{This is the equivalent of SD_C (the standard deviation of the control group) based on a different variance for the student studies or the practitioner studies depending on the subject type of the study with the missing value.}
#' \item{V_Alt}{The variance of the mean difference in the meta-analysis based on SD_C_Alt}
#' \item{SS_C}{The sum of squares of the Control group values. For within subjects studies SS=V_C*(n-1). For between subjects studies SS=V_C*(n_C-1)}
#' \item{n_PBR}{The number of PBR teams}
#' \item{n_C}{The number of Control (CBR or AR) teams}
#' \item{ControlType}{Type of Control treatment: CRB or AR}
#' \item{ParticipantsType}{Type of participants: Engineers or Students}
#' \item{TeamType}{Type of team: Nominal or Real}
#' \item{TwoPersonTeamVsLargerTeam}{Reflects size of the teams: 2-PersonTeam or LargerTeam}
#' \item{ArtefactType}{The type of artefact: Requirements or Other}
#' \item{AssociatedWithBasili}{Whether study is associated with Basili (the forerunner): Yes or No}
#' \item{ControlType_Basili}{Combined ControlType and AssociatedWithBasili: AH_AssociatedWithBasili, CBR_AssociatedWithBasili, CBR_NotAssociatedWithBasili}
#' }
#' @source \url{http://madeyski.e-informatyka.pl/reproducible-research/}
#' @examples
#' MadeyskiKitchenham.MetaAnalysis.PBRvsCBRorAR
#'
"MadeyskiKitchenham.MetaAnalysis.PBRvsCBRorAR"

