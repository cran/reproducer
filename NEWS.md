# CHANGES IN reproducer VERSION 0.4.0

## NEW FEATURES
- Data set: 
`KitchenhamEtAl.CorrelationsAmongParticipants.Madeyski10`,
`KitchenhamEtAl.CorrelationsAmongParticipants.Scanniello17TOSEM`,
`KitchenhamEtAl.CorrelationsAmongParticipants.Ricca10TSE`,
`KitchenhamEtAl.CorrelationsAmongParticipants.Romano18ESEM`,
`KitchenhamEtAl.CorrelationsAmongParticipants.Scanniello14JVLC`,
`KitchenhamEtAl.CorrelationsAmongParticipants.Reggio15SSM`,
`KitchenhamEtAl.CorrelationsAmongParticipants.Gravino15JVLC`,
`KitchenhamEtAl.CorrelationsAmongParticipants.Ricca14TOSEM`,
`KitchenhamEtAl.CorrelationsAmongParticipants.Scanniello14EASE`,
`KitchenhamEtAl.CorrelationsAmongParticipants.Abrahao13TSE`,
`KitchenhamEtAl.CorrelationsAmongParticipants.Torchiano17JVLC`,
`KitchenhamEtAl.CorrelationsAmongParticipants.Scanniello15EMSE`,
`KitchenhamEtAl.CorrelationsAmongParticipants.Scanniello14TOSEM`,


- Functions including computational procedures used to reproduce the main findings in a joint paper
(planned to be submitted): Barbara Kitchenham, Lech Madeyski, Giuseppe Scanniello and 
Carmine Gravino, "The importance of the Correlation between Results from Individual Participants in 
Crossover Experiments": 
`CalculateRLevel1`,
`ExtractGroupSizeData`,
`ConstructLevel1ExperimentRData`,
`ExtractExperimentData`,
`CalculateLevel2ExperimentRData`,
`ExtractSummaryStatisticsRandomizedExp`,
`calculateBasicStatistics`,
`calculateGroupSummaryStatistics`,
`rSimulations`


# CHANGES IN reproducer VERSION 0.3.1

## NEW FEATURES
- Data set: `MadeyskiLewowski.IndustryRelevantGitHubJavaProjects20191022` - over 15\% of entries present in this data set is not present in the previous data set `MadeyskiLewowski.IndustryRelevantGitHubJavaProjects20190324` due to moved time windows for the project creation and last push dates.


## UPDATED FEATURES
- Updated function `searchForIndustryRelevantGitHubProjects` - now supports flexible creation date and last push thresholds (enabling the script to better support researchers interested in gathering evolving data sets).


# CHANGES IN reproducer VERSION 0.3.0

## NEW FEATURES
- Functions including computational procedures used to reproduce the main findings, e.g., tables, in a joint paper: Barbara Kitchenham, Lech Madeyski, Pearl Brereton, "Meta-analysis for Families of Experiments in Software Engineering: A Systematic Review and Reproducibility and Validity Assessment": 
`transformHgtoZr`,
- Function to search for industry relevant software projects available from GitHub related to a report: Lech Madeyski, “Training data preparation method,” tech. rep., code quest (research project NCBiR POIR.01.01.01-00-0792/16), 2019; as well as a paper: Tomasz Lewowski and Lech Madeyski, "Creating evolving project data sets in software engineering", 2019:
`searchForIndustryRelevantGitHubProjects`
- Data set: `MadeyskiLewowski.IndustryRelevantGitHubJavaProjects20190324`


## UPDATED FEATURES
- Updated function `reproduceTablesOfPaperMetaAnalysisForFamiliesOfExperiments`

# CHANGES IN reproducer VERSION 0.2.1

## UPDATED FEATURES
- Fixed version of the package in README
- Added warning to the `ExtractMAStatistics` function: it works with `metafor` version 2.0-0, but changes to metafor's method of providing access to its individual results may introduce errors into the function. 


# CHANGES IN reproducer VERSION 0.2.0

## NEW FEATURES
- functions including computational procedures used to reproduce the main findings, e.g., tables, in a joint paper: Barbara Kitchenham, Lech Madeyski, Pearl Brereton, "Meta-analysis for Families of Experiments: A Systematic Review and Reproducibility Assessment": 
`calculateSmallSampleSizeAdjustment`, 
`constructEffectSizes`, 
`transformRtoZr`, 
`transformZrtoR`, 
`transformHgtoR`, 
`calculateHg`, 
`transformRtoHg`, 
`transformZrtoHgapprox`, 
`transformZrtoHg`, 
`PrepareForMetaAnalysisGtoR`, 
`ExtractMAStatistics`, 
`aggregateIndividualDocumentStatistics`, 
`reproduceTablesOfPaperMetaAnalysisForFamiliesOfExperiments`.
- unit tests
- data sets used in the aforementioned paper:  `KitchenhamMadeyskiBrereton.MetaAnalysisReportedResults`,  `KitchenhamMadeyskiBrereton.ABBAMetaAnalysisReportedResults`,  `KitchenhamMadeyskiBrereton.ReportedEffectSizes`,  `KitchenhamMadeyskiBrereton.ABBAReportedEffectSizes` `KitchenhamMadeyskiBrereton.ExpData`, and `KitchenhamMadeyskiBrereton.DocData`

## UPDATED FEATURES


# CHANGES IN reproducer VERSION 0.1.9

## NEW FEATURES
- added unit tests and assertions (run-time tests)

## UPDATED FEATURES
- updated data set `MadeyskiKitchenham.EUBASdata` and functions `getEffectSizesABBA`, `effectSizeCI`



# CHANGES IN reproducer VERSION 0.1.8

## NEW FEATURES
- added function `getTheoreticalEffectSizeVariancesABBA`

## UPDATED FEATURES
- updated functions `getSimulationData`, `plotOutcomesForIndividualsInEachSequenceGroup`, `getEffectSizesABBA`, `effectSizeCI`



# CHANGES IN reproducer VERSION 0.1.7

## NEW FEATURES
- added function `effectSizeCI` to calculate 95% Confidence Intervals (CI) on Standardised Effect Sizes (d) for cross-over repeated-measures designs

## REMOVED FEATURES
- removed field dIG.Var.CurtinAve from the data frame returned by the `reproduceSimulationResultsBasedOn500Reps1000Obs` function (we agreed to write joint paper with Dr Curtin describing corrections to his equations to calculate effect size variances for continuous outcomes of cross-over clinical trials)


# CHANGES IN reproducer VERSION 0.1.6

## NEW FEATURES
- added functions to reproduce simulation data, to calculate effect sizes and their variance discussed in a paper "Effect Sizes and their Variance for AB/BA Crossover Design Studies" by Lech Madeyski and Barbara Kitchenham.
    - `getSimulationData`
    - `plotOutcomesForIndividualsInEachSequenceGroup`
    - `getEffectSizesABBA`
    - `getEffectSizesABBAIgnoringPeriodEffect`
    - `reproduceSimulationResultsBasedOn500Reps1000Obs`
    - `percentageInaccuracyOfLargeSampleVarianceApproximation`
    - `proportionOfSignificantTValuesUsingCorrectAnalysis`
    - `proportionOfSignificantTValuesUsingIncorrectAnalysis`
- new data set:
    - `KitchenhamMadeyski.SimulatedCrossoverDataSets` backed by functions (`varianceSimulation`, `getSimulatedCrossoverDataSets`) to reproduce the data set.

## REMOVED FEATURES
- removed function:
    - `cloudOfWords`


# CHANGES IN reproducer VERSION 0.1.5
## BUG FIXES
- corrected typo in a surname of one of the data contributors
- corrected the name of the university in Hong Kong
- corrected the boxplotAndDensityCurveOnHistogram function
- updated references to research papers supported by the package


# CHANGES IN reproducer VERSION 0.1.4

## NEW FEATURES
- new data sets:
    - `KitchenhamMadeyskiBudgen16.FINNISH`
    - `KitchenhamMadeyskiBudgen16.PolishSubjects`
    - `KitchenhamMadeyskiBudgen16.SubjectData`
    - `KitchenhamMadeyskiBudgen16.PolishData`
    - `KitchenhamMadeyskiBudgen16.DiffInDiffData`
    - `KitchenhamMadeyskiBudgen16.COCOMO`
- updated to be compatible with the ggplot2 2.0.0
- minor visual changes in figures produced by the following functions:
    - `densityCurveOnHistogram`
    - `boxplotHV`
    - `boxplotAndDensityCurveOnHistogram`


# CHANGES IN reproducer VERSION 0.1.3

## NEW FEATURES
- added general functions
    - `printXTable`
    - `cloudOfWords`
- added functions to reproduce all of the tables, figures and outputs in the paper by Lech Madeyski and Barbara Kitchenham on how variations in experimental designs impact the construction of comparable effect sizes for meta-analysis:
    - `reproduceForestPlotRandomEffects`
    - `reproduceMixedEffectsAnalysisWithEstimatedVarianceAndExperimentalDesignModerator`
    - `reproduceMixedEffectsAnalysisWithExperimentalDesignModerator`
    - `reproduceMixedEffectsForestPlotWithExperimentalDesignModerator`
    - `reproduceTableWithEffectSizesBasedOnMeanDifferences`
    - `reproduceTableWithPossibleModeratingFactors`
    - `reproduceTableWithSourceDataByCiolkowski`
- added data sets:
    - `Ciolkowski09ESEM.MetaAnalysis.PBRvsCBRorAR`
    - `MadeyskiKitchenham.MetaAnalysis.PBRvsCBRorAR`

## BUG FIXES
- corrected `Madeyski15EISEJ.StudProjects$STUD` data set

# CHANGES IN reproducer VERSION 0.1.2

## NEW FEATURES
- The first published version of reproducer. It includes data sets: 
    - `Madeyski15SQJ.NDC`
    - `Madeyski15EISEJ.OpenProjects`
    - `Madeyski15EISEJ.PropProjects`
    - `Madeyski15EISEJ.StudProjects`
and functions (for importing data, visualization and descriptive analyses):
    - `readExcelSheet`
    - `densityCurveOnHistogram`
    - `boxplotHV`
    - `boxplotAndDensityCurveOnHistogram`

See the package homepage (https://madeyski.e-informatyka.pl/reproducible-research/) for documentation and examples.
