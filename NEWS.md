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

See the package homepage (http://madeyski.e-informatyka.pl/reproducible-research/) for documentation and examples.
