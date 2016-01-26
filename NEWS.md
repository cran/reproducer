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
- minor visual changes in figures produced by the following methods:
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
