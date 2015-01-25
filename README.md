# reproducer

The R package **reproducer** is aimed to support reproducible research in software engineering. See the package [homepage](http://madeyski.e-informatyka.pl/reproducible-research/) for details and examples.

## Installation

In a close future you will be able to install the stable version on
[CRAN](http://cran.rstudio.com/package=reproducer):

```r
install.packages('reproducer', dependencies = TRUE)
```

You can use **devtools** to install the development version from my web site:

```r
install.packages("devtools", dependencies = T, repos = "http://cran.rstudio.com/")
library(devtools)
devtools::install_url("http://madeyski.e-informatyka.pl/download/R/reproducer_0.1.2.tar.gz")
library(reproducer)
```

## Motivation
The motivation is to support reproducible research in software engineering via sharing data sets and code behind the published or just submitted papers.
