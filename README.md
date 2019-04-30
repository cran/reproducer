# reproducer

The R package **reproducer** is aimed to support reproducible research in software engineering. See the package [homepage](http://madeyski.e-informatyka.pl/reproducible-research/) for details and examples.

## Installation

One may install the stable version from
[CRAN](https://cran.r-project.org/package=reproducer):

```r
install.packages('reproducer', dependencies = TRUE)
```

You can use **devtools** to install the development version from my web site:

```r
install.packages("devtools", dependencies = T, repos = "http://cran.r-project.org/")
library(devtools)
devtools::install_url("http://madeyski.e-informatyka.pl/download/R/reproducer_0.3.0.tar.gz")
library(reproducer)
```

## Motivation
The motivation is to support reproducible research in software engineering via sharing data sets and code behind the published or just submitted papers.
