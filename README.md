
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![minimal R
version](https://img.shields.io/badge/R%3E%3D-4.2.1-brightgreen.svg)](https://cran.r-project.org/)
[![Licence](https://img.shields.io/github/license/mashape/apistatus.svg)](http://choosealicense.com/licenses/mit/)

## Research compendium for alternative COVID-19 mortality analyses

This is a working R compendium (think R package but for reproducible
analysis). The analysis directory contains R scripts used to generate
the results.

### Installation

    git clone https://github.com/mrc-ide/covid-alternative-mortality.git
    cd covid-alternative-mortality
    open covid-alternative-mortality.Rproj
    devtools::install_deps()

### Overview

The structure within analysis is as follows:

    analysis/
        |
        ├── 01_xxxxx /           # analysis scripts used for generating figures
        |
        ├── figures/              # location of figures produced by the analysis scripts
        |
        ├── data/
        │   ├── DO-NOT-EDIT-ANY-FILES-IN-HERE-BY-HAND
        │   ├── raw_data/       # data obtained from elsewhere
        │   └── derived_data/   # data generated during the analysis

### Compendium DOI:

<https://doi.org/10.5281/zenodo.7185171>

The files at the URL above will generate the results as found in the
publication.

### The R package

This repository is organized as an R package. There are no/negligable R
functions exported in this package - the majority of the R code is in
the analysis. The R package structure is here to help manage
dependencies, to take advantage of continuous integration, and so we can
keep file and data management simple.

To download the package source as you see it on GitHub, for offline
browsing, use this line at the shell prompt (assuming you have Git
installed on your computer):

``` r
git clone https://github.com/mrc-ide/covid-alternative-mortality.git
```

Once the download is complete, open the
`covid-alternative-mortality.Rproj` in RStudio to begin working with the
package and compendium files. We will endeavour to keep all package
dependencies required listed in the DESCRIPTION.

In addition we use `renv` to track package dependencies for
reproducibility. Please use `renv::restore` to restore the state of the
project and see <https://rstudio.github.io/renv/articles/renv.html> for
more information.

### Licenses

Code: [MIT](http://opensource.org/licenses/MIT) year: 2022, copyright
holder: OJ Watson

Data: [CC-0](http://creativecommons.org/publicdomain/zero/1.0/)
attribution requested in reuse
