## mixtools: Tools for Analyzing Finite Mixture Models	<a href='https://github.com/dsy109/mixtools'><img src='man/figures/mixtools.png' align="right" height="138.5" /></a>

[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)	
![CRAN/METACRAN](https://img.shields.io/cran/l/mixtools)	
[![CRAN status](https://www.r-pkg.org/badges/version/mixtools)](https://CRAN.R-project.org/package=mixtools)
[![Dependencies](https://tinyverse.netlify.com/badge/mixtools)](https://cran.r-project.org/package=mixtools)
![GitHub last commit](https://img.shields.io/github/last-commit/dsy109/mixtools)
[![Downloads](https://cranlogs.r-pkg.org/badges/mixtools?color=brightgreen)](https://www.r-pkg.org/pkg/mixtools)
[![JSS](https://img.shields.io/badge/JSS-10.18637%2Fjss.v032.i06-brightgreen)](http://dx.doi.org/10.18637/jss.v032.i06)
[![HAL](https://img.shields.io/badge/hal-00717545-brightgreen)](https://hal.archives-ouvertes.fr/hal-00717545/document)

### Synopsis

The [mixtools package](https://cran.r-project.org/package=mixtools) provides functions for analyzing finite mixture models.  Parametric and semiparametric mixture models are handled.  Various tools are include for determining the number of components.   Visualizations include histograms with the estimated mixture model overlaid, 2D and 3D scatterplots for relevant mixture fits, and the [mixturegram](https://doi.org/10.1080/10618600.2017.1398093).    More details about the package are included in both the original [JSS](http://dx.doi.org/10.18637/jss.v032.i06) article as well as a technical report on the open archive [HAL](https://hal.archives-ouvertes.fr/hal-00717545/document).

Other highlights:

- Includes estimation of various parametric and semiparametric mixtures-of-regressions models.

- Functions to help with determining the number of components, including bootstrapping the likelihood ratio test statistic, mixturegrams, and model selection criteria.

- Functions available for estimating Reliability Mixture Models (RMMs).

- Includes a Metropolis-Hastings algorithm for estimating a mixture-of-linear-regressions model.

- Includes some real datasets for which mixture models have been shown to provide good fits.

### Documentation

The [JSS](http://dx.doi.org/10.18637/jss.v032.i06) article and the technical report on [HAL](https://hal.archives-ouvertes.fr/hal-00717545/document) both provide documentation about the mixtools package.  Moreover, the [help file](https://CRAN.R-project.org/package=mixtools) also documents the references used for each function.

### Examples

Additional examples for the mixtools package are currently being developed for a Shiny app.

### Installation

Released and tested versions of mixtools are available via the
[CRAN](https://cran.r-project.org) network, and can be installed from within R via

```R
install.packages("mixtools")
```

### Support

The [issue tickets at the GitHub repo](https://github.com/dsy109/mixtools/issues)
are the primary bug reporting interface.  As with the other web resources,
previous issues can be searched as well.

### Authors

Derek S. Young, Tatiana Benaglia, Didier Chauveau, David Hunter, Kedai Cheng, Ryan Elmore, Thomas Hettmansperger, Hoben Thomas, Fengjuan Xua

### License

GPL (>= 2)

### Funding Acknowledgment

This package is based upon work supported by the National Science Foundation, Grant Number [SES-0518772](https://www.nsf.gov/awardsearch/showAward?AWD_ID=0518772) and the Chan Zuckerberg Initiative, Grant Number [2020-225193](https://chanzuckerberg.com/eoss/proposals/enhancing-usability-of-mixtools-and-tolerance-for-the-biomedical-community/).

### Code of Conduct

As contributors and maintainers of this project, we pledge to respect all people who 
contribute through reporting issues, posting feature requests, updating documentation, 
submitting pull requests or patches, and other activities.  Both contributors and 
maintainers must consistently demonstrate acceptable behavior, respectful communications, 
and professional conduct.  Project maintainers have the right and responsibility to remove, 
edit, or reject comments, commits, code, wiki edits, issues, and other contributions that 
are not aligned to this Code of Conduct.  Project maintainers who do not follow the 
Code of Conduct may be removed from the project team.  Instances of abusive, harassing, 
or otherwise unacceptable behavior may be reported by opening an issue or contacting one 
or more of the project maintainers.  By contributing to this project, you agree to abide 
by its terms.

We are here for a love of coding and a passion for cultivating knowledge.  Let us enjoy 
this collaboration together!


