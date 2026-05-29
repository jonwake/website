# Software

## Small Area Estimation

- SAE for prevalence software summaries at [SAE4Health](https://wu-thomas.github.io/surveyPrev_website/). The next two products are reachable from here.
- The R package [surveyPrev](https://cran.r-project.org/web/packages/surveyPrev/index.html) carries out small area estimation using DHS data and a variety of statistically rigorous methods. See also [SAE4Health](https://wu-thomas.github.io/surveyPrev_website/).
- The associated [shinyApp](https://rsc.stat.washington.edu/surveyPrevRShiny/) provides a user-friendly interface to the package. See also [SAE4Health](https://wu-thomas.github.io/surveyPrev_website/).
- Methods for smoothing of direct estimates in time and space, and mapping the subsequent estimates, are contained in the SUMMER package at CRAN: [SUMMER](https://cran.r-project.org/web/packages/SUMMER/index.html).
- GitHub version of SUMMER (contains the latest version): [GitHub SUMMER](https://github.com/bryandmartin/SUMMER).
- Shiny App for SUMMER: [SUMMER Shiny App](https://www.dropbox.com/sh/s5i9bn0kzk59nm3/AAAUcR0GHcKeE4FB8tt-cr3na?dl=0).
- ShinyApp for TMB/INLA comparison: [Shiny App for TMB/INLA Simulations](https://rsc.stat.washington.edu/content/34/). TMB/INLA paper: [here](https://arxiv.org/abs/2103.09929).
- TMB Code: [here](https://github.com/aaron-oz/tmb-inla-paper/). TMB/INLA paper: [here](https://arxiv.org/abs/2103.09929).
- Files to accompany Wakefield, Fuglstad, Riebler, Godwin, Wilson and Clark, "Estimating Under 5 Mortality in Space and Time in a Developing World Context." Paper: here. Supplementary Materials: here.
- Reference: Chen, Wakefield and Lumley, "The use of sample weights in Bayesian spatial hierarchical models for small area estimation" is [here](files/software/SpatialSurveyExample.zip).
- Code to implement the simulation in Wakefield, Simpson and Godwin is [here](files/software/intoSpaceCode.R).

## Spatial Epidemiology

- Details of the R package SpatialEpi to perform various methods in Spatial Epidemiology: [SpatialEpi](http://cran.r-project.org/web/packages/SpatialEpi/index.html).
- WinBUGS code for models fitted in Bauer and Wakefield, "Stratified space-time infectious disease modeling: with application to hand, foot and mouth disease in China": [here](files/software/EpidemicWinBUGScode.txt).

## Genetic Epidemiology

- Excel spreadsheet for calculating approximate Bayes factors for genetic epidemiology studies: [here](files/BFDP.xls). R code: [here](files/BFDP.R) and [here](files/BFDP2.R) (the latter for combining two studies also). Example: [here](files/BayesAssocExample.R). The paper "A Bayesian Measure of the Probability of False Discovery in Genetic Epidemiology Studies" is on the publications page.
- Details of the R package to evaluate methods described in the [Biometrics paper](files/HWbiometrics.pdf): [HWEBAYES](http://cran.r-project.org/web/packages/HWEBayes/index.html). WinBUGS code for examining Hardy-Weinberg Equilibrium: [here](files/HWwinbugs.txt). R script to run the four-allele example from the *Biometrics* paper: [here](files/HWEBayesFourAllele.R). Short user manual: [here](files/HWEBayes-UserGuide.pdf).

## Two-Phase Sampling

- Bayes Analysis of Two-Phase Data: Random Sampling Code is [here](files/BayesTwoPhase.R) and example is [here](files/BayesTwoPhase_Example.R). Bayes Analysis of Two-Phase Data: Case-Control Sampling Code is [here](files/BayesTwoPhase_CaseControl.R) and example is [here](files/BayesTwoPhase_CaseControl_Example.R). Accompanying paper appears in *Biometrics*. Code was written by Michelle Ross.
- The R Code to carry out two-phase sampling in an ecological setting as described in "Overcoming ecological bias using the two-phase study design" by Wakefield and Haneuse, *American Journal of Epidemiology* (see the publications page) is [here](files/TPSanalysis.q), and requires the R functions [tps.q](files/tps.q) — the latter is an R version of the Splus code written by Norm Breslow and Nilanjan Chatterjee. The code works with the North Carolina data [infantsAgg.dat](files/infantsAgg.dat).

## Longitudinal Data

- Bayesian population PK/PD modeling developed for BUGS by Jon Wakefield, Dave Lunn, Nicky Best and Dave Spiegelhalter: [PKBUGS](https://www.mrc-bsu.cam.ac.uk/software/bugs/the-bugs-project-pkbugs/).
- BUGS code for the Zhou and Wakefield (2006, *Biometrics*) paper on curve clustering: [here](files/biometrics_code.txt).
