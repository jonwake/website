# SAE Short Course

*Course summary for SAE Conference Short Course, Bucharest June 19, 2026*

Small area estimation is of crucial importance in low- and middle-income countries (LMICs). A modern Bayesian treatment will be presented and illustrated using a range of examples. Area-level (Fay–Herriot) and unit-level models will be presented. Unit-level models for both linear and generalized linear models will be discussed. Fast computation is carried out with the Integrated Nested Laplace Approximation (INLA) method, which is embedded within the SUMMER and surveyPrev R packages. Hyperprior specification is via penalized complexity priors. Between-area variation will be modeled using independent and spatial random effects. For the latter, the Besag, York, Mollié model will be described.

*Short Course Slides*

- Lecture 1: Context and Motivation: [Slides](files/SAE-files/SAE-Motivation.pdf)
- Lecture 2: Introduction to Bayes: [Slides](files/SAE-files/SAE-Bayes.pdf)
- Lecture 3: Area-Level Models: [Slides](files/SAE-files/SAE-Area-Level-Models.pdf)
- Lecture 4: Unit-Level Models: [Slides](files/SAE-files/SAE-Unit-Level-Models.pdf)
- Lecture 5: Further Topics: [Slides](files/SAE-files/SAE-FurtherTopics.pdf)

*R Packages:*

[surveyPrev](https://cran.r-project.org/web/packages/surveyPrev/index.html) with [vignette on prevalence mapping](https://cran.r-project.org/web/packages/surveyPrev/vignettes/vignette-main.html) and [vignette on creating indicators](https://cran.r-project.org/web/packages/surveyPrev/vignettes/vignette-data-preparation.html)

[SUMMER](https://cran.r-project.org/web/packages/SUMMER/index.html) See the cran site for various vignettes 

*DHS Data:*

Demographic and Heath Surveys (DHS) data can be downloaded, after registering for an account [here](https://dhsprogram.com/) 

When requesting specific datasets, remember to request the GPS data (locations of clusters)

*Web Apps*

[DHS R Shiny App](https://rsc.stat.washington.edu/sae4health/)

[MICS R Shiny App (beta test version)](https://rsc.stat.washington.edu/app_direct/sae4mics/) 

[Pre-modeled App](https://sae4lmic.stat.uw.edu/) 

*Papers:*

[Wakefield, Fuglstad, Riebler, Godwin, Wilson, Clark (2019). Estimating under-five mortality in space
and time in a developing world context. SMMR, 28, 2614-2634](files/SAE-files/wakefieldetal19.pdf)

[Wakefield, Ononek, Pedersen (2020). Small Area Estimation for disease prevalence mapping. International Statistical Review, 88, 398-418](files/SAE-files/wakefieldokonekpedersen20.pdf)

[Dong, Wu, Li, Wakefield (2026). Toward a principled workflow for prevalence mapping using household survey data. Journal of Survey Statistics and Methodology, 14, 209–237](files/SAE-files/DongWorkflow2026.pdf)

[Vignette for Dong, Wu. Li, Wakefield (2026). Toward a principled workflow for prevalence mapping using household survey data. Journal of Survey Statistics and Methodology, 14, 209–237](files/SAE-files/DongWorkflow2026Vignette.pdf)

[Wakefield, Gao, Fuglstad, Li (2026). The two cultures for prevalence mapping: small area estimation and model-based geostatistics, To appear (with discussion). Statistical Science](files/SAE-files/Two_Cultures-Main.pdf)

[Additional Materials at the sae4health website](https://sae4health.stat.uw.edu/)

