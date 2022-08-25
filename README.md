# gASCA

## Description

Implementation of a GLM based Analysis of variance – Simultaneous Component Analysis (ASCA or ANOVA–SCA), particularly focused on the analysis of count data.

## Installation

```{r}
## if devtools package is not installed
install.packages("devtools")

## install the last version from github
devtools::install_github("pietrofranceschi/gASCA", build_vignettes=TRUE)

```

## Details

ASCA is a multivariate method which can be used to analyze designed experiments and highlight the variables which are more affecting each specific design factor ([Wikipedia](https://en.wikipedia.org/wiki/ANOVA%E2%80%93simultaneous_component_analysis)).

In the ASCA package implementation the expected values for each factor are estimated trough Generalized Linear Modeling (GLM) this extend the original approach to unbalanced designs and to non normally distributed data like counts.


## References

-  Jansen, J. J., Hoefsloot, H. C., van der Greef, J., Timmerman, M. E., Westerhuis, J. A., & Smilde, A. K.. "ASCA: analysis of multivariate data obtained from an experimental design". *Journal of Chemometrics: A Journal of the Chemometrics Society*, 19(9) (2005), 469-481.

-   Thiel, Michel, Baptiste Feraud, and Bernadette Govaerts. "ASCA+ and APCA+: Extensions of ASCA and APCA in the analysis of unbalanced multifactorial designs." *Journal of Chemometrics* 31.6 (2017): e2895.

-   Martin, Manon, and Bernadette Govaerts. "LiMM‐PCA: Combining ASCA+ and linear mixed models to analyse high‐dimensional designed data." *Journal of Chemometrics* 34.6 (2020): e3232.
