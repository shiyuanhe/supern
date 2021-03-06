# supern
Type Ia Supernova Functional Template Model

## Folder Strucutre

1. The ```data``` folder contains a sample light curve data and a SNIa data table. The data table contains all SNIa samples in this paper with the scores, spectral class and etc.
2. In the ```supern``` folder is an R package. The software provides a web-based user interface for light curve fitting,  intrinsic color estimation, and spectral classes determination. It also provide the probability that the submitted sample belongs to SNIa. This implements the band-vague model.
3. In the ```template``` folder are the SNe template for both the band-vauge and band-specific model. The mean function and the first four functional principal components are provided.

## The supern Package
The supern Package provides a web user interface for light curve fitting. Light curve file can be uploaded for analysis.


To install the package, we need to install the prerequisite packages firstly,

```r
install.packages("R6", dependencies= TRUE)
install.packages("shiny", dependencies= TRUE)
install.packages("quadprog", dependencies= TRUE)
install.packages("devtools", dependencies= TRUE)
```


Then install the ```supern``` package from github by

```r
library(devtools)
install_github("shiyuanhe/supern/supern")
```

At last, start the web user interface on your computer by


```r
library(supern)
runWebInterface()
```


### Submitted file format

1. Provide SNe name, survey, redshift in file head.
2. The light curve data begins with '#Data'.
3. The light curve data is space (or tab) seperated with columns: filter, MJD, Mag, Mag err.
4. The light curves should have standard filters: B, V, R, I.
5. An example file ```SN2008fp_CSP_main.txt``` is in the ```data``` folder.

```
#Name sn2008fp
#Redshift 0.02216
#Survey CSP
#Data
V 54724.323 14.249 0.006
V 54725.320 14.145 0.006
...
```