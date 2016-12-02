# supern
Type Ia Supernova Functional Template Model

## Folder Strucutre

1. In the ```data``` folder contains a sample light curve data and a SNIa data table. The data table contains all SNIa samples in this paper with the scores, spectral class and etc.
2. In the ```supern``` folder is an R package. The software provides a web-based user interface for light curve fitting,  intrinsic color estimation, and spectral classes determination. It also provide the probability that the submitted sample belongs to SNIa. 

## The supern Package
The supern Package provides a web user interface for light curve fitting. Light curve file can be uploaded for analysis.



To install the package, we need to install the prerequisite packages firstly,

```r
install.packages("R6")
install.packages("shiny")
install.packages("quadprog")
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

The web interface is also hosted at [http://heshiyuan.name:3838/snelc/](http://heshiyuan.name:3838/snelc/).

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