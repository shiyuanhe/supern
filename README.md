# supern
Type Ia Supernova Functional Template Model

## Folder Strucutre

1. In the ```data``` folder is a SNIa data table. The data table contains all SNIa samples in this paper with their scores, spectral line strength and spectral class.
2. In the ```supern``` folder is an R package. The software provides a web-based user interface for light curve fitting,  intrinsic color estimation, spctral line strength estimation and spectral classes determination. It also provide the probability that the submitted sample belongs to SNIa. 

## The supern Package

Install the R package with 

```r
library(devtools)
install_github("shiyuanhe/supern/supern")
```

Then start the web user interface on your computer by

```r
library(supern)
runWebInterface()
```

The web interface is also hosted at [http://heshiyuan.name:3838/snelc/](http://heshiyuan.name:3838/snelc/).
