## Level 0, a collection of SNe

heapObsMag_asList <- function(sneHeap, basisObj, weighting, meanC = NULL){
    tmp = lapply(sneHeap, sneObsMatStack_noMax, ##XXX
                 basisObj = basisObj, 
                 weighting = weighting,
                 meanC = meanC)
    tmp = unlist(tmp, recursive = FALSE)
    tmp = lapply(tmp, function(tmat) tmat[,"mag"])
    return(tmp)
}


heapObsM_asMat <- function(sneHeap, basisObj, weighting, meanC = NULL){
    tmp = lapply(sneHeap, sneObsMatStack_noMax, 
                 basisObj = basisObj, 
                 weighting = weighting,
                 meanC = meanC)
    tmp2 = unlist(tmp, recursive = FALSE)
    tmp2 = do.call(rbind, tmp2)
    return(tmp2)
}



heapBasisM_asList_weighted <- function(sneHeap, basisObj){
    tmp = lapply(sneHeap, sneBasisMStack_weighted, basisObj = basisObj)
    tmp = unlist(tmp, recursive = FALSE)
    return(tmp)
}


heapBasisM_asMat_weighted <- function(sneHeap, basisObj){
    tmp = lapply(sneHeap, sneBasisMStack_weighted, basisObj = basisObj)
    tmp2 = unlist(tmp, recursive = FALSE)
    tmp2 = do.call(rbind, tmp2)
    return(tmp2)
}

## Level 1, sne level
sneObsMatStack_noMax <- function(cSNe, basisObj, 
                                 weighting, meanC = NULL){
    tmp = lapply(cSNe$LCurves, lcObsMat_noMax, 
                 basisObj = basisObj,
                 weighting = weighting,
                 meanC = meanC)
    return(tmp)
}


sneBasisMStack_weighted <- function(cSNe, basisObj){
    tmp = lapply(cSNe$LCurves, lcBasisMat, 
                 basisObj = basisObj)
    return(tmp)
}


## Level 2, light curve level
lcObsMat_noMax <- function(cLC, basisObj, 
                           weighting = FALSE,
                           meanC = NULL){
    obsMat = cLC$getObsMatInPhase()
    obsMat = obsMatDeMax(obsMat, cLC)
    obsMat = obsMatClean(obsMat, basisObj)
    if(!is.null(meanC))
        obsMat = obsMatDeMeanC(obsMat, meanC)
    if(weighting)
        obsMat = obsMatAddWeight(obsMat)
    return(obsMat)
}


lcBasisMat <- function(cLC, basisObj){
    obsMat = lcObsMat_noMax(cLC, basisObj, TRUE)
    Basis = basisObj$bD0EvalAt(obsMat[ ,"phase"])
    Basis = basisAddWeight(Basis, obsMat[ ,"sigma"] + 0.01)
    return(Basis)
}



## Level 3, basic functions
## remove mmax
obsMatDeMax = function(obsMat, LCObject){
    obsMat[ ,"mag"] = obsMat[,"mag"] - LCObject$params$mmax
    return(obsMat)
}

obsMatDeMeanC <- function(obsMat, meanC){
    obsMat[ ,"mag"] = obsMat[,"mag"] - 
        meanC(obsMat[,"phase"])
    return(obsMat)
}

## Cut into phase range
## Remove outlier
obsMatClean = function(obsMat, basisObject){
    sel = obsMat[ ,"phase"] > basisObject$phase_min &
        obsMat[ ,"phase"] < basisObject$phase_max &
        obsMat[ ,"mag"] < 10 & 
        obsMat[ ,"mag"] > (-3)
    obsMat = subset(obsMat, sel)
    return(obsMat)
}

obsMatAddWeight = function(obsMat){
    w = 1/(obsMat[, "sigma"] + 0.01)
    obsMat[, "mag"] = obsMat[, "mag"] * w
    return(obsMat)
}

basisAddWeight = function(basisMat, sigma){
    basisMat = sweep(basisMat, MARGIN = 1, ## stat corresponds to row of basisMat
                     STATS = 1/sigma, FUN = '*')
    return(basisMat)
}
