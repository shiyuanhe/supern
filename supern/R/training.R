curveListReduce = function(allCurves){
    n = length(allCurves)
    mCurve = Reduce(f = "+", allCurves) / n
    allCurvesP2 = lapply(allCurves, function(tmp) tmp^2)
    mCurveS2 = Reduce(f = "+", allCurvesP2) / n
    mCurveSD = sqrt(mCurveS2 - mCurve^2)
    return(list(mCurve = mCurve, mCurveSD = mCurveSD))
}



modelTrainingWrap = function(heapList, phaseRange){
    basisObj = basisClass$new(phaseRange[1], phaseRange[2])
    modelObj = sneLCModel$new()
    cvres = modelObj$trainingMean(heapList, basisObj, 
                                  lambda = 100, cv = TRUE)
    # cvres / 1e8
    # modelObj$plotMeanC(heapList)
    modelObj$trainingFPCA(heapList, basisObj, 
                          lambda1 = 1e-5, lambda3 = 1e-3, 
                          tK = 7)
    #modelObj$plotTempC(4)
    return(modelObj)
}


cvFitting = function(heapList, figdir, phaseRange, bootstrap = FALSE){
    rmList <<- numeric(0)
    nbSNe = length(heapList)
    
    pp_grid = seq(phaseRange[1], phaseRange[2], length.out = 500)
    mCurves = phiCurves = list()
    
    for(sneI in 1:nbSNe){
        print(sneI)
        tryCatch({
            modelObj = modelTrainingWrap(heapList[-sneI], phaseRange)
            cSNe = heapList[[sneI]]
            res = sapply(cSNe$LCurves, FUN = modelObj$modelFit)
            if(bootstrap)
                res2 = lapply(cSNe$LCurves, FUN = modelObj$bootFit, nBoot = 50)
            if(!all(res))
                stop("max failed!")
            cSNe$plotFit(figdir)
            
            mm_grid = modelObj$mean_fun(pp_grid)
            phis_grid = modelObj$template_fun(pp_grid)
            mCurves = c(mCurves, list(mm_grid))
            phiCurves = c(phiCurves, list(phis_grid))
            
            
        }, error = function(e){
            rmList <<- c(rmList, sneI)
        })
    }
    
    mCurve = curveListReduce(mCurves)
    phiCurve = curveListReduce(phiCurves)
    
    if(length(rmList) > 0)
        fittedList = heapList[-rmList]
    else
        fittedList = heapList
    
    return(list(fittedList = fittedList, 
                mCurve = mCurve, 
                phiCurve = phiCurve,
                pp_grid = pp_grid))
}



allFitting = function(heapList, figdir, phaseRange, bootstrap = FALSE){
    rmList <<- numeric(0)
    nbSNe = length(heapList)
    
    modelObj = modelTrainingWrap(heapList, phaseRange)
    
    for(sneI in 1:nbSNe){
        print(sneI)
        tryCatch({
            cSNe = heapList[[sneI]]
            res = sapply(cSNe$LCurves, FUN = modelObj$modelFit)
            if(bootstrap)
                res2 = lapply(cSNe$LCurves, FUN = modelObj$bootFit, nBoot = 100)
            if(!all(res))
                stop("max failed!")
            cSNe$plotFit(figdir)
            
        }, error = function(e){
            rmList <<- c(rmList, sneI)
        })
    }
    
    
    if(length(rmList) > 0)
        fittedList = heapList[-rmList]
    else
        fittedList = heapList
    
    return(list(fittedList = fittedList, modelObj = modelObj))
}


createHeapList = function(extCorr, folderpath, tablepath, figdir){
    filterN = 4
    sne_table = read_csv(tablepath)
    sneModel = sneLCModel$new()
    heapList = list()
    for(rowI in 1:dim(sne_table)[1]){
        cat(rowI,dim(sne_table)[1], "\r")
        flush.console()
        try({
            cSNe = sneReadIn(sne_table[rowI, ],
                             filterN = filterN,
                             folderpath = folderpath,
                             extCorr = extCorr)
            
            res = sapply(cSNe$LCurves, 
                         FUN = sneModel$maxFinderSpline)
            if(all(res)){
                cSNe$plotObs(figdir)
                heapList = c(heapList,cSNe)
            }
        })
    }
    return(heapList)
    
}

