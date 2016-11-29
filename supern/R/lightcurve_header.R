paramsClass <- R6Class(
    "paramsClass",
    portable = FALSE,
    public = list(
        phase_min = NA,
        phase_max = NA,
        tmax = NA,
        mmax = NA,
        dM15 = NA,
        betas = NA,
        betasSigma = NA,
        mcurve = NA,
        
        ## Bootstrap
        betas_boot = NA, 
        betas_boot_mean = NA,
        betas_boot_sd = NA,
        
        tmax_boot = NA, 
        mmax_boot = NA,
        mmax_boot_mean = NA, 
        mmax_boot_sd = NA,
        mcurve_boot = NULL, 
        mcurve_bootSD = NULL
    ))



lightcurve <- R6Class(
    "lightcurve",
    portable = FALSE,
    public = list(
        redshift = NA,
        bandName = "",
        mjd = numeric(0),
        phase = numeric(0),
        mag = numeric(0),
        sigma = numeric(0),
        code = character(0),
        params = NULL,
        
        initialize = function(mjd_, mag_, sigma_, 
                              code_, redshift_){
            mjd <<- mjd_
            mag <<- mag_
            sigma <<- sigma_
            code <<- code_
            redshift <<- redshift_
            params <<- paramsClass$new()
        },
        
        computePhase = function(){
            phase <<- (mjd - params$tmax) / (1 + redshift)
        },
        
        addObsPoints = function(pCol){
            points(mjd, mag, pch = 20, col = pCol)
            arrows(mjd, mag + sigma,
                   mjd, mag - sigma, 
                   length = 0.02, angle = 90,
                   col = pCol, code = 3)
            if (!is.na(params$tmax) && 
                !is.na(params$mmax)){
                points(params$tmax, params$mmax, cex = 2,
                       pch = 4, col = pCol, lwd = 2)
            }
        },
        
        addPhaseFit = function(){
            computePhase()
            addMaxLine()
            addPhasePoints()
            addFitLine()
            if(!is.null(params$mcurve_boot))
                addBootConf()
        },
        
        
        getObsMatInPhase = function(){
            computePhase()
            resMat = cbind(phase, mag, sigma)
            return(resMat)
        }
        
    ), # public 
    
    private = list(
        addPhasePoints = function(){
            points(phase, mag, pch = 20)
            arrows(phase, mag + sigma,
                   phase, mag - sigma, 
                   length = 0.02, angle = 90,
                   code = 3)
        },
        
        addFitLine = function(){
            bSeq = seq(0, 1, length.out = 100)
            pSeq = bSeq * (params$phase_max - params$phase_min) + 
                params$phase_min
            magSeq = params$mcurve(pSeq)
            lines(pSeq, magSeq)
        },
        
        addBootConf = function(){
            bSeq = seq(0, 1, length.out = 100)
            pSeq = bSeq * (params$phase_max - params$phase_min) + 
                params$phase_min
            magSeq = params$mcurve_boot(pSeq)
            sdSeq = params$mcurve_bootSD(pSeq)
            lines(pSeq, magSeq, lty = 2)
            lines(pSeq, magSeq + 2*sdSeq, lty = 3)
            lines(pSeq, magSeq - 2*sdSeq, lty = 3)
        },
        
        
        addMaxLine = function(){
            abline(h = params$mmax, lty =2)
            abline(v = 0, lty = 2)
        }
    ))


