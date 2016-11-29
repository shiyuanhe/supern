sneLCModel <- R6Class(
    "sneLCModel",
    portable = FALSE,
    public = list(
        phase_min = NA,
        phase_max = NA,
        mean_fun = NULL,  ## mean functions
        tfunctions = NULL,  ## template functions
        tfunctions2D = NULL,  ## template functions, second derivatives
        SIGMA = NULL,
        SIGMA_FULL = NULL,
        LCurve = NULL,
        params = NULL,
        constr_A = NULL,
        constr_b0 = NULL,
        
        mjdCut = NULL,
        magCut = NULL,
        sigmaCut = NULL,
        phase = NULL,
        mu0Point = NULL,
        PhisPoint = NULL,
        Gamma_Inv = NULL,
        
        
        template_fun = function(phase4basis){
            tmp = lapply(tfunctions, 
                         function(tfun) tfun(phase4basis))
            res = do.call(cbind, tmp)
            return(res)
        },
        
        template2D_fun = function(phase4basis){
            tmp = lapply(tfunctions2D, 
                         function(tfun) tfun(phase4basis))
            res = do.call(cbind, tmp)
            return(res)
        },
        
        modelFit = function(LCurve_){
            setLightcurve(LCurve_)
            succ = lcBasicCheck(2,1)
            succ2 = maxFinderInternal()
            succ = succ2 || succ
            computePhase(phase_min, phase_max)
            computeModel()
            lowRankLSFit()
            LRLSFIT_constraint()
            LRLSFIT_closeFit()
            # maxAdjustAfterFit()
            # computePhase(phase_min, phase_max)
            # computeModel()
            # lowRankLSFit()
            return(succ)
        },
        
        # bootstrap after model fit
        bootFit = function(LCurve_, nBoot){
            bootInitParams(LCurve_, nBoot)
            pp_grid = seq(phase_min, phase_max, length.out = 300)
            allFits = matrix(0, 300, nBoot)

            ylim = c(max(LCurve_$mag),  min(LCurve_$mag))
            plot(LCurve_$mjd, LCurve_$mag, ylim = ylim)

            for(b in 1:nBoot){
                newLC = bootNewLC(LCurve_)
                newLC$addObsPoints(pCol = "grey")
                modelFit(newLC)
                bootUpdateParams(LCurve_, newLC, b)
                allFits[, b] = newLC$params$mcurve(pp_grid)

            }
            yy_grid = apply(allFits, 1, mean)
            LCurve_$params$mcurve_boot = approxfun(pp_grid, yy_grid)
            yy_grid = apply(allFits, 1, sd)
            LCurve_$params$mcurve_bootSD = approxfun(pp_grid, yy_grid)
            bootClose(LCurve_)
        },

                
        maxFinder = function(LCurve_){
            setLightcurve(LCurve_)
            succ = maxFinderInternal()
            return(succ)
        },
        
        maxFinderSpline = function(LCurve_){
            setLightcurve(LCurve_)
            mjd = LCurve_$mjd
            mag = LCurve_$mag
            t0 = mjd[which.min(mag)]
            tSeq = seq(t0 - 5, t0 + 5, by = 0.1)
            ySeq = sapply(tSeq, splineLocalFit, 
                          mjd = mjd, mag = mag,
                          sigma = LCurve_$sigma)
            bestPos = which.min(ySeq)
            params$mmax <<- ySeq[bestPos]
            params$tmax <<- tSeq[bestPos]
            succ = maxFinderCheck()
            return(succ)
        },
        
        splineLocalFit = function(t0, mjd, mag, sigma){
            Xtmp = data.frame(mjd = mjd - t0, 
                              mjd2 = (mjd - t0)^2,
                              mag = mag)
            w = exp(-abs(Xtmp$mjd)/2.5) / sigma^2
            localfit = lm(mag ~ mjd + mjd2, Xtmp, weights = w)
            Ypred = coef(localfit)[1]
            return(Ypred)
        },
        
        
        trainingMean = function(heapList, basisObj,
                                lambda, cv = TRUE){
            obsMatStack = heapObsM_asMat(heapList, basisObj, 
                                         weighting = TRUE)
            magStack = obsMatStack[ ,"mag"]
            basisStack = heapBasisM_asMat_weighted(heapList, basisObj)
            residuals = trainingMeanInternal(magStack, basisStack, 
                                             basisObj$OmegaCC, lambda, cv)
            calMeanFun(basisObj)
            constrVector(basisObj)
            return(residuals)
        },
        
        
        ## functional reduced rank model
        trainingFPCA = function(heapList, basisObj,
                                lambda1, lambda3, tK){
            magList = heapObsMag_asList(heapList, basisObj, TRUE, mean_fun)
            basisList = heapBasisM_asList_weighted(heapList, basisObj)
            fpcaInternal(magList, basisList, 
                         lambda1, lambda3,
                         basisObj$OmegaCC, tK)
            calTempFun(basisObj, tK)
            constrMatrix(basisObj)
        },
        
        constrMatrix = function(basisObj){
            knots_phase = basisObj$bknots * (phase_max - phase_min) + phase_min
            knots_sel1 = subset(knots_phase, knots_phase < 0 & knots_phase > phase_min)
            knots_sel2 = subset(knots_phase, knots_phase > 35 & knots_phase < phase_max)
            constr_A <<- rbind(-basisObj$bD1EvalAt(knots_sel1) %*% svdU,
                               basisObj$bD2EvalAt(knots_sel1) %*% svdU,
                               basisObj$bD1EvalAt(knots_sel2) %*% svdU)
            constr_A <<- t(constr_A)
        },
        
        constrVector = function(basisObj){
            knots_phase = basisObj$bknots * (phase_max - phase_min) + phase_min
            knots_sel1 = subset(knots_phase, knots_phase < 0 & knots_phase > phase_min)
            knots_sel2 = subset(knots_phase, knots_phase > 35 & knots_phase < phase_max)
            constr_b0 <<- c(basisObj$bD1EvalAt(knots_sel1) %*% thetaHat,
                            -basisObj$bD2EvalAt(knots_sel1) %*% thetaHat,
                            -basisObj$bD1EvalAt(knots_sel2) %*% thetaHat)
        },
        
        
        plotMeanC = function(heapList){
            tmp = heapObsM_asMat(heapList, basisObj, weighting = FALSE)
            plot(tmp[,1], tmp[,2], xlab = "Phase",
                 ylab = "Mag", pch = 20, col ="grey")
            tt = seq(phase_min, phase_max, length.out = 100)
            lines(tt, mean_fun(tt))
        },
        
        plotTempC = function(tJ){
            tmp = heapObsM_asMat(heapList, basisObj, 
                                 weighting = FALSE, mean_fun)
            plot(tmp[,1], tmp[,2], 
                 xlab = "Phase", ylab = "Mag", 
                 pch = 20, col = "grey")
            tt = seq(phase_min, phase_max, length.out = 100)
            yy = tfunctions[[tJ]](tt)
            lines(tt,  sqrt(SIGMA[tJ,tJ]) * 2 * yy)
        }
        
    ),## public
    
    
    private = list(
        thetaHat = NA,
        svdU = NA,
        ## functional reduced rank model
        fpcaInternal = function(magList, basisList,
                                lambda1, lambda3, Omega, 
                                nK){
            n = length(magList)
            fpca_res = nuclearFPCA(magList, basisList, Omega,
                                   lambda1, lambda3)
            svdS = svd(fpca_res$S)
            scores = sweep(svdS$v, MARGIN = 2, svdS$d, "*")
            SIGMA <<- var(scores)[1:nK, 1:nK]
            SIGMA_FULL <<- var(scores)
            svdU <<- svdS$u[,1:nK]
            
            ## Sign Adjustment, positive first component
            for(k in 1:nK){
                if(svdU[1,k]<0)
                    svdU[,k] <<- -svdU[,k]
            }
        },
        
        calTempFun = function(basisObj, nTmpBasis){
            phase_min <<- basisObj$phase_min
            phase_max <<- basisObj$phase_max
            xx = seq(0, 1, length.out = 200)
            tt = xx * (basisObj$phase_max - basisObj$phase_min) + 
                basisObj$phase_min
            tfunctions <<- list()
            tfunctions2D <<- list()
            for(t in 1:nTmpBasis){
                yy = basisObj$bD0EvalAt(tt) %*% svdU[,t]
                tfunctions[[t]] <<- approxfun(tt, yy)
                yy = basisObj$bD2EvalAt(tt) %*% svdU[,t]
                tfunctions2D[[t]] <<- approxfun(tt, yy)
            }
        },
        
        
        ### Training Mean Fuction 
        trainingMeanInternal = function(y, X, Omega, lambda, cv){
            numObs = length(y)
            thetaHat <<- calThetaHat(y, X, Omega, lambda)
            
            residuals = 0
            if(cv){
                for(i in 1:100){
                    bootSel = sample(1:numObs, numObs, replace = T)
                    thetaCV = calThetaHat(y[bootSel], X[bootSel,], 
                                          Omega, lambda)
                    oobSel = 1:numObs %in% bootSel
                    oobSel = which(!oobSel)
                    residuals = residuals + 
                        calRSSHat(y[oobSel], X[oobSel, ], thetaCV)
                }
            }
            return(residuals)
        },
        
        calThetaHat = function(y, X, Omega, lambda){
            tmp = t(X) %*% X + lambda * Omega
            theta = solve(tmp) %*% t(X) %*% y
            return(theta)
        },
        
        
        calRSSHat = function(y, X, thetaHat){
            resids = y - X %*% thetaHat
            rss = sum(resids^2)
            return(rss)
        },
        
        calMeanFun = function(basisObj){
            phase_min <<- basisObj$phase_min
            phase_max <<- basisObj$phase_max
            xx = seq(0,1,length.out = 100)
            tt = xx * (basisObj$phase_max - basisObj$phase_min) + 
                basisObj$phase_min
            yy = basisObj$bD0EvalAt(tt) %*% thetaHat
            mean_fun <<- approxfun(tt, yy)
        },
        
        
        ## peak quality check before fitting
        lcBasicCheck = function(cut1, cut2){
            peakRow = which.min(LCurve$mag)
            peakT = LCurve$mjd[peakRow]
            sel1 = (LCurve$mjd > peakT - 5)  & 
                (LCurve$mjd < peakT + 5)
            sel2 = (LCurve$mjd > peakT - 10) & 
                (LCurve$mjd < peakT)
            sel3 = (LCurve$mjd > peakT) & 
                (LCurve$mjd < peakT + 10)
            if (sum(sel1) < cut1 || 
                sum(sel2) < cut2 || 
                sum(sel3) < cut2){
                return(FALSE)
            }else{
                return(TRUE)
            }
        },
        
        maxAdjustAfterFit = function(){
            pSeq = seq(-5, 5, length.out = 1000)
            bSeq = (pSeq - phase_min) / (phase_max - phase_min)
            ySeq = params$mcurve(pSeq)
            k = which.min(ySeq)
            params$mmax <<- ySeq[k]
            params$tmax <<- params$tmax + 
                pSeq[k] * (1 + LCurve$redshift) 
        },
        
        
        setLightcurve = function(LCurve_){
            LCurve <<- LCurve_
            params <<- LCurve$params
            params$phase_min <<- phase_min
            params$phase_max <<- phase_max
        },
        
        computeModel = function(){
            phase4basis = phase
            mu0Point <<- mean_fun(phase4basis)
            PhisPoint <<- template_fun(phase4basis)
            Gamma = PhisPoint %*% SIGMA %*% t(PhisPoint)
            diag(Gamma) = diag(Gamma) + sigmaCut^2
            Gamma_Inv <<- solve(Gamma)
        },
        
        computeMax = function(){
            ytilde = magCut - mu0Point
            ## mmax_hat is (1^T*GammaInv*ytilde) / (1^T*GammaInv*1)
            miHat = sum(Gamma_Inv %*% ytilde) / sum(Gamma_Inv)
            ytilde = ytilde - miHat
            nCut = length(magCut)
            aveResid = 1/nCut * 
                as.numeric(t(ytilde) %*% Gamma_Inv %*% ytilde)
            return(list(miHat = miHat, 
                        aveResid = aveResid))
        },
        
        computePhase = function(cut_min, cut_max){
            phase <<- (LCurve$mjd - params$tmax) / 
                (1 + LCurve$redshift)
            sel <- (phase > cut_min) & (phase < cut_max)
            phase <<- subset(phase, sel)
            mjdCut <<- subset(LCurve$mjd, sel)
            magCut <<- subset(LCurve$mag, sel)
            sigmaCut <<- subset(LCurve$sigma, sel)
        },
        
        templateSeq = function(){
            #qSeq = seq(0, 1, length.out = 1000)
            pSeq = seq(phase_min, phase_max, length.out = 1000)
            muSeq = mean_fun(pSeq)
            PhisSeq = template_fun(pSeq)
            return(cbind(pSeq, muSeq, PhisSeq))
        },
        
        lowRankLSFit = function(){
            ## compute betas_hat
            yDM <- magCut - mu0Point - params$mmax
            XC = PhisPoint
            betas <- SIGMA %*% t(XC)%*% Gamma_Inv %*% yDM
            params$betas <<- betas
            ## compute betas_var
            W2_Inv <- diag(1/sigmaCut^2)
            var_betas_Inv <- solve(SIGMA) +
                t(PhisPoint) %*% W2_Inv %*% PhisPoint
            params$betasSigma <<- solve(var_betas_Inv)
        },
        
        LRLSFIT_constraint = function(){
            Dmat = solve(params$betasSigma)
            dvec = Dmat %*% params$betas
            solve_res = solve.QP(Dmat, dvec, 
                                 constr_A, constr_b0)
            params$betas <<- solve_res$solution
        }, 
        
        LRLSFIT_closeFit = function(){
            tpSeq <- templateSeq()
            ySeq <- params$mmax + 
                tpSeq[,-1] %*% c(1, params$betas)
            params$mcurve <<- approxfun(tpSeq[,1], ySeq)
            params$dM15 <<- params$mcurve(15)- params$mcurve(0)
        },
        
        getTSeq = function(){
            peakGuess = which.min(LCurve$mag)
            ti0 = LCurve$mjd[peakGuess]
            tseq = seq(from = ti0 - 3, to = ti0 + 3, by = 0.1)
            return(tseq)
        },
        
        maxFinderInternal = function(){
            tseq = getTSeq()
            mi =  rep(0, length(tseq)) # peak magnitude
            tmax_gof = rep(1e10, length(tseq)) # goodness of fit
            for (j in 1:length(tseq)){
                try({
                    params$tmax <<- tseq[j]
                    computePhase(-10, 20)
                    if (length(mjdCut)<5)
                        next()
                    computeModel()
                    res = computeMax()
                    tmax_gof[j] = res$aveResid
                    mi[j] = res$miHat
                })
            }
            bestPos = which.min(tmax_gof)
            if (tmax_gof[bestPos]>1000) 
                return(FALSE)
            params$mmax <<- mi[bestPos]
            params$tmax <<- tseq[bestPos]
            succ = maxFinderCheck()
            return(succ)
        },
        
        maxFinderCheck = function(){
            computePhase(-10, 20)
            cnt1 = sum(self$phase >= (-5) & self$phase <= 0)
            cnt2 = sum(self$phase >= 0  & self$phase <= 5)
            cnt3 = sum(self$phase >= (-2.5)  & self$phase <= 2.5)
            if (!(cnt1 > 0 && cnt2 > 0 && cnt3 > 0))
                return(FALSE)
            return(TRUE)
        },
        
        ## Bootstrap Internal
        
        bootInitParams = function(LC, nBoot){
            nPCs = length(LC$params$betas)
            LC$params$betas_boot = matrix(0, nPCs, nBoot)
            LC$params$tmax_boot = rep(0, nBoot)
            LC$params$mmax_boot = rep(0, nBoot)
        },
        
        bootNewLC = function(LC){
            mjd = LC$mjd
            sigma = LC$sigma
            code = LC$code
            phase = (LC$mjd - LC$params$tmax) / (LC$redshift + 1)
            
            sigma1 = sigma
            sel = sigma1 > 0.1
            sigma1[sel] = 0.1

            nPoints = length(mjd)
            mag_new = LC$params$mcurve(phase) + 
                rnorm(nPoints, 0, sigma1)
            sel = !is.na(mag_new)
            newLC = lightcurve$new(mjd[sel], mag_new[sel], 
                                   sigma[sel],
                                   code[sel], LC$redshift)
            newLC$bandName = LC$bandName
            return(newLC)
        },
        
        
        bootUpdateParams = function(LC, newLC, b){
            LC$params$betas_boot[,b] = newLC$params$betas
            LC$params$tmax_boot[b] = newLC$params$tmax
            LC$params$mmax_boot[b] = newLC$params$mmax
        },
        
        bootClose = function(LC){
            LC$params$mmax_boot_mean = mean(LC$params$mmax_boot)
            LC$params$mmax_boot_sd = sd(LC$params$mmax_boot)
            LC$params$betas_boot_mean = apply(LC$params$betas_boot, 1, mean)
            LC$params$betas_boot_sd = apply(LC$params$betas_boot, 1, sd)
        }
        
        
        
    )
)## sneLCModel

