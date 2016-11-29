
supernova <- R6Class(
    "supernova",
    portable = FALSE,
    public = list(
        sneName = "",
        obsProj = "",
        redshift = numeric(0),
        extinction = numeric(0), # milky way extinction value
        numBands = 0,
        LCBandNames = character(0),
        LCurves = list(),
        
        initialize = function(sneName_, obsProj_, redshift_){
            sneName <<- sneName_
            obsProj <<- obsProj_
            redshift <<- redshift_
        },
        
        readFileInternal = function(filePath, selFilters){
            sneObs <- read.csv(filePath, header=TRUE, 
                               stringsAsFactors = F)
            sel = sneObs$Mag < 30
            sneObs = subset(sneObs, sel)
            for(j in 1:length(selFilters)){
                sel <- (sneObs$Filter == selFilters[j])
                if(sum(sel) < 2) next()
                addBand(sneObs$MJD[sel], sneObs$Mag[sel],
                        sneObs$Sigma[sel], sneObs$Code[sel],
                        selFilters[j])
            }
        },
        
        extractField = function(sneObs, fieldName){
            selLine = grepl(fieldName, sneObs, ignore.case = TRUE)
            info = sneObs[selLine]
            info = gsub(fieldName, "", info, ignore.case = T)
            info = gsub("\\s","", info, perl = TRUE)
            return(info)
        },
        
        extractData = function(sneObs){
            selLine = grepl("#data", sneObs, ignore.case = TRUE)
            lineNum = which(selLine)
            lineB = lineNum + 1
            lineE = length(sneObs)
            sneData = sneObs[lineB:lineE]
            sneData = paste(sneData,  collapse = "\r")
            sneData = read.table(header = FALSE, text = sneData, 
                               stringsAsFactors = FALSE)
            colnames(sneData) = c("Filter", "MJD", "Mag", "Sigma")
            return(sneData)
        },
        
        
        parameterTable = function(){
            res1 = tablePartI()
            res2 = tablePartII()
            res3 = tablePartIII()
            res = rbind(res1, res2, res3)
            colnames(res) = cSNe$LCBandNames
            return(res)
        },
        
        readFile = function(filePath){
            #filePath = "~/Downloads/SN2008fp_CSP_main.csv"
            con = file(filePath)
            sneObs = readLines(con)
            close(con)
            sneName <<- extractField(sneObs, "#name")
            obsProj <<- extractField(sneObs, "#survey")
            redshift <<- extractField(sneObs, "#redshift")
            redshift <<- as.numeric(redshift)

            sneData = extractData(sneObs)
            sneData$Filter = toupper(sneData$Filter)
            sneData$Filter = substr(sneData$Filter,1,1)
            for(ff in c("B", "V", "R", "I")){
                sel <- (sneData$Filter == ff)
                if(sum(sel) < 2) next()
                addBand(sneData$MJD[sel], sneData$Mag[sel],
                        sneData$Sigma[sel], "CodeNA",
                        ff)
            }
        },
        
        plotObs = function(plotFolder){
            setupPlot(plotFolder,"Obs")
            initObsPlot()
            cols = c("red","blue","green","black","purple")
            for(i in 1:numBands){
                LCurves[[i]]$addObsPoints(cols[i])
            }
            closePlot(plotFolder)
        },
        
        plotFit = function(plotFolder, plotID = NULL){
            if(is.null(plotID)){
                plotID = 1:numBands
            }
            for(i in plotID){
                setupPlot(plotFolder, LCBandNames[i])
                initPhasePlot(LCBandNames[i])
                LCurves[[i]]$addPhaseFit()
                closePlot(plotFolder)
            }
        },
        
        corrExt = function(evs){
            extinction <<- evs
            res = mapply(function(onelc, evalue){
                            onelc$mag = onelc$mag - evalue
                            }, LCurves, evs)
        }
        
    ), ##public
    private = list(
        yObsRange = function(){
            ymin = 1e100
            ymax = -1e100
            for(i in 1:(self$numBands)){
                ymin  = min(ymin, LCurves[[i]]$mag)
                ymax  = max(ymax, LCurves[[i]]$mag)
            }
            ylimObs = c(ymax, ymin)
            return(ylimObs)
        },
        
        xObsRange = function(){
            xmin = 1e100
            xmax = -1e100
            for(i in 1:(self$numBands)){
                xmin  = min(xmin, LCurves[[i]]$mjd)
                xmax  = max(xmax, LCurves[[i]]$mjd)
            }
            xlimObs = c(xmin, xmax)
            return(xlimObs)
        },
        
        xPhaRange = function(){
            xmin = LCurves[[1]]$params$phase_min
            xmax = LCurves[[1]]$params$phase_max
            xlimObs = c(xmin, xmax)
            return(xlimObs)
        },
        
        initObsPlot = function(){
            ylimObs = yObsRange()
            xlimObs = xObsRange()
            plot(0,0, xlab="MJD",
                 ylab = "Magnitude", type="n",
                 xlim = xlimObs, ylim = ylimObs,
                 main = paste0(sneName,"-",obsProj))
        },
        
        initPhasePlot = function(filter){
            ylimObs = yObsRange()
            xlimPhase = xPhaRange()
            plot(0, 0, xlab="Phase",
                 ylab = "Magnitude", type="n",
                 xlim = xlimPhase, ylim = ylimObs,
                 main = paste0(sneName,"-",obsProj,
                               "-", filter))
        },
        
        
        setupPlot = function(plotFolder, nameSufix){
            if (nchar(plotFolder)>0){
                outfile = paste0(plotFolder,sneName,"-",
                                 obsProj ,"-", 
                                 nameSufix, ".png")
                png(outfile, height=800,width = 900)
            }
        },
        
        closePlot = function(plotFolder){
            if (nchar(plotFolder)>0)
                dev.off()
        },
        
        # New band with light curve
        addBand = function(mjd, mag, sigma, 
                           code, bandName){
            LCBandNames <<- c(LCBandNames, bandName)
            numBands <<- numBands + 1
            newLC = lightcurve$new(mjd, mag, sigma,
                                   code,redshift)
            newLC$bandName = bandName
            LCurves[[numBands]] <<- newLC
        },
        
        tablePartI = function(){
            res = sapply(LCurves, 
                         FUN = function(LC){
                             c(LC$params$tmax,
                               LC$params$mmax, 
                               LC$params$mmax_boot_sd)})
            rownames(res) = c("MJD_max","Mag_max", 
                              "Mag_max_SD")
            return(res)
        },
        
        tablePartII = function(){
            tmp = SNE_BVRIatBmax(LCurves, redshift)
            res = rbind(tmp$resMean, tmp$resSD)
            tmp = SNE_computeDM15(LCurves)
            res = rbind(res, tmp$resV, tmp$resSD)
            rownames(res) = c("Mag_at_Bmax", "mag_at_Bmax_SD",
                              "DeltaM15", "DeltaM15_SD")
            return(res)
        }, 
        
        tablePartIII = function(){
            res = sapply(LCurves, 
                         FUN = function(LC){
                             c(LC$params$betas,  ##betas_boot_mean or betas
                               LC$params$betas_boot_sd) })
            res = matrix(as.vector(res), nrow = 14)
            rownames(res) = c(paste0("beta", 1:7),
                              paste0("beta", 1:7,"_SD"))
            seqR1 = seq(1,14, by = 2)
            seqR2 = seq(2,14, by = 2)
            tmp1 = as.matrix(res[1:7,])
            tmp2 = as.matrix(res[8:14,])
            res[seqR1,] = tmp1
            rownames(res)[seqR1] = rownames(tmp1)
            res[seqR2,] = tmp2
            rownames(res)[seqR2] = rownames(tmp2)
            return(res)
        }
        
    )#private
)#supernova class

