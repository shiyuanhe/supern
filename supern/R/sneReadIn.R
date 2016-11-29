sneReadIn <- function(tableRow, filterN,
                      folderpath = "./data/data_raw_convert/", 
                      extCorr = FALSE){
    suffix = "_main.csv"
    if(filterN == 2)
        filterNames = c("B","V")
    if(filterN == 4)
        filterNames = c("B","V","R","I","r","i","r'","i'")
    sneName = tableRow$SN
    sneSurvey = tableRow$Survey
    redshift = tableRow$Zcmb
    filepath = paste0(folderpath,
                      sneName,"_",sneSurvey,
                      suffix)
    if(!file.exists(filepath))
        stop(paste("Files Not Found:", filepath))
    cSNe = supernova$new(sneName, sneSurvey, redshift)
    cSNe$readFile(filepath, filterNames)
    if(extCorr){
        #evs = getExtVector(tableRow, cSNe$LCBandNames)
        bn = cSNe$LCBandNames
        for(i in 1:length(bn)){
            if(substr(bn[i],1,1) == "r")
                bn[i] = "R"
            if(substr(bn[i],1,1) == "i")
                bn[i] = "I"
        }
        evs = getExtVector(tableRow, bn)
        cSNe$corrExt(evs)
    }
    return(cSNe)
}



getExtVector <- function(tableRow, filterNames){
    resVec = numeric(0)
    for(i in 1:length(filterNames)){
        tmp = switch (
            substr(filterNames[i],1,1),
            B = tableRow$AB,
            V = tableRow$AV,
            R = tableRow$AR,
            I = tableRow$AI,
            r = tableRow$sdssr,
            i = tableRow$sdssi
        )
        resVec[i] = tmp
    }
    return(resVec)
}
