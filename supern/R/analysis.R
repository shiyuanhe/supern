SNE_matchColor = function(bandN = NULL){
    bandNList = c("B", "V", "R", "I")
    cols = c("red","blue","green","black","purple")
    pos = match(bandN, bandNList)
    pos[is.na(pos)] = 5
    return(cols[pos])
}


SNE_matchShape = function(bandN = NULL){
    bandNList = c("B", "V", "R", "I")
    shapes = 1:5
    pos = match(bandN, bandNList)
    pos[is.na(pos)] = 5
    return(shapes[pos])
}

SNE_BVRIatBmax = function(LCurves,redshift){
    tBmax = LCurves[[1]]$params$tmax
    resMean = numeric(0)
    resSD = numeric(0)
    for(j in 1:length(LCurves)){
        t_max = LCurves[[j]]$params$tmax
        phase = (tBmax - t_max)/(1 + redshift)
        mCurve = LCurves[[j]]$params$mcurve
        sdCurve = LCurves[[j]]$params$mcurve_bootSD
        resMean = c(resMean, mCurve(phase))
        resSD = c(resSD,  sdCurve(phase))
    }
    res = list(resMean = resMean, resSD = resSD)
    return(res)
}


SNE_computeDM15 = function(LCurves){
    resV = numeric(0)
    resSD = numeric(0)
    for(j in 1:length(LCurves)){
        t_max = LCurves[[j]]$params$tmax
        #phase = (15 - t_max)/(1 + redshift)
        phase = 15
        mCurve = LCurves[[j]]$params$mcurve
        sdCurve = LCurves[[j]]$params$mcurve_bootSD
        dM15 = mCurve(phase) - mCurve(0)
        dM15SD = sdCurve(phase)#^2 + sdCurve(0)^2
        resV = c(resV, dM15)
        resSD = c(resSD, dM15SD)
    }
    res = list(resV = resV, resSD  = resSD)
    return(res)
}

SNE_getNewBeta = function(beta1, beta2, beta1sd, beta2sd, bandN){
    beta12 = c(beta1, beta2, beta1sd, beta2sd)
    if(bandN == "B"){
        res = SNE_getNewBeta_CORE(beta12, 
                                  manifoldMappingList[[1]]$x,
                                  manifoldMappingList[[1]]$y)
    }else{
        res = SNE_getNewBeta_CORE(beta12, 
                                  manifoldMappingList[[2]]$x,
                                  manifoldMappingList[[2]]$y)
    }
    return(res)
}

SNE_getNewBeta_CORE = function(beta12, beta1Seq, beta2Seq){
    dist = (beta12[1] - beta1Seq)^2 + (beta12[2] - beta2Seq)^2
    dist = sqrt(dist)
    kmin = which.min(dist)
    newbeta1 = kmin/length(beta1Seq)
    newbeta2 = min(dist) * sign(beta12[2] - beta2Seq[kmin])
    
    k = which.min(dist)
    if(k == length(dist)) k = k - 1
    deriv = (beta2Seq[k+1] - beta2Seq[k]) / (beta1Seq[k+1] - beta1Seq[k])
    theta = atan(deriv)
    sigma1 = sqrt(cos(theta)^2 * beta12[3]^2 + sin(theta)^2* beta12[4]^2)
    sigma1 = sigma1 / 14
    sigma2 = sqrt(cos(theta)^2 * beta12[4]^2 + sin(theta)^2* beta12[3]^2)
    newbeta = c(newbeta1, sigma1, newbeta2,  sigma2)
    names(newbeta) = NULL
    return(newbeta)
}




# pairMatch = function(allScores, cutOff){
#     dist_mat = as.matrix(dist(allScores))
#     np = dim(dist_mat)[1]
#     diag(dist_mat) = 1e5
#     sne_unused = rep(TRUE, np)
#     k = 0
#     pair = rep(-1, np)
#     pair_dist =rep(1e5, np)
#     while(sum(sne_unused)>1){
#         m = min(dist_mat)
#         if(m > cutOff) break()
#         p = which.min(dist_mat) - 1
#         i = p %% np + 1
#         j = p %/% np + 1
#         dist_mat[i,] = 1e10
#         dist_mat[,j] = 1e10
#         dist_mat[j,] = 1e10
#         dist_mat[,i] = 1e10
#         sne_unused[i] = FALSE
#         sne_unused[j] = FALSE
#         k = k + 1
#         pair[i] = k
#         pair[j] = k
#         pair_dist[i] = m
#         pair_dist[j] = m
#     }
#     res = cbind(pair, pair_dist)
#     return(res)
# }
# 
# 
