combineAsVector = function(colnames){
    X = allInfo[, colnames]
    X = as.matrix(X); X = as.vector(X)
    return(X)
}


getXLim = function(X, Xnew = NULL){
    xmin = min(X, na.rm = TRUE)
    xmax = max(X, na.rm = TRUE)
    if(!is.null(Xnew)){
        xnew_min = min(Xnew, na.rm = TRUE)
        xnew_max = max(Xnew, na.rm = TRUE)
        xmin = min(xmin, xnew_min)
        xmax = max(xmax, xnew_max)
    }
    xlim = c(xmin, xmax)
    return(xlim)
}

plotXYSe = function(X, Y, xlab, ylab,
                    xyCol,
                    Xse = NULL, Yse = NULL, 
                    newX  = NULL, newY = NULL, 
                    newXse = NULL, newYse  = NULL, 
                    newCols = NULL){
    xlim = getXLim(X, newX)
    ylim = getXLim(Y, newY)
    plot(0, 0, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, type = "n")
    if(!is.null(Xse)){
        arrows(X - Xse, Y, X + Xse, Y, length = 0.02, 
               angle = 90, code = 3, col = "grey")
    }
    if(!is.null(Yse)){
        arrows(X, Y - Yse, X, Y + Yse, length = 0.02, 
               angle = 90, code = 3, col = "grey")
    }
    if(!is.null(newCols)){
        newCols = SNE_matchColor(newCols)
    }else if(!is.null(newX)){
        newCols = rep("black", length(newX))
    }
    
    if(!is.null(newXse)){
        arrows(newX - newXse, newY, 
               newX + newXse, newY, 
               length = 0.02, angle = 90, 
               code = 3, col = newCols)
    }
    if(!is.null(newYse)){
        arrows(newX, newY - newYse, 
               newX, newY + newYse, 
               length = 0.02, angle = 90, 
               code = 3, col = newCols)
    }
    points(X, Y, col = "purple", pch = SNE_matchShape(xyCol))
    if(!is.null(newX)){
        points(newX, newY, pch = 19, 
               col = newCols, cex = 2)
    }
}


SNE_plotSpcClass = function(plotdata, xlab, ylab, newdata = NULL){
    if(!is.null(newdata)){
        plotdata = rbind(newdata, plotdata)
        pColor = scale_colour_manual(values = 
                                         c("black","red","blue", "green", "purple"))
        pShape = scale_shape_manual(values = c(19, 1, 2, 0, 6))
    }else{
        pColor = scale_colour_manual(values = 
                                         c("red","blue", "green", "purple"))
        pShape = scale_shape_manual(values = c(1, 2, 0, 6))
    }
    p = ggplot(plotdata, aes(s1, s2,shape = type, col = type)) +
        geom_point(size = 3)
    p = p + pColor + pShape
    p = p + geom_errorbar(aes(ymin = s2 - s2SD, ymax = s2 + s2SD))
    p = p + geom_errorbarh(aes(xmin = s1 - s1SD, xmax = s1 + s1SD))
    p = p + theme_bw() + xlab(xlab) + ylab(ylab)
    p
}





