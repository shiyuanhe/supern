#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(supern)
library(ggplot2)

NLcurves = 4

# Define UI for application that draws a histogram
ui <- shinyUI(bootstrapPage(
    
    # Application title
    titlePanel("Type Ia Supernova Light Curve Fitting"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            h3("Upload File"),
            fileInput('file1', 'Choose file to upload',
                      accept = c(
                          'text/tab-separated-values',
                          'text/plain',
                          '.txt',
                          '.dat'
                      )),
            h3("File Format"),
            HTML("
<ol><li>Provide SNe name, survey, redshift in file head.</li>
<li>The light curve data begins with '#Data'.</li>
<li>The light curve data is space (or tab) seperated with columns: filter, MJD, Mag, Mag_err.</li>
<li>The light curves should have standard filters: B, V, R, I.</li>
<li> Download an example file <a href='https://raw.githubusercontent.com/shiyuanhe/supern/master/data/SN2008fp_CSP_main.txt'>SN2008fp_CSP_main.txt</a>.</li>
</ol>
<br/><br/>
<b>Example:</b><br/>
#Name sn2008fp<br/>
#Redshift 0.02216<br/>
#Survey CSP<br/>#Data<br/>
V 54724.323 14.249 0.006<br/>
 V 54725.320 14.145 0.006<br/>
...<br/><br/>
                 ")
        ),
        mainPanel(
            h3("Observation"),
            plotOutput('obsPlot'),
            h3("Light curve fit"),
            uiOutput("fitPlots"),
            h3("Parameters"),
            p(paste0(
                "This table contains the information of each band:",
                " (1) Peak MJD, (2) Peak Magnitude",
                "(3) Magnitude at B Max (4) Delta M_15 and (5) 7 model scores."
            )),
            tableOutput("parameter"),
            h3("First Two Scores"),
            plotOutput("beta12Plots"),
            plotOutput("beta12Legend", height = "100%"),
            br(),
            p("Nonlinear scores via projecting to beta_1 and beta_2 curves."),
            tableOutput("nonlinearBeta12"),
            p("The probability that this is a Type Ia Supernova."),
            tableOutput("sniaProb"),
            h3("Intrinsic Color"),
            plotOutput("colorPlots"),
            p("Estimation of the intrinsic color:"),
            tableOutput("colorTable"),
            h3("Spectral Class Classification"),
            p("The sepctral classes are from the work of
              Benetti et al. (2005), Branch et al. (2009) 
              and Wang et al. (2009). The new observation 
              is plotted as black filled circle."),
            plotOutput("Spc1"),
            plotOutput("Spc2"),
            plotOutput("Spc3")
        )
        
    )
    
))


# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {
    
    output$obsPlot <- renderPlot({
        inFile <- input$file1
        if (is.null(inFile)) return(NULL)
        
        sneName = "x"
        sneSurvey = "y"
        redshift = 0.001
        
        cSNe <<- supernova$new(sneName, sneSurvey, redshift)
        cSNe$readFile(inFile$datapath)
        res = sapply(cSNe$LCurves, FUN = modelObj$modelFit)
        res2 = lapply(cSNe$LCurves, FUN = modelObj$bootFit, nBoot = 10)
        
        parameterTable <<- cSNe$parameterTable()
        
        nonlinearParam = numeric(0)
        for(ff in cSNe$LCBandNames){
            tmp = SNE_getNewBeta(parameterTable["beta1",ff],
                                 parameterTable["beta2",ff], 
                                 parameterTable["beta1_SD",ff],
                                 parameterTable["beta2_SD",ff],
                                 ff)
            nonlinearParam = cbind(nonlinearParam, tmp)
        }
        colnames(nonlinearParam) = cSNe$LCBandNames
        rownames(nonlinearParam) = c("nonlinear_beta1", "nonlinear_beta1_SD",
                                     "nonlinear_beta2", "nonlinear_beta2_SD")
        nonlinearParam <<- nonlinearParam
        
        NLcurves <<- length(cSNe$LCurves)
        cSNe$plotObs("")
    })
    
    
    output$parameter <- renderTable({
        inFile <- input$file1
        if(is.null(inFile)) return(NULL)
        
        return(parameterTable)
    }, rownames = TRUE, digits = 3)
    
    output$fitPlots <- renderUI({
        inFile <- input$file1
        if (is.null(inFile)) return(NULL)
        
        plot_output_list <- lapply(1:4, function(i) {
            plotname <- paste("plot", i, sep="")
            plotOutput(plotname)
        })
        
        do.call(tagList, plot_output_list)
    })
    
    output$beta12Plots <- renderPlot({
        inFile <- input$file1

        
        bandN = c("B", "V", "R", "I")
        nSample = dim(allInfo)[1]
        eCols = as.vector(sapply(bandN, function(tmp) rep(tmp, nSample)))
        
        X = combineAsVector(paste0(bandN, "_scores1"))
        Y = combineAsVector(paste0(bandN, "_scores2"))
        Xse = combineAsVector(paste0(bandN, "_scoresSD1"))
        Yse = combineAsVector(paste0(bandN, "_scoresSD2"))
        
        if (is.null(inFile)){
            plotXYSe(X, Y, xlab = "beta1", ylab = "beta2",
                     xyCol = eCols,
                     Xse = Xse, Yse = Yse)
            return(NULL)         
        }
        
        newX = parameterTable["beta1",]
        newY = parameterTable["beta2",]
        newXse = parameterTable["beta1_SD",]
        newYse = parameterTable["beta2_SD",]
        
        xlab = expression(beta^{(1)})
        ylab = expression(beta^{(2)})
        
        plotXYSe(X, Y, xlab = xlab, ylab = ylab,
                 xyCol = eCols, Xse = Xse, Yse = Yse,
                 newX = newX, newY = newY, 
                 newXse = newXse, newYse = newYse,
                 newCols = cSNe$LCBandNames)
    })
    
    
    output$beta12Legend <- renderPlot({
        par(mar = rep(0,4))
        plot(0,0,type = "n", xaxt = "n", yaxt = "n", bty = "n")
        legend(-0.9,0.6, legend = c("B", "V", "R", "I"), 
               pch = 1:4, col = "purple", horiz = TRUE)
        legend(0,0.6, legend = paste0(c("B", "V", "R", "I"), "-new"),
               pch = 20, col = c("red","blue","green", "black"),horiz = TRUE)
    }, height = 80)
    output$nonlinearBeta12 <- renderTable({
        inFile <- input$file1
        if (is.null(inFile)) return(NULL)
        
        nonlinearParam
    }, rownames = TRUE, digits = 3)
    
    output$sniaProb <- renderTable({
        
        inFile <- input$file1
        if(is.null(inFile)) return(NULL)
        
        newX = parameterTable["beta1",]
        newY = parameterTable["beta2",]
        newBeta12 = c(newX, newY)
        Z = newBeta12 - beta12_mu
        Zscalar = t(Z) %*% solve(beta12_cov) %*% Z *0.9
        probIa =  2 * pnorm(Zscalar, lower.tail = FALSE)
        #probIa = sqrt(probIa)
        probIa = matrix(probIa, 1,1)
        colnames(probIa) = "Probability"
        rownames(probIa) = "Type Ia"
        return(probIa)
    }, rownames = TRUE, digits = 3)
    
    output$colorPlots <- renderPlot({
        inFile <- input$file1
        
        X = allInfo$R_beta1
        Y = allInfo$BminusV
        Xse = allInfo$R_beta1SD
        Yse = allInfo$BminusVSD
        
        if (is.null(inFile)){
            plotXYSe(X, Y, xlab = "nonlinear R beta_1", ylab = "Color",
                     xyCol = "R",
                     Xse = Xse, Yse = Yse)
            rseq = seq(0.01,0.8, length.out = 100)
            lines(rseq, Renvelope(rseq))
            return(NULL)
        }

        newX = nonlinearParam["nonlinear_beta1", "R"]
        newXse = nonlinearParam["nonlinear_beta1_SD", "R"]
        newY = parameterTable["Mag_at_Bmax", "B"] - 
            parameterTable["Mag_at_Bmax", "V"]
        newYse = parameterTable["mag_at_Bmax_SD", "B"]^2 +
            parameterTable["mag_at_Bmax_SD", "V"]^2
        newYse = sqrt(newYse)
        
        BminusV <<- newY
        IntrinsicC <<- Renvelope(newX)
        EBV <<- newY - Renvelope(newX)
        
        plotXYSe(X, Y, xlab = "nonlinear R beta_1", ylab = "Color",
                 xyCol = "R",
                 Xse = Xse, Yse = Yse, 
                 newX = newX, newY = newY, 
                 newXse = newXse, newYse = newYse,
                 newCols = "R")
        rseq = seq(0.01,0.8, length.out = 100)
        lines(rseq, Renvelope(rseq))
    })
    
    output$colorTable <- renderTable({ 
        inFile <- input$file1
        if (is.null(inFile))  return(NULL)
        
        colTable = matrix(0, 3, 1)
        colTable[1,1] = round(BminusV, 3)
        colTable[2,1] = round(IntrinsicC, 3)
        colTable[3,1] = round(EBV, 3)
        colnames(colTable) = "Color"
        rownames(colTable) = c("Observed", "Intrinsic", "E(B-V)")
        
        return(colTable)
    }, rownames = TRUE, digits = 3)
    
    output$Spc1 <- renderPlot({
        inFile <- input$file1
        if (is.null(inFile)){
            newdata = NULL
        }else{
            newdata = data.frame(
                s1 = parameterTable["beta1", "I"],
                s2 = parameterTable["beta2", "I"],
                s1SD = parameterTable["beta1_SD", "I"], 
                s2SD = parameterTable["beta2_SD", "I"],
                type = "new")
        }
        
        xlab = expression(beta[I]^{(1)})
        ylab = expression(beta[I]^{(2)})
        SNE_plotSpcClass(spcType1, xlab, ylab, newdata)
    })
    
    output$Spc2 <- renderPlot({
        inFile <- input$file1
        if (is.null(inFile)){
            newdata = NULL
        }else{
            newdata = data.frame(
                s1 = parameterTable["beta1", "B"],
                s2 = parameterTable["beta4", "R"],
                s1SD = parameterTable["beta1_SD", "B"], 
                s2SD = parameterTable["beta4_SD", "R"],
                type = "new")
        }
        
        xlab = expression(beta[B]^{(1)})
        ylab = expression(beta[R]^{(4)})
        SNE_plotSpcClass(spcType2, xlab, ylab,newdata)
    })
    
    
    output$Spc3 <- renderPlot({
        inFile <- input$file1
        if (is.null(inFile)){
            newdata = NULL
        }else{
            newdata = data.frame(
                s1 = parameterTable["beta1", "R"],
                s2 = parameterTable["beta3", "R"],
                s1SD = parameterTable["beta1_SD", "R"], 
                s2SD = parameterTable["beta3_SD", "R"],
                type = "new")
        }
        
        xlab = expression(beta[R]^{(1)})
        ylab = expression(beta[R]^{(3)})
        SNE_plotSpcClass(spcType3, xlab, ylab,newdata)
    })
    
    # Call renderPlot for each one. Plots are only actually generated when they
    # are visible on the web page.
    for (i in 1:NLcurves) {
        # Need local so that each item gets its own number. Without it, the value
        # of i in the renderPlot() will be the same across all instances, because
        # of when the expression is evaluated.
        local({
            my_i <- i
            plotname <- paste("plot", my_i, sep="")
            output[[plotname]] <- renderPlot({
                inFile <- input$file1
                if (is.null(inFile))
                    return(NULL)
                cSNe$plotFit("", plotID = my_i)
            })
        })
    }
    
    
}
)

# Run the application 
shinyApp(ui = ui, server = server)

