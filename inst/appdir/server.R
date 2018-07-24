function(input, output) {
    output$ProbSuccesses <- renderPlot({
      p0=input$p0
      p1=input$p1
      pF=input$pF
     shape1F=input$shape1F
      shape2F=input$shape2F
      n.range_1=as.numeric(input$n.range_1)

      par(mfrow=c(1,2))
      plotBDP2(x="n",y="Prob0Successes",n=c(n.range_1[1],n.range_1[2]),p0=p0,p1=p1)
      plotBDP2(x="n",y="PostProb0or1Successes",n=c(n.range_1[1],n.range_1[2]),pF=pF,shape1F=shape1F,shape2F=shape2F)
    })
    
    
    cE.vs.PEcall= reactive({
      pF=input$pF
      pE=input$pE
      p0=input$p0
      p1=input$p1
      shape1F=input$shape1F
      shape2F=input$shape2F
      shape1E=input$shape1E
      shape2E=input$shape2E

      cF=input$cF
      cE.range_1=as.numeric(input$cE.range_1)
      cEaccuracy=0.001
      cEvec=seq(from=cE.range_1[1],to=cE.range_1[2],by=cEaccuracy)
      interims.at= c(input$firstInterim.at,as.numeric(unlist(strsplit(input$furtherInterims.at," "))))
#     vn.int=input$firstInterim.at
      n=input$nfinal

      res=plotBDP2(x="cE",y="PEcall",n=n,interim.at=interims.at,pF=pF,cF=cF,pE=pE,cE=cEvec,p0=p0,p1=p1,
                   shape1F=shape1F,shape2F=shape2F,shape1E=shape1E,shape2E=shape2E,col=c("green","red"),cex.lab=1.4,show=FALSE)
      return(res)
  
    })
    
    
    output$cE.vs.PEcall <- renderPlot({
      plot.cE_vs_pEcall(cE.vs.PEcall())
      abline(v=input$cE,col="gray",lty="dashed")
    })

    output$OCs.selected.cE <- renderText({
      cE=input$cE
      res=cE.vs.PEcall()
      paste0("Selected cE leads to a type I error of ", round(res$y.p0[which.min(abs(res$x.p0-cE))[1]],3)," at p0=",input$p0,
             " and a power of ",round(res$y.p1[which.min(abs(res$x.p1-cE))[1]],3)," at p1=",input$p1,".")

    })


    output$n.vs.bFbE <- renderPlot({
      pF=input$pF
      pE=input$pE
      shape1F=input$shape1F
      shape2F=input$shape2F
      shape1E=input$shape1E
      shape2E=input$shape2E

      cF=input$cF
      cE=input$cE
      vn.int=input$firstInterim.at
      n=input$n.range_4[2]
      plotBDP2(x="n",y="bFbE",n=n,pF=pF,cF=cF,pE=pE,cE=cE,
                shape1F=shape1F,shape2F=shape2F,shape1E=shape1E,shape2E=shape2E,col=c("red","green"))
    })


    output$n.vs.PFstopEcall <- renderPlot({
      pF=input$pF
      pE=input$pE
      p0=input$p0
      p1=input$p1
      cF=input$cF
      cE=input$cE
      shape1F=input$shape1F
      shape2F=input$shape2F
      shape1E=input$shape1E
      shape2E=input$shape2E
      n.range_2=as.numeric(input$n.range_2)
        nvec=c(n.range_2[1]:n.range_2[2])
      interims.at= c(input$firstInterim.at,as.numeric(unlist(strsplit(input$furtherInterims.at," "))))

      par(mfrow=c(1,2))
      plotBDP2(x="n",y="PFstopEcall",n=nvec,interim.at=interims.at,pF=pF,cF=cF,pE=pE,cE=cE,ptrue=p0,
                         shape1F=shape1F,shape2F=shape2F,shape1E=shape1E,shape2E=shape2E, progress=T,
               main="Type I error and probability of true stopping for varying n",cex.lab=1.4)
      plotBDP2(x="n",y="PFstopEcall",n=nvec,interim.at=interims.at,pF=pF,cF=cF,pE=pE,cE=cE,ptrue=p1,
                         shape1F=shape1F,shape2F=shape2F,shape1E=shape1E,shape2E=shape2E,progress=T,
               main="Power and probability of false stopping for varying n",cex.lab=1.4,cex.legend = 1)
    })


    output$ptrue.vs.PEcall <- renderPlot({
      pF=input$pF
      pE=input$pE
      p0=input$p0
      p1=input$p1
      cF=input$cF
      cE=input$cE
      shape1F=input$shape1F
      shape2F=input$shape2F
      shape1E=input$shape1E
      shape2E=input$shape2E
      nvec= as.numeric(unlist(strsplit(input$nfinal.vec," ")))


      pvec=seq(input$ptrue.range_1[1],input$ptrue.range_1[2],by=.01)

      interims.at= c(input$firstInterim.at,as.numeric(unlist(strsplit(input$furtherInterims.at," "))))

      par(mfrow=c(1,2))
      # x="ptrue",y="PEcall"
        n=nvec[1]
        # vn.int=seq(0,n,by=input$interim.atEvery)[-1]
        # vn.int=vn.int[-length(vn.int)]
       vn.int=interims.at[interims.at<n]
        plotBDP2(x="ptrue",y="PEcall",n=n,interim.at=vn.int,pF=pF,cF=cF,pE=pE,cE=cE,ptrue=pvec,
                          shape1F=shape1F,shape2F=shape2F,shape1E=shape1E,shape2E=shape2E,col=1,cex.lab=1.4)
        abline(v=p0,col="grey")
        abline(v=p1,col="grey")

        for (jj in 2:length(nvec)) {
          n=nvec[jj]
          # vn.int=seq(0,n,by=input$interim.atEvery)[-1]
          # vn.int=vn.int[-length(vn.int)]
        vn.int=interims.at[interims.at<n]
         plotBDP2(x="ptrue",y="PEcall",n=n,interim.at=vn.int,pF=pF,cF=cF,pE=pE,cE=cE,ptrue=pvec,
                            shape1F=shape1F,shape2F=shape2F,shape1E=shape1E,shape2E=shape2E,add=TRUE,col=jj)

        }
        legend("bottomright",title="final analysis at ",legend=nvec,text.col = 1:length(nvec))

      # (x="ptrue",y="ExpectedNumber"
         n=nvec[1]
        # vn.int=seq(0,n,by=input$interim.atEvery)[-1]
        # vn.int=vn.int[-length(vn.int)]
        vn.int=interims.at[interims.at<n]
        plotBDP2(x="ptrue",y="ExpectedNumber",n=n,interim.at=vn.int,pF=pF,cF=cF,pE=pE,cE=cE,ptrue=pvec,shape1F=shape1F,shape2F=shape2F,col=1,ylim=c(0,max(nvec)),cex.lab=1.4)
        abline(v=p0,col="grey")
        abline(v=p1,col="grey")

        for (jj in 2:length(nvec)) {
          n=nvec[jj]
          # vn.int=seq(0,n,by=input$interim.atEvery)[-1]
          # vn.int=vn.int[-length(vn.int)]
          vn.int=interims.at[interims.at<n]
          plotBDP2(x="ptrue",y="ExpectedNumber",n=n,interim.at=vn.int,pF=pF,cF=cF,pE=pE,cE=cE,ptrue=pvec,shape1F=shape1F,shape2F=shape2F,col=jj,add=TRUE)
        }
        legend("bottomright",title="final analysis at ",legend=nvec,text.col = 1:length(nvec))
    })


  output$report <- downloadHandler(
      # For PDF output, change this to "report.pdf"
       filename = function() {
      paste('report', sep = '.', switch(
        input$format, PDF = 'pdf', HTML = 'html', Word = 'docx'
      ))
},
      content = function(file) {
        # Copy the report file to a temporary directory before processing it, in
        # case we don't have write permissions to the current working dir (which
        # can happen when deployed).
        tempReport <- file.path(tempdir(), "report.Rmd")
      #  file.copy(paste0(normalizePath( find.package("BDP2"),winslash = "/"),"/exdata/report.Rmd"), tempReport, overwrite = TRUE)
        file.copy(file.path(system.file("appdir", package = "BDP2"), "report.Rmd"), tempReport, overwrite = TRUE)

        # Set up parameters to pass to Rmd document
        params <- list(
                        pF=input$pF,
                        pE=input$pE,
                        p0=input$p0,
                        p1=input$p1,
                        shape1F=input$shape1F,
                        shape2F=input$shape2F,
                        shape1E=input$shape1E,
                        shape2E=input$shape2E,
                        n.range_1=input$n.range_1,
                        n.range_2=input$n.range_2,
                        n.range_4=input$n.range_4,
                        cF=input$cF,
                        cE.range_1=input$cE.range_1,
                        cE=input$cE,
                        firstInterim.at=input$firstInterim.at,
                        furtherInterims.at=input$furtherInterims.at,
                        nfinal=input$nfinal,
                        nfinal.vec= input$nfinal.vec,
                        ptrue.range_1=input$ptrue.range_1
                      )

        # Knit the document, passing in the `params` list, and eval it in a
        # child of the global environment (this isolates the code in the document
        # from the code in this app).
        # render(tempReport, output_file = file,
        #   params = params,
        #   envir = new.env(parent = globalenv())
        # )
        out <- render(tempReport, switch(
        input$format,
        PDF = pdf_document(), HTML = html_document(), Word = word_document()
      ),params = params,
          envir = new.env(parent = globalenv())
      )
      file.rename(out, file)
      }
    )

}




