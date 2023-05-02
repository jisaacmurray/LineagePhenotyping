library(ggplot2)
library(plotly)

library(plot3D)
source("functions.R")
nob1_peak=ReadPeakExpression("data/CA20140825_nob-1_JIM284_L3.csv")

Name="nob-1_mutant"

PlotDeviationsList <- function(Name,exp=nob1_peak,t=200){

    
    WTPositions=read.table("Richard_et_al_plus_comma_WT/Richard_et_al_plus_comma_WTpositions.txt", header=T,stringsAsFactors=F,sep="\t")
    CellTimes=paste(WTPositions[,1],WTPositions[,2],sep=":")
    rownames(WTPositions)=CellTimes

    Cells=unique(WTPositions[,1])
    Times=unique(WTPositions[,2])
    
    wtX=WTPositions[,substr(colnames(WTPositions),start=1,stop=1)=="X"]
    wtEmbryos=colnames(wtX)
    wtEmbryos=sub("X_","",wtEmbryos)
    wtY=WTPositions[,substr(colnames(WTPositions),start=1,stop=1)=="Y"]
    wtZ=WTPositions[,substr(colnames(WTPositions),start=1,stop=1)=="Z"]
    colnames(wtX)=wtEmbryos
    colnames(wtY)=wtEmbryos
    colnames(wtZ)=wtEmbryos
    

    wtXMeans=rowMeans(wtX,na.rm=T)
    wtYMeans=rowMeans(wtY,na.rm=T)
    wtZMeans=rowMeans(wtZ,na.rm=T)

    wtXCellMeans=tapply(wtXMeans,WTPositions[,1],mean,na.rm=T)
    wtYCellMeans=tapply(wtYMeans,WTPositions[,1],mean,na.rm=T)
    wtZCellMeans=tapply(wtZMeans,WTPositions[,1],mean,na.rm=T)


    wtXsds=apply(wtX,1,sd,na.rm=T)
    wtYsds=apply(wtY,1,sd,na.rm=T)
    wtZsds=apply(wtZ,1,sd,na.rm=T)
   
    Counts=rowSums(!is.na(wtX))
    startDir = getwd()
    WTDir=paste(startDir,"Richard_et_al_plus_comma_WT",sep="/")


    MutantPositions=read.table(paste(Name,"/",Name,"positions.txt",sep=""), header=T,stringsAsFactors=F,sep="\t")
    #rownames should be directly comparable to WT with major caveat that this is not aware of potentially major CC changes
    #in mutants (e.g. if cell X divides much later than in WT, its positions are being analyzed relative to embryos that are inappropriately young
    rownames(MutantPositions)=paste(MutantPositions[,1],MutantPositions[,2],sep=":")

    mutantX=MutantPositions[,substr(colnames(MutantPositions),start=1,stop=1)=="X"]
    mutantEmbryos=colnames(mutantX)
    mutantEmbryos=sub("X_","",mutantEmbryos)
    mutantY=MutantPositions[,substr(colnames(MutantPositions),start=1,stop=1)=="Y"]
    mutantZ=MutantPositions[,substr(colnames(MutantPositions),start=1,stop=1)=="Z"]
    colnames(mutantX)=mutantEmbryos
    colnames(mutantY)=mutantEmbryos
    colnames(mutantZ)=mutantEmbryos
    

    mutantXMeans=rowMeans(mutantX,na.rm=T)
    mutantYMeans=rowMeans(mutantY,na.rm=T)
    mutantZMeans=rowMeans(mutantZ,na.rm=T)

    mutantXsds=apply(mutantX,1,sd,na.rm=T)
    mutantYsds=apply(mutantY,1,sd,na.rm=T)
    mutantZsds=apply(mutantZ,1,sd,na.rm=T)


    PlotDeviationsSingle(wtXMeans, wtYMeans, wtZMeans,
                         mutantXMeans,mutantYMeans,mutantZMeans,
                         outfile=paste(Name,"/",Name,"_mean_arrows.pdf",sep=""))
    PlotDeviationsSingle(wtXMeans, wtYMeans, wtZMeans,
                         mutantXMeans,mutantYMeans,mutantZMeans,
                         outfile=paste(Name,"/",Name,"_mean_arrows",sep=""),jpg=T)
    PlotDeviationsSingle(wtXMeans, wtYMeans, wtZMeans,
                         mutantXMeans,mutantYMeans,mutantZMeans,
                         outfile=paste(Name,"/",Name,"_mean_arrows_exp",sep=""),jpg=T,peakExpression=exp)
    PlotDeviationsSingle(wtXMeans, wtYMeans, wtZMeans,
                         mutantXMeans,mutantYMeans,mutantZMeans,
                         outfile=paste(Name,"/",Name,"_mean_arrows_exp.pdf",sep=""),peakExpression=exp)


    ##provide color variable with same names as 'x; e.g. wtXMeans
    Founders = sapply(WTPositions[,1],GetFounder)
    names(Founders) <- names(wtXMeans)
    PlotDeviationsSingle(wtXMeans, wtYMeans, wtZMeans,
                         mutantXMeans,mutantYMeans,mutantZMeans,
                         outfile=paste(Name,"/",Name,"_mean_arrows_lineage.pdf",sep=""),color=Founders,colname="Founders")


    
    for(i in mutantEmbryos){
        theseX=MutantPositions[,paste("X",i,sep="_")]
        theseY=MutantPositions[,paste("Y",i,sep="_")]
        theseZ=MutantPositions[,paste("Z",i,sep="_")]
        names(theseX) <- rownames(MutantPositions)
        names(theseY) <- rownames(MutantPositions)
        names(theseZ) <- rownames(MutantPositions)
        theseX <- theseX[!is.na(theseX)]
        theseY <- theseY[!is.na(theseY)]
        theseZ <- theseZ[!is.na(theseZ)]

        
        PlotDeviationsSingle(wtXMeans, wtYMeans, wtZMeans,
                             theseX,theseY,theseZ,dlim=15,
                             outfile=paste(Name,"/",i,"_arrows.pdf",sep=""),color=Founders,colname="Founders")
    }

    times=WTPositions[,2]
    Devs=data.frame(x=wtXMeans-mutantXMeans[names(wtXMeans)],y=wtYMeans-mutantYMeans[names(wtXMeans)],z=wtZMeans-mutantZMeans[names(wtXMeans)])
    TimeMeanDevs = aggregate(Devs,by=list(times),FUN=mean,na.rm=T)

    pdf(paste(Name,"/",Name,"_posDirectionPlots.pdf",sep=""))
    plot(TimeMeanDevs[,1],sqrt(rowSums(TimeMeanDevs[2:4]^2)),xlab="Time(min)",ylab="Mean Position Bias(microns)")
    plot(TimeMeanDevs[,1],TimeMeanDevs$x,xlab="Time(min)",ylab="Mean Deviation (black:x,red:y,green:z)")
    points(TimeMeanDevs[,1],TimeMeanDevs$y,col=2)
    points(TimeMeanDevs[,1],TimeMeanDevs$z,col=3)
    abline(h=0)
    plot(wtZMeans,Devs$x,xlab="Z Pos",ylab=paste("X Deviation, cor =",round(cor(wtZMeans,Devs$x,use="pairwise.complete.obs"),2)),cex=0.1)
    plot(wtYMeans,Devs$x,xlab="Y Pos",ylab=paste("X Deviation, cor =",round(cor(wtYMeans,Devs$x,use="pairwise.complete.obs"),2)),cex=0.1)
    plot(wtXMeans,Devs$x,xlab="X Pos",ylab=paste("X Deviation, cor =",round(cor(wtXMeans,Devs$x,use="pairwise.complete.obs"),2)),cex=0.1)
    

        
    dev.off()
    #PlotlyDeviations(wtXMeans, wtYMeans, wtZMeans,
    #                      mutantXMeans,mutantYMeans,mutantZMeans,time=250)
    

    return(Devs)
}




PlotDeviationsSingle <- function(x,y,z,mx,my,mz,dlim=10,outfile,jpg=FALSE,peakExpression=NULL,elim=1000,color=NULL,colname=NA,phi=0,theta=0){
    message(paste("Plotting Deviations -", outfile))

    WTPositions=read.table("Richard_et_al_plus_comma_WT/Richard_et_al_plus_comma_WTpositions.txt", header=T,stringsAsFactors=F,sep="\t")
    CellTimes=paste(WTPositions[,1],WTPositions[,2],sep=":")
    rownames(WTPositions)=CellTimes



    commonCells = intersect(names(x),names(mx))
    data=data.frame(x=x[commonCells],y=y[commonCells],z=z[commonCells],mx=mx[commonCells],my=my[commonCells],mz=mz[commonCells],
                    time=WTPositions[commonCells,2],cell=WTPositions[commonCells,1],size=1/WTPositions[commonCells,2])
    data$length <- sqrt((data$x-data$mx)^2+(data$y-data$my)^2+(data$z-data$mz)^2)
    data$length[data$length>dlim] <-  dlim
    if(!is.null(peakExpression)){
        data$exp = pmin(elim,peakExpression[data$cell])
        elim=max(data$exp,na.rm=T)
        emin=min(-300,data$exp,na.rm=T)
    }else if(!is.null(color)){
        data$color=color[commonCells]
    }

    if(jpg){        
        dir.create(outfile)
    }else{
        pdf(outfile)
    }
    for(i in sort(unique(data$time))){
        theseData = data[data$time==i,]
        
        if(jpg){
            jpeg(paste(outfile,"/",i,".jpg",sep=""),width=1024,height=1024)
        }
        if(!is.null(peakExpression)){
            arrows3D(theseData$x,theseData$y,theseData$z,theseData$mx,theseData$my,theseData$mz,
                     xlim=c(-25,25),ylim=c(-25,25),zlim=c(-25,25),
                     xlab="AP",ylab="LR",zlab="DV",
                     colvar=theseData$exp,clim=c(emin,elim), clab="Expression",
                     phi=phi,theta=theta,labels=theseData$cell,
                     main=paste("Time ",i, " minutes", sep=""),
                     lwd=1,length=.1)
            
        }else if(!is.null(color)){
            par(mar=c(1,1,1,1))
            colorScheme = rainbow(length(levels(as.factor(theseData$color))))
            names(colorScheme)=levels(as.factor(theseData$color))
            arrows3D(theseData$x,theseData$y,theseData$z,theseData$mx,theseData$my,theseData$mz,
                     xlim=c(-25,25),ylim=c(-25,25),zlim=c(-25,25),
                     xlab="AP",ylab="LR",zlab="DV",
                     colvar=as.numeric(as.factor(theseData$color)),col=colorScheme,colkey=FALSE,
                     phi=phi,theta=theta,labels=theseData$cell,
                     main=paste("Time ",i, " minutes", sep=""),
                     lwd=1,length=.1
                     )
            legend(.3,.05,names(colorScheme),col=colorScheme,pch=1)
            
        }else{
            arrows3D(theseData$x,theseData$y,theseData$z,theseData$mx,theseData$my,theseData$mz,
                     xlim=c(-25,25),ylim=c(-25,25),zlim=c(-25,25),
                     xlab="AP",ylab="LR",zlab="DV",
                     colvar=theseData$length,clim=c(0,dlim), clab="Dev(microns)",
                     phi=phi,theta=theta,labels=theseData$cell,
                     main=paste("Time ",i, " minutes", sep=""),
                     lwd=1,length=.1)
        }
        
        if(jpg){
            dev.off()
        }
    }
    if(!jpg){dev.off()}
}
    
PlotlyDeviations <- function(x,y,z,mx,my,mz,time){
    commonCells = intersect(names(x),names(mx))
    data=data.frame(x=x[commonCells],y=y[commonCells],z=z[commonCells],mx=mx[commonCells],my=my[commonCells],mz=mz[commonCells],
                    time=WTPositions[commonCells,2],cell=WTPositions[commonCells,1],size=1/WTPositions[commonCells,2])
    theseData = data[data$time==time,]
    
    ##could allow plotting only a subset by dev length threshold
##    data$length <- sqrt((data$x-data$mx)^2+(data$y-data$my)^2+(data$z-data$mz)^2)
##    data$length[data$length>dlim] <-  dlim
    segments = NULL
    for(cell in theseData$cell){
        cellData = theseData[theseData$cell==cell,]
        
        deviation=sqrt((cellData$x-cellData$mx)^2 +
                       (cellData$y-cellData$my)^2 +
                       (cellData$z-cellData$mz)^2)
            
             thisLine = data.frame(x=c(cellData$x,cellData$mx),y=c(cellData$y,cellData$my),z=c(cellData$z,cellData$mz),cellName=cell,dev=c(deviation,NA),exp=c(cellData$exp,NA),col=c("black","white"))
             segments=rbind(segments,thisLine)
         }
    fig <- plot_ly(data=segments, x = ~x, y = ~y, z = ~z,  mode = 'lines',
                   line = list(width = 6,reverscale = FALSE),
                   color = ~cellName, text=~cellName, group=~cellName, hoverinfo="dev")    
    fig
}



GetFounder <- function(cell){
    if(grepl("ABala",cell)){
        return("ABala")
    }else if(grepl("ABalp",cell)){
        return("ABalp")
    }else if(grepl("ABara",cell)){
        return("ABara")
    }else if(grepl("ABarp",cell)){
        return("ABarp")
    }else if(grepl("ABpla",cell) | grepl("ABpra",cell)){
        return("ABpxa")
    }else if(grepl("ABplp",cell) | grepl("ABprp",cell)){
        return("ABpxp")
    }else if(grepl("MSaa",cell) | grepl("MSpa",cell)){
        return("MSxa")
    }else if(grepl("MSap",cell) | grepl("MSpp",cell)){
        return("MSxp")
    }else if(grepl("Ea",cell)){
        return("Ea")
    }else if(grepl("Ep",cell)){
        return("Ep") 
   }else if(grepl("Caa",cell) | grepl("Cpa",cell)){
        return("Cxa")
    }else if(grepl("Cap",cell) | grepl("Cpp",cell)){
        return("Cxp")
    }else if(grepl("D",cell)){
        return("D")
    }else{
        return("P")
    }
}
     



     




##     data <- read.csv('https://raw.githubusercontent.com/plotly/datasets/master/3d-line1.csv')
##     data$color <- as.factor(data$color)

##     data=data.frame(x=wtXMeans[commonCells],y=wtYMeans[commonCells],z=wtZMeans[commonCells],
##                     mx=mutantXMeans[commonCells],my=mutantYMeans[commonCells],mz=mutantZMeans[commonCells],
##                     time=WTPositions[commonCells,2],cell=WTPositions[commonCells,1],
##                     exp=exp[WTPositions[commonCells,1]],size=1/WTPositions[commonCells,2])
##     data$exp[is.na(data$exp)]<-0

##     data$length <- sqrt((data$x-data$mx)^2+(data$y-data$my)^2+(data$z-data$mz)^2)
##     pdf(paste(Name,"/",Name,"_mean_arrows.pdf",sep=""))
##     for(i in sort(unique(data$time))){
##         message(i)
##         theseData = data[data$time==i,]

##         arrows3D(theseData$x,theseData$y,theseData$z,theseData$mx,theseData$my,theseData$mz,xlim=c(-25,25),ylim=c(-25,25),zlim=c(-25,25),xlab="AP",ylab="LR",zlab="DV",colvar=theseData$length,phi=0,theta=0,clim=c(0,10),clab="Mean Deviation(microns)",main=paste("Time ",i, " minutes", sep=""))
##     }
##     dev.off()
## #         arrows3D (x0, y0, z0, x1 = x0, y1 = y0, z1 = z0, ...,  
## #              colvar = NULL, phi = 40, theta = 40,
## #              col = NULL, NAcol = "white", breaks = NULL,
## #              colkey = NULL, panel.first = NULL,
## #              clim = NULL, clab = NULL, bty = "b", type = "triangle", 
## #              add = FALSE, plot = TRUE)
    
        
##         ##simple plot of cell positions (WT) in plotly
##         fig <- plot_ly(theseData,x=~x, y=~y, z=~z, type="scatter3d", mode="markers", color=~exp, size=I(1000000/i))
##         fig


##         ##one way to plot segments in plotly - make a series of segment ends in a data frame, use grouping by cell to separate by color
##         segments = NULL
##         for(cell in theseData$cell){
##             cellData = theseData[theseData$cell==cell,]

##             deviation=sqrt((cellData$x-cellData$mx)^2 +
##                  (cellData$y-cellData$my)^2 +
##                  (cellData$z-cellData$mz)^2)
            
##             thisLine = data.frame(x=c(cellData$x,cellData$mx),y=c(cellData$y,cellData$my),z=c(cellData$z,cellData$mz),cellName=cell,dev=c(deviation,NA),exp=c(cellData$exp,NA),col=c("black","white"))
##             segments=rbind(segments,thisLine)
##         }
##         fig <- plot_ly(data=segments, x = ~x, y = ~y, z = ~z,  mode = 'lines',
##                        line = list(width = 6,reverscale = FALSE),
##                        color = ~cellName, text=~cellName, group=~cellName, hoverinfo="dev")

##         fig <- fig %>% add_trace(segments,x=~x, y=~y, z=~z, type="scatter3d", mode="markers", size=I(1000000/i))

        
##         fig

## color=~col
##         type = 'scatter3d',
##         plot(fig)
    
       

    
    
 
## }
    
##     length(mutantXMeans[commonCells])
##     length(wtXMeans[commonCells])
## }




## ##sample code draws a happy face
## line1 <- data.frame(x=seq(3.5,4.5,len=NP), y=rep(2.5,NP), text="hello")
## line2 <- data.frame(x=seq(3,3.5,len=NP), y=seq(3,2.5,len=NP), text="mouth")
## line3 <- data.frame(x=seq(4.5,5,len=NP), y=seq(2.5,3,len=NP), text="mouth")
## line4 <- data.frame(x=rep(4,NP), y=seq(2.75,3.5,len=NP), text="nose")
## rect1 <- data.frame(x=c(seq(2,6,len=NP), rep(6,NP), seq(6,2,len=NP), rep(2,NP)),
##                     y=c(rep(2,NP), seq(2,4.5,len=NP), rep(4.5,NP), seq(4.5,2,len=NP)),
##                     text="head")
## rect2 <- data.frame(x=c(seq(2.5,3.5,len=NP), rep(3.5,NP), seq(3.5,2.5,len=NP), rep(2.5,NP)),
##                     y=c(rep(3.5,NP), seq(3.5,4,len=NP), rep(4,NP), seq(4,3.5,len=NP)),
##                     text="left eye")
## rect3 <- data.frame(x=c(seq(4.5,5.5,len=NP), rep(5.5,NP), seq(5.5,4.5,len=NP), rep(4.5,NP)),
##                     y=c(rep(3.5,NP), seq(3.5,4,len=NP), rep(4,NP), seq(4,3.5,len=NP)),
##                     text="right eye")

## trace_dat <- rbind(line1, line2, line3, line4, rect1, rect2, rect3)

## plot_ly(data=trace_dat, x=~x, y=~y, mode="lines", hoverinfo="text", text=~text, group = ~text,color=~text)

