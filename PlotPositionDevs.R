library(ggplot2)
library(gridExtra)
library(tidyr)
source("functions.R")

PlotDeviationsList <- function(Name, exp=NULL, t=200, skipIndividual=FALSE, skipArrows=FALSE,
                               data_dir=".", output_dir=NULL, wt_ref_dir=NULL){
    if (is.null(output_dir)) output_dir <- file.path(data_dir, Name)
    if (is.null(wt_ref_dir)) wt_ref_dir <- file.path(data_dir, "Richard_et_al_plus_comma_WT")
    wt_prefix <- basename(wt_ref_dir)
    mutant_input_dir <- file.path(data_dir, Name)

    message(paste("Plotting Deviations List for", Name))

    WTPositions=read.table(file.path(wt_ref_dir, paste0(wt_prefix, "positions.txt")), header=T,stringsAsFactors=F,sep="\t", check.names=FALSE)
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
    WTDir = wt_ref_dir


    # Load Rotated Mutant Data (written by AnalyzePositions to output_dir)
    message("Loading rotated mutant position data...")
    mutantX = read.csv(file.path(output_dir, paste0(Name, "_rotatedX.csv")), row.names=1, check.names=FALSE)
    mutantY = read.csv(file.path(output_dir, paste0(Name, "_rotatedY.csv")), row.names=1, check.names=FALSE)
    mutantZ = read.csv(file.path(output_dir, paste0(Name, "_rotatedZ.csv")), row.names=1, check.names=FALSE)
    
    mutantEmbryos = colnames(mutantX)
    # The columns in rotated csvs already have the correct names (might have X prefix if check.names was T, 
    # but we used FALSE in AnalyzePositions writing? No, AnalyzePositions didn't specify. 
    # Let's be robust.
    colnames(mutantX) = sub("^X", "", colnames(mutantX))
    colnames(mutantY) = sub("^X", "", colnames(mutantY))
    colnames(mutantZ) = sub("^X", "", colnames(mutantZ))
    mutantEmbryos = colnames(mutantX)
    

    mutantXMeans=rowMeans(mutantX,na.rm=T)
    mutantYMeans=rowMeans(mutantY,na.rm=T)
    mutantZMeans=rowMeans(mutantZ,na.rm=T)

    mutantXsds=apply(mutantX,1,sd,na.rm=T)
    mutantYsds=apply(mutantY,1,sd,na.rm=T)
    mutantZsds=apply(mutantZ,1,sd,na.rm=T)


    if(!skipArrows){
        arrows_dir <- .plots_dir(output_dir, "position/arrows")
        PlotDeviationsSingle(wtXMeans, wtYMeans, wtZMeans,
                             mutantXMeans,mutantYMeans,mutantZMeans,
                             outfile=file.path(arrows_dir, paste0(Name,"_mean_arrows.pdf")),
                             wt_ref_dir=wt_ref_dir)
        PlotDeviationsSingle(wtXMeans, wtYMeans, wtZMeans,
                             mutantXMeans,mutantYMeans,mutantZMeans,
                             outfile=file.path(arrows_dir, paste0(Name,"_mean_arrows")),jpg=T,
                             wt_ref_dir=wt_ref_dir)
        PlotDeviationsSingle(wtXMeans, wtYMeans, wtZMeans,
                             mutantXMeans,mutantYMeans,mutantZMeans,
                             outfile=file.path(arrows_dir, paste0(Name,"_mean_arrows_exp")),jpg=T,peakExpression=exp,
                             wt_ref_dir=wt_ref_dir)
        PlotDeviationsSingle(wtXMeans, wtYMeans, wtZMeans,
                             mutantXMeans,mutantYMeans,mutantZMeans,
                             outfile=file.path(arrows_dir, paste0(Name,"_mean_arrows_exp.pdf")),peakExpression=exp,
                             wt_ref_dir=wt_ref_dir)


        ##provide color variable with same names as 'x; e.g. wtXMeans
        Founders = sapply(WTPositions[,1],GetFounder)
        names(Founders) <- names(wtXMeans)
        PlotDeviationsSingle(wtXMeans, wtYMeans, wtZMeans,
                             mutantXMeans,mutantYMeans,mutantZMeans,
                             outfile=file.path(arrows_dir, paste0(Name,"_mean_arrows_lineage.pdf")),color=Founders,colname="Founders",
                             wt_ref_dir=wt_ref_dir)

        # Phase 5.2D: clean up the per-frame jpg directories now that the
        # corresponding PDFs are sealed. Set keep_arrow_jpgs=TRUE in
        # YAML/options to preserve them (e.g. for movie generation).
        if (!isTRUE(getOption("LineagePhenotyping.keep_arrow_jpgs", FALSE))) {
            for (sub in c("_mean_arrows", "_mean_arrows_exp")) {
                d <- file.path(arrows_dir, paste0(Name, sub))
                if (dir.exists(d)) unlink(d, recursive = TRUE)
            }
        }
    }


    
    if(!skipIndividual){
        for(i in mutantEmbryos){
            theseX=mutantX[,i]
            theseY=mutantY[,i]
            theseZ=mutantZ[,i]
            names(theseX) <- rownames(mutantX)
            names(theseY) <- rownames(mutantX)
            names(theseZ) <- rownames(mutantX)
            theseX <- theseX[!is.na(theseX)]
            theseY <- theseY[!is.na(theseY)]
            theseZ <- theseZ[!is.na(theseZ)]

            
            PlotDeviationsSingle(wtXMeans, wtYMeans, wtZMeans,
                                 theseX,theseY,theseZ,dlim=15,
                                 outfile=file.path(.plots_dir(output_dir, "position/arrows", per_embryo = TRUE), paste0(i,"_arrows.pdf")),color=Founders,colname="Founders",
                                 wt_ref_dir=wt_ref_dir)
        }
    }

    time=WTPositions[,2]
    Devs=data.frame(x=wtXMeans-mutantXMeans[names(wtXMeans)],y=wtYMeans-mutantYMeans[names(wtXMeans)],z=wtZMeans-mutantZMeans[names(wtXMeans)])
    TimeMeanDevs = aggregate(Devs,by=list(Time=time),FUN=mean,na.rm=T)

    pdf(file.path(.plots_dir(output_dir, "position/cell"), paste0(Name,"_posDirectionPlots.pdf")), width=10, height=8)
    
    # 1. Total Mean Bias vs Time
    TimeNumeric <- as.numeric(TimeMeanDevs[,1])
    TotalBias <- sqrt(rowSums(TimeMeanDevs[2:4]^2))
    p_bias <- ggplot(data.frame(Time=TimeNumeric, Bias=TotalBias), aes(x=Time, y=Bias)) +
        geom_line() + geom_point(size=0.5) + theme_bw() + 
        labs(x="Time (min)", y="Mean Position Bias (µm)", title="Total Position Bias over Time")
        
    # 2. XYZ Mean Deviation vs Time
    time_dev_long <- TimeMeanDevs %>%
        gather(key="Axis", value="Deviation", x, y, z)
    time_dev_long$Time <- as.numeric(time_dev_long$Time)
    
    p_time_xyz <- ggplot(time_dev_long, aes(x=Time, y=Deviation, color=Axis)) +
        geom_line() + geom_point(size=0.5) + 
        geom_hline(yintercept=0, linetype="dashed") +
        scale_color_manual(values=c("x"="black", "y"="red", "z"="green")) +
        theme_bw() +
        labs(x="Time (min)", y="Mean Deviation (µm)", title="Mean Deviation per Axis")
        
    # 3. Correlations (X Deviation vs WT Position)
    cor_df <- data.frame(
        WT_X = wtXMeans,
        WT_Y = wtYMeans,
        WT_Z = wtZMeans,
        Dev_X = Devs$x
    )
    
    get_cor_label <- function(x,y) paste("R =", round(cor(x, y, use="pairwise.complete.obs"), 2))
    
    # Use geom_bin2d for better legibility of large datasets
    p_xz <- ggplot(cor_df, aes(x=WT_Z, y=Dev_X)) + 
        geom_bin2d(bins=100) + scale_fill_gradient(low="lightgrey", high="blue") + theme_bw() +
        labs(x="WT Z Position", y="X Deviation", title=paste("Z vs X Dev,", get_cor_label(cor_df$WT_Z, cor_df$Dev_X)))
        
    p_xy <- ggplot(cor_df, aes(x=WT_Y, y=Dev_X)) + 
        geom_bin2d(bins=100) + scale_fill_gradient(low="lightgrey", high="blue") + theme_bw() +
        labs(x="WT Y Position", y="X Deviation", title=paste("Y vs X Dev,", get_cor_label(cor_df$WT_Y, cor_df$Dev_X)))
        
    p_xx <- ggplot(cor_df, aes(x=WT_X, y=Dev_X)) + 
        geom_bin2d(bins=100) + scale_fill_gradient(low="lightgrey", high="blue") + theme_bw() +
        labs(x="WT X Position", y="X Deviation", title=paste("X vs X Dev,", get_cor_label(cor_df$WT_X, cor_df$Dev_X)))
        
    grid.arrange(p_bias, p_time_xyz, nrow=2)
    grid.arrange(p_xx, p_xy, p_xz, nrow=2, top="Correlations: WT Position vs Mutant X Deviation")
        
    dev.off()
    #PlotlyDeviations(wtXMeans, wtYMeans, wtZMeans,
    #                      mutantXMeans,mutantYMeans,mutantZMeans,time=250)
    

    return(Devs)
}




PlotDeviationsSingle <- function(x,y,z,mx,my,mz,dlim=10,outfile,jpg=FALSE,peakExpression=NULL,elim=1000,color=NULL,colname=NA,phi=0,theta=0,
                                 wt_ref_dir="Richard_et_al_plus_comma_WT"){
    message(paste("Plotting Deviations -", outfile))

    wt_prefix <- basename(wt_ref_dir)
    WTPositions=read.table(file.path(wt_ref_dir, paste0(wt_prefix, "positions.txt")), header=T,stringsAsFactors=F,sep="\t")
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
        emin=min(data$exp,na.rm=T) # Fixed min for expression scaling
    }else if(!is.null(color)){
        data$color=color[commonCells]
    }
    data <- data[complete.cases(data[,c("x","y","z","mx","my","mz")]),]

    # Phase 5.3G: Compute global axis ranges across ALL timepoints up front.
    # Each per-frame ggplot used to auto-fit its xlim/ylim to that frame's
    # cells, which made the embryo visually grow/shrink as the cell count
    # changed across time and obscured real movement. Now every frame uses
    # the same fixed limits, so motion through a stable frame is visible.
    # Includes both WT and mutant endpoints so the arrow tips never clip.
    .pad_range <- function(rng, pad_frac = 0.05) {
        if (any(!is.finite(rng))) return(c(-1, 1))
        span <- diff(rng); if (span == 0) span <- 1
        c(rng[1] - pad_frac * span, rng[2] + pad_frac * span)
    }
    xlim_g <- .pad_range(range(c(data$x, data$mx), na.rm=TRUE))
    ylim_g <- .pad_range(range(c(data$y, data$my), na.rm=TRUE))
    zlim_g <- .pad_range(range(c(data$z, data$mz), na.rm=TRUE))

    # Setup plot output
    if(jpg){
        if(!dir.exists(outfile)) dir.create(outfile)
    }else{
        pdf(outfile, width=12, height=6)
    }

    times <- sort(unique(data$time))

    for(i in times){
        theseData = data[data$time==i,]

        if(jpg){
             jpeg(file.path(outfile, paste0(i,".jpg")),width=1200,height=600)
        }

        # Base plot function
        create_proj <- function(data, x_col, y_col, x_lab, y_lab, mx_col, my_col,
                                 xlim_panel, ylim_panel) {
            p <- ggplot(data) +
                theme_bw() +
                labs(x=x_lab, y=y_lab, title=paste(x_lab, "vs", y_lab)) +
                coord_fixed(xlim = xlim_panel, ylim = ylim_panel)

            if(!is.null(peakExpression)){
                p <- p + geom_segment(aes_string(x=x_col, y=y_col, xend=mx_col, yend=my_col, color="exp"), arrow=arrow(length=unit(0.1,"cm"))) +
                         geom_point(aes_string(x=x_col, y=y_col, color="exp")) +
                         scale_color_gradient(low="blue", high="red", limits=c(emin, elim))
            } else if(!is.null(color)){
                p <- p + geom_segment(aes_string(x=x_col, y=y_col, xend=mx_col, yend=my_col, color="color"), arrow=arrow(length=unit(0.1,"cm"))) +
                         geom_point(aes_string(x=x_col, y=y_col, color="color"))
            } else {
                 p <- p + geom_segment(aes_string(x=x_col, y=y_col, xend=mx_col, yend=my_col, color="length"), arrow=arrow(length=unit(0.1,"cm"))) +
                          geom_point(aes_string(x=x_col, y=y_col, color="length")) +
                          scale_color_viridis_c(limits=c(0, dlim), name="Dev (µm)")
            }
            return(p)
        }

        p_xy <- create_proj(theseData, "x", "y", "AP (x)", "LR (y)", "mx", "my",
                             xlim_g, ylim_g)
        p_xz <- create_proj(theseData, "x", "z", "AP (x)", "DV (z)", "mx", "mz",
                             xlim_g, zlim_g)
        p_yz <- create_proj(theseData, "y", "z", "LR (y)", "DV (z)", "my", "mz",
                             ylim_g, zlim_g)

        grid.arrange(p_xy, p_xz, p_yz, nrow=1, top=paste("Time", i, "min -", nrow(theseData), "cells"))

        if(jpg){
            dev.off()
        }
    }
    if(!jpg){dev.off()}
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

