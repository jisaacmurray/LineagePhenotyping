source("functions.R")
library(rgl)
library(gdata)
library(plotrix)

library(ggplot2)


source("GetPeak.R")

##Note that while the default is to compare position defects with ceh-36 expression, the function allows expression peak values for any gene to be passed in

ceh36_exp=read.csv("data/CA20120714_JIM136_L3.csv",row.names=2)["blot"]
ceh36_peak=sapply(Cells,function(X){
        GetPeak(ceh36_exp,X)
        })  


tmpCellNames=read.csv("CellNames.csv")
CellNames=as.vector(tmpCellNames[,2])
names(CellNames)=tmpCellNames[,1]

write.csv(ceh36_peak,file="ceh36_peak.csv")

PlotExpVsDev <- function(Name,outfile,exp=ceh36_peak, type="mean", eGain=1000,ylim=c(0,2),xlim=c(0,10)){ 

    pdf(paste(Name,outfile,sep="/"))

    devs = read.csv(paste(Name,"/",Name, "_Cell", type, "PositionDevs.csv", sep=""),row.names=1)
    devs[,1] <- NULL
    NNs =  read.csv(paste(Name,"/",Name, "_NN_Scores.csv", sep=""),row.names=1)
    NNs[,1] <-  NULL

    meanDevs =rowMeans(devs,na.rm=T)
    meanNNs = rowMeans(NNs,na.rm=T)

    data = data.frame(dev = meanDevs, NN = meanNNs[names(meanDevs)], e = pmax(-200,pmin(exp[names(meanDevs)],eGain)))
    
    p <- ggplot(data, aes(x=dev,y=NN)) +
        geom_point(aes(color = e), size = 1) +
        theme_classic() +        
        scale_color_gradient(name="Expression",low="blue", high="red")+
        ggtitle(paste(Name,"type"),"Mean of Embryos")+
        labs(x="Mean Deviation (microns)", y="Mean Neighbor Deviation Score") +
        xlim(xlim)+
        ylim(ylim)
    
    plot(p)
##        scale_colour_gradientn(name="Expression",colours = viridisLite::viridis(100)) +

    for(i in colnames(devs)){
        
        data = data.frame(dev = devs[names(meanDevs),i], NN = NNs[names(meanDevs),i], e = pmin(exp[names(meanDevs)],eGain))
        p <- ggplot(data, aes(x=dev,y=NN)) +
            geom_point(aes(color = e), size = 1) +
            theme_classic() +
            scale_color_gradient(name="Expression",low="blue", high="red")+
            ggtitle(i) +
            labs(x="Deviation (microns)", y="Mean Neighbor Deviation Score")+
            xlim(xlim)+
            ylim(ylim)

        plot(p)
    }

    dev.off()
}


#START POSITIONS
AnalyzePositions <- function(Name, CalculateNeighbors=FALSE, Expression=ceh36_peak,PercentCellsRequired=0.25, expCutoff=500){
    
   message(paste("Analyzing Positions -" , Name, sep=" "))

    ##should be deprecated but leaving in case I missed renaming anywhere
    ceh36_peak=Expression

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

#    message(Name)
    setwd(Name)

    MutantPositions=read.table(paste(Name,"positions.txt",sep=""), header=T,stringsAsFactors=F,sep="\t")
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

    #Compute mean/stdev statistics for WT
    WTXdevs=wtX-wtXMeans
    wtX_SD=apply(WTXdevs,1,sd,na.rm=T)
    WTYdevs=wtY-wtYMeans
    wtY_SD=apply(WTYdevs,1,sd,na.rm=T)
    WTZdevs=wtZ-wtZMeans
    wtZ_SD=apply(WTZdevs,1,sd,na.rm=T)
    
    WTDdevs=sqrt(WTXdevs^2+WTYdevs^2+WTZdevs^2)
    WTDmeans=apply(WTDdevs,1,mean,na.rm=T)
    WTDsds=apply(WTDdevs,1,sd,na.rm=T)
    
    WTDcounts=apply(WTDdevs,1,function(X){sum(!is.na(X))})
    pdf(paste(Name,"_WT_stats.pdf",sep=""))
    plot(WTDmeans,WTDsds, xlab="WT mean deviation (microns)", ylab="STDEV of WT deviations (microns)", cex=0.1)
    plot(WTDsds,jitter(WTDcounts), xlab="STDEV of WT deviations (microns)", ylab="WT observed counts", cex=0.1)
    boxplot(WTDsds~WTDcounts,xlab="WT Observed Counts",ylab="STDEV of WT deviations (microns)")
    boxplot(WTDmeans~WTDcounts,xlab="WT Observed Counts",ylab="Mean WT deviation (microns)")
    dev.off()
    WTDsds=(WTDsds[WTDcounts>2])[names(WTDcounts)]


    message("Calculating neighbors")
    ##for debugging purposes
    ##    CalculateNeighbors=FALSE
    if(CalculateNeighbors == TRUE){
     NN_distance_summaries=data.frame(Cells)
#     Find defects embryo by embryo - UNCOMMENT TO PERFORM - ~30 minutes per embryo?
      for(i in mutantEmbryos){
        #Try single deviation calculation.  Various options - integrate over all times, try targeted time points per cell, use average positions
         #for each Cell
        CellDeviationScores=sapply(Cells,function(X){
                    #to pool across times that cell is present - this is ~1 minute per cell (e.g. many hours per embryo)
                    #CellTimeScores=sapply(MutantPositions[MutantPositions[,1]==X,2],function(Y){

                    
                    myTimes=MutantPositions[MutantPositions[,1]==X,2]
                    #This uses four time points for each cell
                    testTimes=c(myTimes[0.1*length(myTimes)],myTimes[0.35*length(myTimes)],myTimes[0.7*length(myTimes)],myTimes[0.9*length(myTimes)])

                    #require two time points (but this is in SCD)
                    if(!(length(myTimes) > 1)){return(NA)}

                    median(sapply(testTimes,function(Y){
                                        #Neighbor score at that time
                                        Neighbors=rownames(MutantPositions)[MutantPositions[,2]==Y]
                                        theseNeighborDistances=sqrt(
                                            (mutantX[paste(X,Y,sep=":"),i]-mutantX[Neighbors,i])^2+
                                            (mutantY[paste(X,Y,sep=":"),i]-mutantY[Neighbors,i])^2+
                                            (mutantZ[paste(X,Y,sep=":"),i]-mutantZ[Neighbors,i])^2)
                                        
                                        theseWTNeighborDistances=sqrt(
                                            (wtXMeans[paste(X,Y,sep=":")]-wtXMeans[Neighbors])^2+
                                            (wtYMeans[paste(X,Y,sep=":")]-wtYMeans[Neighbors])^2+
                                            (wtZMeans[paste(X,Y,sep=":")]-wtZMeans[Neighbors])^2)
                                        
                                        #THIS SCORE IS THE AVERAGE OF THE 10 MOST DEVIANT CELLS FOR LOG2 FOLD CHANGE IN DISTANCE
                                        mean(tail(sort(abs(log2(theseNeighborDistances/theseWTNeighborDistances))),10),na.rm=T)
                                        }),na.rm=T)
                    })
        
        oldColnames=colnames(NN_distance_summaries)
        NN_distance_summaries=cbind(NN_distance_summaries,CellDeviationScores)
        colnames(NN_distance_summaries)=c(oldColnames,i)
        }
    write.csv(NN_distance_summaries,file=paste(Name,"_NN_Scores.csv",sep=""))
     }



    MutantDevs=MutantPositions[,1:2]
    MutantDevs=cbind(MutantDevs,WTDmeans[rownames(MutantPositions)],WTDsds[rownames(MutantPositions)],WTDcounts[rownames(MutantPositions)])
    MutantCellDevMaxs=data.frame(Cells)
    MutantCellDevMeans=data.frame(Cells)
    colnames(MutantDevs)=c("Cell","Time","WT D mean", "WT D SD", "WT D count")
    
    RotatedX=NULL
    RotatedY=NULL
    RotatedZ=NULL
    message("optimizing rotation")
    for(i in mutantEmbryos){

        print(i)
        pdf(paste(Name,i,"PositionPlots.pdf",sep="_")) 
        
        Thetas=NULL
        Distances=NULL
        UnrotatedDistances=NULL
        RotatedMutant=NULL
        for(j in sort(Times)){

            #Find optimum rotation around x axis
            Pos=data.frame(mutantX[MutantPositions[2]==j,i],mutantY[MutantPositions[2]==j,i],mutantZ[MutantPositions[2]==j,i],row.names=MutantPositions[MutantPositions[2]==j,1])
            colnames(Pos)=c("X","Y","Z")

             if(length(Pos[,1])==0){
                 BestTheta=0
                 MeanDist=NA
                 UnrotatedMeanDist=NA
                 Thetas=rbind(Thetas,c(j,BestTheta))
                 rownames(Thetas)=Thetas[,1]
                 Distances=rbind(Distances,c(j,MeanDist))
                 rownames(Distances)=Distances[,1]
                 UnrotatedDistances=rbind(UnrotatedDistances,c(j,UnrotatedMeanDist))
                 rownames(UnrotatedDistances)=UnrotatedDistances[,1]
                 next
                 }
            WTMeanPos=data.frame(wtXMeans[WTPositions[2]==j],wtYMeans[WTPositions[2]==j],wtZMeans[WTPositions[2]==j],row.names=WTPositions[WTPositions[2]==j,1])
            colnames(WTMeanPos)=c("X","Y","Z")
            WTMeanPos=WTMeanPos[rownames(Pos),]
            WTMutantDevs=WTMeanPos-Pos
            WTMutantDists=sqrt(rowSums(WTMutantDevs^2))
            MeanDist=mean(WTMutantDists,na.rm=T)
            UnrotatedMeanDist=MeanDist
            BestTheta=0
            if(j>2){
                BestTheta=Thetas[as.character(j-1),2]
                }
            if(is.na(BestTheta)){
                BestTheta=0
                }
            prevTheta=BestTheta
            testRange=pi*(-60:60)/60
            DistOpt=c(0,MeanDist)
            bestPos=Pos
            for(k in testRange){
                
                #y' = y*cos q - z*sin q
                #z' = y*sin q + z*cos q
                #thesePos=data.frame(Pos[,"X"],Pos[,"Y"]*cos(k)-Pos[,"Z"]*sin(k),Pos[,"Z"]*cos(k)+Pos[,"Y"]*sin(k))
                thesePos=rbind(rotate3d(as.matrix(Pos[,1:3]),k,1,0,0))

                colnames(thesePos)=c("X","Y","Z")
                rownames(thesePos)=rownames(Pos)
                
                theseWTMutantDevs=WTMeanPos-thesePos
                theseWTMutantDists=sqrt(rowSums(theseWTMutantDevs^2))
                thisMeanDist=mean(sqrt(rowSums(theseWTMutantDevs^2)),na.rm=T)
                
                DistOpt=rbind(DistOpt,c(k,thisMeanDist))
                if(!is.na(thisMeanDist)&&!is.na(MeanDist)){
                    if(thisMeanDist<MeanDist){
                        BestTheta=k
                        bestPos=thesePos
                        WTMutantDevs=theseWTMutantDevs
                        WTMutantDists=theseWTMutantDists
                        MeanDist=thisMeanDist
                        }
                            
                    }
                }
            Thetas=rbind(Thetas,c(j,BestTheta))
            rownames(Thetas)=Thetas[,1]
            Distances=rbind(Distances,c(j,MeanDist))
            rownames(Distances)=Distances[,1]
            UnrotatedDistances=rbind(UnrotatedDistances,c(j,UnrotatedMeanDist))
            rownames(UnrotatedDistances)=UnrotatedDistances[,1]
            
            #Update rotated positions
            rownames(bestPos)=paste(rownames(bestPos),j,sep=":")
            RotatedMutant=rbind(RotatedMutant,bestPos)

#            if(sum(!is.na(DistOpt[,2]))>0){
#                plot(DistOpt,main=j)
#                abline(v=BestTheta,col="red")
#                abline(v=prevTheta,col="blue")
#                }
            #TODO/FIXME: Should also optimize x/y/z scaling after rotation!

            }
        plot(Thetas,main=i,xlab="Time", ylab="Least square distance Theta (rotation)",cex=0.5)
        plot(UnrotatedDistances,main=i, xlab="Time", ylab="Mean deviations (Std black, rotated red)", cex=0.5)
        points(Distances,col=2, cex=0.5)

        #version in same order as mutantX
        RotatedMutant=RotatedMutant[rownames(mutantX),]
        
        rownames(RotatedMutant)=rownames(mutantX)
 
        #Plot # of WT and mutant cells present per time point for comparison with previous plots
        WTCellCounts=sapply(Times,function(X){
                sum(!is.na(wtXMeans[WTPositions[,2]==X]))
                })
        plot(Times,WTCellCounts,xlab="Time", ylab="Cell counts - WT model (black), This embryo (red)")
        MutantCellCounts=sapply(Times,function(X){
                sum(!is.na(RotatedMutant[MutantPositions[,2]==X,1]))
                })

        points(Times,MutantCellCounts,col=2)


        #Eliminate time points where <25% the WT cell count is present and replot
        Keep=MutantCellCounts/WTCellCounts>PercentCellsRequired
        
        Keep2=Keep[MutantPositions[rownames(RotatedMutant),2]]
        names(Keep2)=rownames(RotatedMutant)
        RotatedMutant[Keep2==FALSE&!is.na(RotatedMutant[,2]),]<-c(NA,NA,NA)


        WTCellCounts=sapply(Times,function(X){
                sum(!is.na(wtXMeans[WTPositions[,2]==X]))
                })
        plot(Times,WTCellCounts,xlab="Time", ylab="Cell counts - WT model (black), this embryo after filtering 25% cell count (red)")
        MutantCellCounts=sapply(Times,function(X){
                sum(!is.na(RotatedMutant[MutantPositions[,2]==X,1]))
                })

        points(Times,MutantCellCounts,col=2)
 
        plot(wtXMeans[rownames(mutantX)],mutantX[,i], cex=0.1, xlab="WT mean X position", ylab="This embryo X position")


        plot(wtYMeans[rownames(mutantX)],mutantY[,i], cex=0.1, xlab="WT mean Y position", ylab="this embryo Y position (red=rotated)")
        points(wtYMeans[rownames(mutantX)],RotatedMutant[rownames(mutantX),2],col=2, cex=0.1)


        plot(wtZMeans[rownames(mutantX)],mutantZ[,i], cex=0.1, xlab="WT mean Z position", ylab="this embryo Z position (red=rotated)")
        points(wtZMeans[rownames(mutantX)],RotatedMutant[rownames(mutantX),3],col=2, cex=0.1)


        # Mutant_CellCell_Distances = sapply(rownames(Pos),function(X){
        #         sqrt((Pos[,"X"]-Pos[X,"X"])^2+(Pos[,"Y"]-Pos[X,"Y"])^2+(Pos[,"Z"]-Pos[X,"Z"])^2)
        #             })  
        # rownames(Mutant_CellCell_Distances)=rownames(Pos)
        # colnames(Mutant_CellCell_Distances)=rownames(Pos)
        
        # WT_CellCell_Distances = sapply(rownames(Pos),function(X){
        #         sqrt((WTMeanPos[,"X"]-WTMeanPos[X,"X"])^2+(WTMeanPos[,"Y"]-WTMeanPos[X,"Y"])^2+(WTMeanPos[,"Z"]-WTMeanPos[X,"Z"])^2)
        #         })  
        # rownames(WT_CellCell_Distances)=rownames(Pos)
        # colnames(WT_CellCell_Distances)=rownames(Pos)

       #Commented out because full data plots are too cumbersome in general - instead use cell means
        #plot(wtXMeans[CellTimes],mutantX[CellTimes,i])
        #plot(wtYMeans[CellTimes],mutantY[CellTimes,i])
        #plot(wtZMeans[CellTimes],mutantZ[CellTimes,i])
        
        xCellMeans=tapply(mutantX[,i],MutantPositions[,1],mean,na.rm=T)
        yCellMeans=tapply(mutantY[,i],MutantPositions[,1],mean,na.rm=T)
        zCellMeans=tapply(mutantZ[,i],MutantPositions[,1],mean,na.rm=T)

        yRCellMeans=tapply(RotatedMutant[rownames(MutantPositions),2],MutantPositions[,1],mean,na.rm=T)
        zRCellMeans=tapply(RotatedMutant[rownames(MutantPositions),3],MutantPositions[,1],mean,na.rm=T)

        RotatedX=cbind(RotatedX,RotatedMutant[rownames(MutantPositions),1])
        RotatedY=cbind(RotatedY,RotatedMutant[rownames(MutantPositions),2])
        RotatedZ=cbind(RotatedZ,RotatedMutant[rownames(MutantPositions),3])

        plot(wtXCellMeans[Cells],xCellMeans[Cells],main=paste(i,"Rsq",round(summary(lm(wtXCellMeans[Cells]~xCellMeans[Cells]))$adj.r.squared,5),sep=":"))
        plot(wtYCellMeans[Cells],yCellMeans[Cells],main=paste(i,"Rsq +/- rotation",round(summary(lm(wtYCellMeans[Cells]~yCellMeans[Cells]))$adj.r.squared,5),
                                                              round(summary(lm(wtYCellMeans[Cells]~yRCellMeans[Cells]))$adj.r.squared,5),sep=":"))
        points(wtYCellMeans[Cells],yRCellMeans[Cells],col=2)
        plot(wtZCellMeans[Cells],zCellMeans[Cells],main=paste(i,"Rsq +/- rotation",round(summary(lm(wtZCellMeans[Cells]~zCellMeans[Cells]))$adj.r.squared,5),
                                                              round(summary(lm(wtZCellMeans[Cells]~zRCellMeans[Cells]))$adj.r.squared,5),sep=":"))
        points(wtZCellMeans[Cells],zRCellMeans[Cells],col=2)

        print("Calculating deviations")

        #Calculate X,Y,Z deviations from expectation (distance and Z Score)
#        xDevs=mutantX[,i]-wtXMeans[rownames(MutantPositions)]
        xDevs=RotatedMutant[rownames(MutantPositions),1]-wtXMeans[rownames(MutantPositions)]

        #This only works inside this loop because tapply only works on vectors.  outside the loop 
        #something like by(xZ,MutantPositions[,1],colMeans,na.rm=T) might work but output is awkward
        xDevMeans=tapply(xDevs,MutantPositions[,1],mean,na.rm=T)
        xDevMax=tapply(xDevs,MutantPositions[,1],max,na.rm=T)
        xDevMaxTime=tapply(xDevs,MutantPositions[,1],function(X){
                match(max(X),X)
                })

        xZ=abs(xDevs)/wtX_SD[rownames(MutantPositions)]
        xZMeans=tapply(xZ,MutantPositions[,1],mean,na.rm=T)
        xZMax=tapply(xZ,MutantPositions[,1],max,na.rm=T)
        xZMaxTime=tapply(xZ,MutantPositions[,1],function(X){
                match(max(X),X)
                })
        #SAME FOR Y
#        yDevs=mutantY[,i]-wtYMeans[rownames(MutantPositions)]
        yDevs=RotatedMutant[rownames(MutantPositions),2]-wtYMeans[rownames(MutantPositions)]
        yDevMeans=tapply(yDevs,MutantPositions[,1],mean,na.rm=T)
        yDevMax=tapply(yDevs,MutantPositions[,1],max,na.rm=T)
        yDevMaxTime=tapply(yDevs,MutantPositions[,1],function(X){
                match(max(X),X)
                })
        yZ=abs(yDevs)/wtY_SD[rownames(MutantPositions)]
        yZMeans=tapply(yZ,MutantPositions[,1],mean,na.rm=T)
        yZMax=tapply(yZ,MutantPositions[,1],max,na.rm=T)
        yZMaxTime=tapply(yZ,MutantPositions[,1],function(X){
                match(max(X),X)
                })
        
        #SAME FOR Z
#        zDevs=mutantZ[,i]-wtZMeans[rownames(MutantPositions)]
        zDevs=RotatedMutant[rownames(MutantPositions),3]-wtZMeans[rownames(MutantPositions)]
        zDevMeans=tapply(zDevs,MutantPositions[,1],mean,na.rm=T)
        zDevMax=tapply(zDevs,MutantPositions[,1],max,na.rm=T)
        zDevMaxTime=tapply(zDevs,MutantPositions[,1],function(X){
                match(max(X),X)
                })
        zZ=abs(zDevs)/wtZ_SD[rownames(MutantPositions)]
        zZMeans=tapply(zZ,MutantPositions[,1],mean,na.rm=T)
        zZMax=tapply(zZ,MutantPositions[,1],max,na.rm=T)
        zZMaxTime=tapply(zZ,MutantPositions[,1],function(X){
                match(max(X),X)
                })

        #Pure Distances    
        #WTDmeans=apply(WTDdevs,1,mean,na.rm=T)
        #WTDsds=apply(WTDdevs,1,sd,na.rm=T)

        dDevs=sqrt(xDevs^2+yDevs^2+zDevs^2)
        dDevMeans=tapply(dDevs,MutantPositions[,1],mean,na.rm=T)
        dDevMax=tapply(dDevs,MutantPositions[,1],max,na.rm=T)
        dDevMaxTime=tapply(dDevs,MutantPositions[,1],function(X){
                match(max(X),X)
                })

        dZ=(dDevs-WTDmeans[rownames(MutantPositions)])/WTDsds[rownames(MutantPositions)]
        dZMeans=tapply(dZ,MutantPositions[,1],mean,na.rm=T)
        dZMax=tapply(dZ,MutantPositions[,1],max,na.rm=T)
        dZMaxTime=tapply(dZ,MutantPositions[,1],function(X){
                match(max(X),X)
                })

        #Some plots - need to label with which series
        #deviations in 2 dimensions
        plot(xDevMeans,yDevMeans,main=i)
        plot(xDevMeans,zDevMeans,main=i)
        plot(yDevMeans,zDevMeans,main=i)

        #cell mean vs cell max deviations  - FIXME the "max" function needs to do absolute max for x/y/z
        plot(xDevMeans,xDevMax,main=i)
        plot(yDevMeans,yDevMax,main=i)
        plot(zDevMeans,zDevMax,main=i)
        plot(dDevMeans,dDevMax,main=i)

        #cell max deviation vs cell max z score
        plot(xDevMax,xZMax,main=i)
        plot(yDevMax,yZMax,main=i)
        plot(zDevMax,zZMax,main=i)
        plot(dDevMax,dZMax,main=i,ylim=c(-2,20))
        plot(dDevMeans,dZMeans,main=i)


        plot(dDevMeans,abs(xDevMeans),main=i)
        plot(dDevMeans,abs(yDevMeans),main=i)
        plot(dDevMeans,abs(zDevMeans),main=i)

        oldCols=colnames(MutantDevs)
        MutantDevs=cbind(MutantDevs,dDevs)
        colnames(MutantDevs)=c(oldCols,i)

        oldCDcols=colnames(MutantCellDevMaxs)
        test=(dDevMax[Cells])
        names(test)=Cells
        MutantCellDevMaxs=cbind(MutantCellDevMaxs,test)
        test=(dDevMeans[Cells])
        names(test)=Cells
        MutantCellDevMeans=cbind(MutantCellDevMeans,test)
        colnames(MutantCellDevMaxs)=c(oldCDcols,i)
        colnames(MutantCellDevMeans)=c(oldCDcols,i)

        #END LOOP THROUGH MUTANT SERIES
        dev.off()
        }

    write.csv(MutantCellDevMaxs,file=paste(Name,"CellMaxPositionDevs.csv",sep="_"))
    write.csv(MutantCellDevMeans,file=paste(Name,"CellMeanPositionDevs.csv",sep="_"))
    write.csv(MutantDevs,file=paste(Name,"PositionDevs.csv",sep="_"))
    
    #FIXME output MutantDevs as CD file so can be plotted in lineage form

    #list d dev per series per cells 

    print("Reading NN Scores")
    print(paste(Name,"_NN_Scores.csv",sep=""))

    MC_NN_scores=read.csv(paste(Name,"_NN_Scores.csv",sep=""),row.names=1)
    MC_NN_scores[,1]<-NULL
    
    print("Reading WT NN Scores")
    print(paste(WTDir,"Richard_et_al_plus_comma_WT_NN_Scores.csv",sep="/"))

    WT_NN_scores=read.csv(paste(WTDir,"Richard_et_al_plus_comma_WT_NN_Scores.csv",sep="/"),row.names=1)
    WT_NN_scores[,1]<-NULL

    print("Reading mean dev Scores")
    print(paste(Name,"_CellMeanPositionDevs.csv",sep=""))

    MC_Devs=read.csv(paste(Name,"_CellMeanPositionDevs.csv",sep=""),row.names=1)
    MC_Devs[,1]<-NULL

    print("Reading mean dev Scores")
    print(paste(WTDir,"Richard_et_al_plus_comma_WT_CellMeanPositionDevs.csv",sep="/"))

    WT_Devs=read.csv(paste(WTDir,"Richard_et_al_plus_comma_WT_CellMeanPositionDevs.csv",sep="/"),row.names=1)
    WT_Devs[,1]<-NULL

    colnames(MC_NN_scores)=mutantEmbryos    
    colnames(MC_Devs)=mutantEmbryos
    colnames(WT_NN_scores)=wtEmbryos
    colnames(WT_Devs)=wtEmbryos

    print("Outputting defect summary")

    pdf(paste(Name,"positionDefects.pdf", sep="_"))
    Cells=rownames(MC_NN_scores)
    NN_wilcox_p=sapply(Cells,function(X){
            try(return(wilcox.test(as.numeric(MC_NN_scores[X,]),as.numeric(WT_NN_scores[X,]))$p.value),silent=T)
            return(NA)
            })
    Dev_wilcox_p=sapply(Cells,function(X){
            try(return(wilcox.test(as.numeric(MC_Devs[X,]),as.numeric(WT_Devs[X,]))$p.value),silent=T)
            return(NA)
            })
    
    

    #Create a function to generate a continuous color palette
    rbPal <- colorRampPalette(c("black","purple","blue","green","red"))
    
    #This adds a column of color values
    # based on the exp values
    Colors=rbPal(20)[as.numeric(cut(Expression[Cells],breaks = 20))]
    
    
    #Colors=rainbow(length(Cells))
    names(Colors)=Cells
    
    plot(jitter(log10(Dev_wilcox_p)),log10(NN_wilcox_p),col=Colors[Cells],xlab="log(p value: deviation vs WT) - Wilcoxon", ylab="log(p value: deviation vs WT) - Wilcoxon", cex=0.2)
    plot(c(min(log10(Dev_wilcox_p),na.rm=T),max(log10(Dev_wilcox_p),na.rm=T)),
         c(min(log10(NN_wilcox_p),na.rm=T),max(log10(NN_wilcox_p),na.rm=T)))
    
    text(log10(Dev_wilcox_p),log10(NN_wilcox_p),labels=Cells,col=Colors[Cells])
    abline(h=-2,v=-2,col=2)
    
    
    combinedP=NN_wilcox_p*Dev_wilcox_p
    plot(Expression[Cells],log10(combinedP[Cells]), xlab="Expression (peak)", ylab="log(combined deviation P value", cex=0.1)
    
    maxP=sapply(Cells,function(X){max(NN_wilcox_p[X],Dev_wilcox_p[X])})
    plot(Expression[Cells],log10(maxP[Cells]), xlab="Expression (peak)", ylab="log(maximum P value) (NN or deviation)", cex=0.1)
    
    minP=sapply(Cells,function(X){min(NN_wilcox_p[X],Dev_wilcox_p[X])})
    plot(Expression[Cells],log10(minP[Cells]), xlab="Expression (peak)", ylab="log(minimum P value) (NN or deviation)", cex=0.1)
    
    
    
    Expressing=Expression[Cells]>expCutoff
    
    par(mfrow=c(1,3))
    boxplot(log10(minP)~Expressing,main="min",ylim=c(-15,0))
    boxplot(log10(maxP)~Expressing,main="max",ylim=c(-15,0))
    boxplot(log10(combinedP)~Expressing,main="combined",ylim=c(-15,0))
    par(mfrow=c(1,1))
    
    
    #wilcox.test(combinedP[Expressing==T],combinedP[Expressing==F])
    #wilcox.test(maxP[Expressing==T],maxP[Expressing==F])
    
    
    #Now find, count and plot individual outliers
    NN_Z=(MC_NN_scores[Cells,]-apply(WT_NN_scores[Cells,],1,mean,na.rm=T))/apply(WT_NN_scores[Cells,],1,sd,na.rm=T)
    WT_NN_Z=(WT_NN_scores[Cells,]-apply(WT_NN_scores[Cells,],1,mean,na.rm=T))/apply(WT_NN_scores[Cells,],1,sd,na.rm=T)
    
    #Fixme should do by cross validation and allow alternative control/WT sets
    max_wt_nn_Z=apply(WT_NN_Z,1,max,na.rm=T)
    max_nn_Z=apply(NN_Z,1,max,na.rm=T)

    Dev_Z=(MC_Devs[Cells,]-apply(WT_Devs[Cells,],1,mean,na.rm=T))/apply(WT_Devs[Cells,],1,sd,na.rm=T)
    WT_Dev_Z=(WT_Devs[Cells,]-apply(WT_Devs[Cells,],1,mean,na.rm=T))/apply(WT_Devs[Cells,],1,sd,na.rm=T)
    
    max_wt_dev_Z=apply(WT_Dev_Z,1,max,na.rm=T)
    max_dev_Z=apply(Dev_Z,1,max,na.rm=T)


    plotnn<-max(max_nn_Z)*1.2
    plotdev<-max(max_dev_Z)*1.2
    plot(max_nn_Z[Cells[!Expressing]],max_dev_Z[Cells[!Expressing]],xlim=c(-5,plotnn),ylim=c(-5,plotdev), xlab="max NN Z score", ylab="max dev Z score; red = expressing cells")
    points(max_nn_Z[Cells[Expressing]],max_dev_Z[Cells[Expressing]],col=2)

    plot(Expression[Cells],log2(rowMeans(NN_Z,na.rm=T)), xlab="Expression (peak)", ylab="mean NN Z score")

    Max_MaxZ=max_nn_Z+max_dev_Z
    hist(Max_MaxZ[Cells[!Expressing]],col="#0000ff99",freq=F,breaks=20, main="Red=expressing Blue=not expressing", xlab="Max max Z score")
    hist(Max_MaxZ[Cells[Expressing]],add=T,col="#ff000099",freq=F,breaks=20)

    hist(Max_MaxZ[Cells[!Expressing]],col="#0000ff99",freq=T,breaks=20,main="Red=expressing Blue=not expressing", xlab="Max Z score")
    hist(Max_MaxZ[Cells[Expressing]],add=T,col="#ff000099",freq=T,breaks=20)

    sumZ=NN_Z+Dev_Z
    mean_sumZ=rowMeans(sumZ,na.rm=T)
    max_sumZ=apply(sumZ,1,max,na.rm=T)


    hist(mean_sumZ[Cells[!Expressing]],col="#0000ff99",freq=T,breaks=20, main="Red=expressing Blue=not expressing", xlab="Max summed Z score")
    hist(mean_sumZ[Cells[Expressing]],add=T,col="#ff000099",freq=T,breaks=20)

    hist(max_sumZ[Cells[!Expressing]],col="#0000ff99",freq=T,breaks=20, main="Red=expressing Blue=not expressing", xlab="Max summed Z score")
    hist(max_sumZ[Cells[Expressing]],add=T,col="#ff000099",freq=T,breaks=20)


    dev.off()
    #End position deviation plots - start trajectory plots

    
    rownames(RotatedX)=rownames(MutantPositions)
    rownames(RotatedY)=rownames(MutantPositions)
    rownames(RotatedZ)=rownames(MutantPositions)
    colnames(RotatedX)=mutantEmbryos
    colnames(RotatedY)=mutantEmbryos    
    colnames(RotatedZ)=mutantEmbryos

    write.csv(RotatedX,file=paste(Name,"_rotatedX.csv",sep=""))
    write.csv(RotatedY,file=paste(Name,"_rotatedY.csv",sep=""))
    write.csv(RotatedZ,file=paste(Name,"_rotatedZ.csv",sep=""))

    RotatedX2=data.frame(RotatedX)[rownames(WTPositions),]
    RotatedY2=data.frame(RotatedY)[rownames(WTPositions),]
    RotatedZ2=data.frame(RotatedZ)[rownames(WTPositions),]
    colnames(RotatedX2)=mutantEmbryos
    colnames(RotatedY2)=mutantEmbryos    
    colnames(RotatedZ2)=mutantEmbryos
    

    plotDeviation <- function(Cell,Movie,Label){
        #fixme - top set uses unrotated, bottom uses rotated
        # plot(wtXMeans[WTPositions[,1]==Cell],wtYMeans[WTPositions[,1]==Cell],xlim=c(-30,30),ylim=c(-20,20),main=Label,cex=0.5)
        # points(mutantX[MutantPositions[,1]==Cell,Movie],mutantY[MutantPositions[,1]==Cell,Movie],col="red",cex=0.5)
        # points(mutantX[MutantPositions[,1]==GetParent(Cell),Movie],mutantY[MutantPositions[,1]==GetParent(Cell),Movie],col="pink",cex=0.5)
        # points(wtXMeans[WTPositions[,1]==GetParent(Cell)],wtYMeans[WTPositions[,1]==GetParent(Cell)],col="grey",cex=0.5)
#        message(paste(Cell,Movie))
        plot(wtXMeans[WTPositions[,1]==Cell&!is.na(RotatedX2[,Movie])],
             wtYMeans[WTPositions[,1]==Cell&!is.na(RotatedX2[,Movie])],
             xlim=c(-30,30),ylim=c(-30,30),main=Label,cex=0.5,xlab="X",ylab="Y")
        points(wtXMeans[WTPositions[,1]==GetParent(Cell)&!is.na(RotatedX2[,Movie])],wtYMeans[WTPositions[,1]==GetParent(Cell)&!is.na(RotatedX2[,Movie])],col="grey",cex=0.5)

        points(RotatedX2[WTPositions[,1]==Cell,Movie],RotatedY2[WTPositions[,1]==Cell,Movie],col="red",cex=0.5)
        points(RotatedX2[WTPositions[,1]==GetParent(Cell),Movie],RotatedY2[WTPositions[,1]==GetParent(Cell),Movie],col="pink",cex=0.5)
        draw.ellipse(0,0,a=25,b=15)


        plot(wtYMeans[WTPositions[,1]==Cell&!is.na(RotatedX2[,Movie])],wtZMeans[WTPositions[,1]==Cell&!is.na(RotatedX2[,Movie])],
             xlim=c(-20,20),ylim=c(-20,20),main=Label,cex=0.5,xlab="Y",ylab="Z")
        points(wtYMeans[WTPositions[,1]==GetParent(Cell)&!is.na(RotatedX2[,Movie])],wtZMeans[WTPositions[,1]==GetParent(Cell)&!is.na(RotatedX2[,Movie])],col="grey",cex=0.5)

        points(RotatedY2[WTPositions[,1]==Cell,Movie],RotatedZ2[WTPositions[,1]==Cell,Movie],col="red",cex=0.5)
        points(RotatedY2[WTPositions[,1]==GetParent(Cell),Movie],RotatedZ2[WTPositions[,1]==GetParent(Cell),Movie],col="pink",cex=0.5)
        draw.ellipse(0,0,a=15,b=15)
        }
    

    pdf(paste(Name,"deviantCellTrajectories.pdf",sep=""))
    par(mfrow=c(2,2))
    

    for(j in Cells){
        if(max_sumZ[j]>10){
            par(mfrow=c(2,2))
            for(i in mutantEmbryos){
                if(!is.na(sumZ[j,i])&&sumZ[j,i]>10){
                    Round=round(c(sumZ[j,i],MC_NN_scores[j,i],MC_Devs[j,i]),1)
                    Label=paste(j," ", CellNames[j],"
",
                                i,"
Z:",
                                Round[1]," NN:",Round[2]," Dev:",Round[3]," Exp:",Expression[j],sep="")
                    plotDeviation(j,i,Label)
                    }
                }
            }
        par(mfrow=c(1,1))
        }
    dev.off()

    pdf(paste(Name,"CellTrajectories.pdf",sep=""))
    
    for(j in Cells){
        par(mfrow=c(2,2))
        for(i in mutantEmbryos){
            if(!is.na(sumZ[j,i])){
                Round=round(c(sumZ[j,i],MC_NN_scores[j,i],MC_Devs[j,i]),1)
                Label=paste(j," ", CellNames[j],"
",
                                i,"
Z:",
                                Round[1]," NN:",Round[2]," Dev:",Round[3]," Exp:",Expression[j],sep="")
                plotDeviation(j,i,Label)
                }
            }
            
        par(mfrow=c(1,1))
        }
    dev.off()

       
    #Also find and plot worst outlier for most significant cells
    
    #AND calculate not just position defects but also trajectory/migration defects
    setwd(startDir)

    return(MutantDevs)
    #GOALS:
        #1) Identify deviant cells at level of absolute A-P, L-R, D-V position or Total distance from predicted position (NEED TO CHECK AND FIX Z DEV ISSUE - POSSIBLE BUG?)
         #A - Update rotation correction for each time point, calculate absolute difference on optimally rotated embryos
        #2) Identify deviant cells based on position relative to neighbors
        #3) Output concise lists of deviant cells.  Also trajectory plots for deviant cells.
        #4) Think harder about use of Max and Mean as metrics to collapse cells.  
        #5) Consider applying time point by time point. E.g. predict position of a cell at a given time point based on other cells' positions in that time point?
        #6) Try aligning on polar coordinates *** Note that the current code to optimize rotations with WT mean should solve this


    #12/11/2013 action list
      #1 - X/Y/Z scaling?
      #1.5 TROUBLESHOOT WHY ARE MEDIAN DISTANCES LARGER FOR MUTANTS? - MAYBE need to start with earliest defects.  If this is a result of defective cells then altering the positions of other cells non-autonomously, would predict that the earlier time points where ceh-36 has not yet been expressed would be normal? But ceh-36 may be acting early (by 50 cell stage)
      #2 - Combined identification of defective cells. And MOST defective cells
      #3 - Visualization of defective cells' trajectories. 



}



