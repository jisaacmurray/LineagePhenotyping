source("functions.R")

#This is to get a WT Cell list
WTDivTimes=read.table("Richard_et_al_plus_comma_WT/Richard_et_al_plus_comma_WTDivTimeNorm.tsv", header=T,row.names=1,stringsAsFactors=F)
    #Remove sulston which is included in these tables for historical reasons
WTDivTimes[,1] <- NULL
Cells=rownames(WTDivTimes)

tmpCellNames=read.csv("CellNames.csv")
CellNames=as.vector(tmpCellNames[,2])
names(CellNames)=tmpCellNames[,1]
CellTimes=as.vector(tmpCellNames[,3])
names(CellTimes)=tmpCellNames[,1]

#old version based on one series - version below uses ACD average from multiple series
#ceh36_exp=read.csv("data/CA20120714_JIM136_L3.csv",row.names=2)["blot"]
#ceh36_peak=sapply(Cells,function(X){
#        GetPeak(ceh36_exp,X)
#        })  



nhr67_peak=ReadPeakExpression("data/CA20130527_nhr-67_JIM202_L2.csv")
ceh36_peak=ReadPeakExpression("data/CA20120714_JIM136_L3.csv")
mls2_peak=ReadPeakExpression("data/CA20110209_UP2051_mls-2_L2.csv")
ceh13_peak=ReadPeakExpression("data/CA20120305_ceh-13_L3.csv")
unc30_peak=ReadPeakExpression("data/CA20121115_unc-30_18_L2.csv")
nob1_peak=ReadPeakExpression("data/CA20140825_nob-1_JIM284_L3.csv")
elt1_peak=ReadPeakExpression("data/CA20120130_elt-1_17_L2.csv")

#See old version at EOF



#for(i in c("unc-30", "ceh-36_unc-30_rescued", "ceh-36_rescued", "ceh-36", "ceh-36_unc-30", "nob-1", "mls-2_mutants", "mls2ceh36.txt","Richard_et_al_plus_comma_WT")){
#for(i in c("SM2373")){

##PlotDefectSummaries("nob-1_mutant", nob1_peak,sig=5,microns=5)

##PlotDefectSummaries("ceh-13_mutant", ceh13_peak,sig=5,microns=5)


PlotDefectSummaries <- function(i, ExpPeak,sig=3,dev=5, expCutoff=200, NN_sig=3.5,max_sig=3.5,mean_sig=3.5,microns = 5, NN_cutoff=0.8,minDivTime=0){
    ExpPeak[is.na(ExpPeak)] <- 0
    CellsLineageOrder=read.csv("Cells_350min_lineageOrder.csv")

    WTMax=read.csv("Richard_et_al_plus_comma_WT/Richard_et_al_plus_comma_WT_CellMaxPositionDevs.csv",row.names=1)

    #clear "-Inf" values and remove duplicate name column
    WTMax[,1]<-NULL
    WTMax[WTMax<(-100)] <- NA

    print(i)
    print(2)

    WTMean=read.csv("Richard_et_al_plus_comma_WT/Richard_et_al_plus_comma_WT_CellMeanPositionDevs.csv",row.names=1)
    WTMean[,1]<-NULL
    WTMean[WTMean<(-100)] <- NA

    WTNN=read.csv("Richard_et_al_plus_comma_WT/Richard_et_al_plus_comma_WT_NN_Scores.csv",row.names=1)
    WTNN[,1]<-NULL


    WTMaxMean=rowMeans(WTMax,na.rm=T)
    WTMaxSD=apply(WTMax,1,sd,na.rm=T)
    WTMeanMean=rowMeans(WTMean,na.rm=T)
    WTMeanSD=apply(WTMean,1,sd,na.rm=T)
    WTNNMean=rowMeans(WTNN,na.rm=T)
    WTNNMax = apply(WTNN,1,max,na.rm=T)
    WTNNSD=apply(WTNN,1,sd,na.rm=T)

 
    ##old values 
    ##sig=3
    ##dev=5

    ##sig=5
    ##dev=6
    
    startDir=getwd()
    setwd(i)

    print(i)
    pdf(paste(i,"summary.pdf",sep="_"))
    #1)Load CC devs
    theseCCDevs=read.csv(paste(i,"_ccDevs.csv",sep=""),row.names=1)

#    theseCCDevs=theseCCDevs[theseCCDevs$Div.time > minDivTime,]
   
    Early=theseCCDevs[,"CC.Z.score"]<(-sig) & theseCCDevs[,"CC.Deviation"]<(-dev)
    names(Early)=rownames(theseCCDevs)
    Early[!is.na(theseCCDevs["WT.CC"]) & !is.na(theseCCDevs["Terminal.Cell.Length"]) | !is.na(theseCCDevs["CC"])] <-0
    Early[theseCCDevs[,"CC.Z.score"]<(-sig) & theseCCDevs[,"CC.Deviation"]<(-dev)]<-1


    Late=theseCCDevs[,"CC.Z.score"]>sig & theseCCDevs[,"CC.Deviation"]>dev
    names(Late)=rownames(theseCCDevs)
    Late[!is.na(theseCCDevs["WT.CC"]) & !is.na(theseCCDevs["Terminal.Cell.Length"]) | !is.na(theseCCDevs["CC"])] <-0
    Late[theseCCDevs[,"CC.Z.score"]>sig & theseCCDevs[,"CC.Deviation"]>dev ] <- 1
    
    Missed=theseCCDevs[,"Terminal.cell.delta"]>dev & theseCCDevs[,"Terminal.Cell.Z"]>sig
    names(Missed)=rownames(theseCCDevs)
    Missed[!is.na(theseCCDevs["WT.CC"]) & !is.na(theseCCDevs["Terminal.Cell.Length"]) | !is.na(theseCCDevs["CC"])] <-0
    Missed[theseCCDevs[,"Terminal.cell.delta"]>dev & theseCCDevs[,"Terminal.Cell.Z"]>sig]<-1
    
    Ectopic=is.na(theseCCDevs[,"WT.Div.time"]) & !is.na(theseCCDevs[,"Div.time"])
    names(Ectopic)=rownames(theseCCDevs)
    Ectopic[is.na(theseCCDevs["WT.CC"]) & !is.na(theseCCDevs["Terminal.Cell.Length"]) | !is.na(theseCCDevs["CC"])] <-0
    Ectopic[is.na(theseCCDevs[,"WT.Div.time"]) & !is.na(theseCCDevs[,"Div.time"])]<-1

    Ectopic[theseCCDevs[,"Div.time"] < minDivTime] <- 0
    Missed[theseCCDevs[,"Div.time"] < minDivTime] <- 0
    Early[theseCCDevs[,"Div.time"] < minDivTime] <- 0
    Late[theseCCDevs[,"Div.time"] < minDivTime] <- 0
    
    ###This reports how many cells in the dataset have each defect class (based on the selected thresholds). For example if ABplpppp divided late in 2 embryos, it would contribute 2 to the "late" count
    DividingCells = unique(theseCCDevs[!is.na(theseCCDevs$CC) | !is.na(theseCCDevs$WT.CC),1])
    print(paste("Number of cells dividing in WT or mutant: ",length(DividingCells),
                  " Expressing: ", sum(ExpPeak[DividingCells]>=expCutoff),
                  " Nonexpressing: ", sum(ExpPeak[DividingCells]<expCutoff)))

    print("NUMBER OF CELL-EMBRYO DEFECTS")
    print(paste("Early:",sum(Early,na.rm=T)))
    print(paste("Late:",sum(Late,na.rm=T)))
    print(paste("Missed:",sum(Missed,na.rm=T)))
    print(paste("Ectopic:",sum(Ectopic,na.rm=T)))
    

    #CombinedDefects=Early==1 | Late==1 | Missed==1 | Ectopic==1
    CombinedDefects=pmax(Early,Late,Missed,Ectopic)
    
    names(CombinedDefects)=rownames(theseCCDevs)
 
    Cells=as.vector(unique(theseCCDevs[,1]))
    Embs=as.vector(unique(theseCCDevs[,2]))
    CCDefectArray=NULL
    MissedArray=NULL
    LateArray=NULL
    EctopicArray=NULL
    EarlyArray=NULL
    for(j in Embs){
        CCDefectArray=cbind(CCDefectArray,CombinedDefects[paste(Cells,j,sep="_")])
        MissedArray=cbind(MissedArray,Missed[paste(Cells,j,sep="_")])
        LateArray=cbind(LateArray,Late[paste(Cells,j,sep="_")])
        EctopicArray=cbind(EctopicArray,Ectopic[paste(Cells,j,sep="_")])
        EarlyArray=cbind(EarlyArray,Early[paste(Cells,j,sep="_")])
        
        }
   
    rownames(CCDefectArray)=Cells
    rownames(MissedArray)=Cells
    rownames(EarlyArray)=Cells
    rownames(LateArray)=Cells
    rownames(EctopicArray)=Cells
    
    CCDefectCount=rowSums(CCDefectArray,na.rm=T)



    print("NUMBER OF CELLS DEFECTIVE BY EXPRESSION")
    print("EXPRESSING")
    print(paste("Early: ", sum(EarlyArray[ExpPeak[Cells]>=expCutoff,],na.rm=T)))
    print(paste("Late: ", sum(LateArray[ExpPeak[Cells]>=expCutoff,],na.rm=T)))
    print(paste("Missed: ", sum(MissedArray[ExpPeak[Cells]>=expCutoff,],na.rm=T)))
    print(paste("Ectopic: ", sum(EctopicArray[ExpPeak[Cells]>=expCutoff,],na.rm=T)))
    print(paste("ANY DEFECT: ", sum(CCDefectArray[ExpPeak[Cells]>=expCutoff,],na.rm=T)))
    print("NOT EXPRESSING")
    print(paste("Early: ", sum(EarlyArray[ExpPeak[Cells]<expCutoff,],na.rm=T)))
    print(paste("Late: ", sum(LateArray[ExpPeak[Cells]<expCutoff,],na.rm=T)))
    print(paste("Missed: ", sum(MissedArray[ExpPeak[Cells]<expCutoff,],na.rm=T)))
    print(paste("Ectopic: ", sum(EctopicArray[ExpPeak[Cells]<expCutoff,],na.rm=T)))
    print(paste("ANY DEFECT: ", sum(CCDefectArray[ExpPeak[Cells]<expCutoff,],na.rm=T)))
    print("NUMBER OF CELLS DEFECTIVE IN ANY EMBRYO BY EXPRESSION ")
    print("EXPRESSING")
    print(paste("Early: ", sum(rowSums(EarlyArray[ExpPeak[Cells]>=expCutoff,],na.rm=T)>0)))
    print(paste("Late: ", sum(rowSums(LateArray[ExpPeak[Cells]>=expCutoff,],na.rm=T)>0)))
    print(paste("Missed: ", sum(rowSums(MissedArray[ExpPeak[Cells]>=expCutoff,],na.rm=T)>0)))
    print(paste("Ectopic: ", sum(rowSums(EctopicArray[ExpPeak[Cells]>=expCutoff,],na.rm=T)>0)))
    print(paste("ANY DEFECT: ", sum(rowSums(CCDefectArray[ExpPeak[Cells]>=expCutoff,],na.rm=T)>0)))
    print("NOT EXPRESSING")
    print(paste("Early: ", sum(rowSums(EarlyArray[ExpPeak[Cells]<expCutoff,],na.rm=T)>0)))
    print(paste("Late: ", sum(rowSums(LateArray[ExpPeak[Cells]<expCutoff,],na.rm=T)>0)))
    print(paste("Missed: ", sum(rowSums(MissedArray[ExpPeak[Cells]<expCutoff,],na.rm=T)>0)))
    print(paste("Ectopic: ", sum(rowSums(EctopicArray[ExpPeak[Cells]<expCutoff,],na.rm=T)>0)))
    print(paste("ANY DEFECT: ", sum(rowSums(CCDefectArray[ExpPeak[Cells]<expCutoff,],na.rm=T)>0)))
    print("NUMBER OF CELLS DEFECTIVE IN AT LEAST TWO EMBRYOS BY EXPRESSION ")
    print("EXPRESSING")
    print(paste("Early: ", sum(rowSums(EarlyArray[ExpPeak[Cells]>=expCutoff,],na.rm=T)>1)))
    print(paste("Late: ", sum(rowSums(LateArray[ExpPeak[Cells]>=expCutoff,],na.rm=T)>1)))
    print(paste("Missed: ", sum(rowSums(MissedArray[ExpPeak[Cells]>=expCutoff,],na.rm=T)>1)))
    print(paste("Ectopic: ", sum(rowSums(EctopicArray[ExpPeak[Cells]>=expCutoff,],na.rm=T)>1)))
    print(paste("ANY DEFECT: ", sum(rowSums(CCDefectArray[ExpPeak[Cells]>=expCutoff,],na.rm=T)>1)))
    print("NOT EXPRESSING")
    print(paste("Early: ", sum(rowSums(EarlyArray[ExpPeak[Cells]<expCutoff,],na.rm=T)>1)))
    print(paste("Late: ", sum(rowSums(LateArray[ExpPeak[Cells]<expCutoff,],na.rm=T)>1)))
    print(paste("Missed: ", sum(rowSums(MissedArray[ExpPeak[Cells]<expCutoff,],na.rm=T)>1)))
    print(paste("Ectopic: ", sum(rowSums(EctopicArray[ExpPeak[Cells]<expCutoff,],na.rm=T)>1)))
    print(paste("ANY DEFECT: ", sum(rowSums(CCDefectArray[ExpPeak[Cells]<expCutoff,],na.rm=T)>1)))


    ##List repeatedly defective cells
    print("Repeat defect cells (expressing): ")
    print(rowSums(CCDefectArray[rowSums(CCDefectArray,na.rm=T)>1 & ExpPeak[Cells]>= expCutoff,],na.rm=T))
    print("Repeat defect cells (non expressing): ")
    print(rowSums(CCDefectArray[rowSums(CCDefectArray,na.rm=T)>1 & ExpPeak[Cells] <expCutoff,],na.rm=T))

    #code to determine for each cell if that cell or an ancestor was defective
    AncestorHadDefect <- function(Cell,DefectArray){
        Defects=DefectArray[Cell,]
        Defects[is.na(Defects)]<-0
        Parent=GetParent(Cell)
#        print(Parent)
        if(Parent=="P" || is.na(Parent)){
            return(Defects)
            }else{
            return(Defects | AncestorHadDefect(Parent,DefectArray))
            }
        }
    print(.5)
    AncestralDefects=data.frame()
    for(j in Cells){
#        print("CELL")
#        print(i)
        if(!is.na(j)){
            AncestralDefects=rbind(AncestralDefects,AncestorHadDefect(j,CCDefectArray))
            }
        }

    rownames(AncestralDefects)=Cells
    colnames(AncestralDefects)=Embs
#
    ExpDef=sum(ExpPeak[Cells]>expCutoff & CCDefectCount>0,na.rm=T)
    ExpOK=sum(ExpPeak[Cells]>expCutoff & CCDefectCount==0,na.rm=T)
    NotExpDef=sum(ExpPeak[Cells]<expCutoff & CCDefectCount>0,na.rm=T)
    NotExpOK=sum(ExpPeak[Cells]<expCutoff & CCDefectCount==0,na.rm=T)
##    DivisionsAssessed =sum(
        
    plot(ExpPeak[Cells],jitter(CCDefectCount),xlab="Peak Expression", ylab=paste(i,"CC Defect Count"),
         main=paste("Sigma Cutoff: ", sig, "Dev cutoff: ", dev, "\nExpressing:", ExpDef, ExpOK, round(ExpDef/(ExpDef+ExpOK),3), "\nNot Expressing:", NotExpDef, NotExpOK, round(NotExpDef/(NotExpDef+NotExpOK),3), "\n"), cex=0.1)   

    cat(paste("CELL CYCLE\nSigma Cutoff: ", sig,
              "Dev cutoff (min): ", dev,
              "\nExpressing:", ExpDef, ExpOK, round(ExpDef/(ExpDef+ExpOK),3),
              "\nNot Expressing:", NotExpDef, NotExpOK, round(NotExpDef/(NotExpDef+NotExpOK),3), "\n"))


    ##2)Load position devs

    
    print(i)
    paste(i,"_CellMaxPositionDevs.csv",sep="")
    theseMaxDevs=read.csv(paste(i,"_CellMaxPositionDevs.csv",sep=""),row.names=1)
    theseDivTimes=read.table(paste(i,"DivTime.tsv",sep=""),row.names=1,sep="\t")
    theseDivTimes <- subset(theseDivTimes, select=-1)
    theseDivTimes <- theseDivTimes[Cells,]
        
    theseMeanDevs=read.csv(paste(i,"_CellMeanPositionDevs.csv",sep=""),row.names=1)
    theseNNs=read.csv(paste(i,"_NN_Scores.csv",sep=""),row.names=1)
    #because these csv files have 2 columns with cell name duplicated
    theseMaxDevs[,1] <-NULL
    theseMeanDevs[,1] <-NULL
    theseNNs[,1]<-NULL
    theseMaxDevs[theseMaxDevs<(-100)]<-NA
    theseMeanDevs[theseMeanDevs<(-100)]<-NA

    
    theseMaxZ=(theseMaxDevs[Cells,]-WTMaxMean[Cells])/WTMaxSD[Cells]
    theseMeanZ=(theseMeanDevs[Cells,]-WTMeanMean[Cells])/WTMeanSD[Cells]
    theseNNZ=(theseNNs[Cells,]-WTNNMean[Cells])/WTNNSD[Cells]
    
    theseMaxMean=rowMeans(theseMaxDevs,na.rm=T)
    theseMaxSD=apply(theseMaxDevs,1,sd,na.rm=T)

    theseMeanMean=rowMeans(theseMeanDevs,na.rm=T)
    theseMeanSD=apply(theseMeanDevs,1,sd,na.rm=T)

    print(3)

    theseNNMean=rowMeans(theseNNs,na.rm=T)
    theseNNSD=apply(theseNNs,1,sd,na.rm=T)
    theseNNMax=apply(theseNNs,1,max,na.rm=T)

    cellObservedCounts=rowSums(!is.na(theseMaxDevs[Cells,]),na.rm=T)
    
    #plots vs SD
    plot(theseMaxMean,theseMaxSD,xlab=paste(i,"max position dev"),ylab="max position dev SD", main = "Black = mutant, Red = WT")
    points(WTMaxMean,WTMaxSD,col=2)

    plot(theseMeanMean,theseMeanSD,xlab=paste(i,"mean position dev"),ylab="mean position dev SD", main = "Black = mutant, Red = WT")
    points(WTMeanMean,WTMeanSD,col=2)

    plot(theseNNMean,theseNNSD,xlab=paste(i,"NN score"),ylab="NN SD", main = "Black = mutant, Red = WT")
    points(WTNNMean,WTNNSD,col=2)
    
    #mean plots vs WT
    plot(WTMaxMean[Cells],theseMaxMean[Cells],main=i, xlab="Mean deviation, WT", ylab="Mean deviation, mutant")
    plot(WTMeanMean[Cells],theseMeanMean[Cells],main=i, xlab="Mean max deviation, WT", ylab="Mean max deviation, mutant")
    plot(WTNNMean[Cells],theseNNMean[Cells],main=i, xlab="Mean NN score, WT", ylab="Mean NN score, mutant")

    plot(theseMaxMean[Cells],theseNNMean[Cells], xlab="max position dev", ylab="Mean NN score", main="Black = mutant, Red = WT")
    points(WTMaxMean[Cells],WTNNMean[Cells],col=2)
    
    #outlier analysis
    
    PosOutliers=(theseMaxZ>max_sig | theseMeanZ>mean_sig) & theseNNZ>NN_sig & theseMaxDevs[Cells,]>microns & ( is.na(theseDivTimes) | theseDivTimes > minDivTime)
    PosOutliers[PosOutliers==T]<-1
    PosOutliers[PosOutliers==F]<-0

    rownames(PosOutliers)=Cells
    PosDefectCount=rowSums(PosOutliers,na.rm=T)
    
    plot(jitter(PosDefectCount[Cells]),jitter(CCDefectCount[Cells]),
         main=paste(i, "R=", round(cor(PosDefectCount[Cells],CCDefectCount[Cells],use="pairwise.complete.obs",method="spearman"),digits=4)),
         xlab="Position defects", ylab="CC defects")
    plot(jitter(PosDefectCount[Cells]), jitter(CCDefectCount[sapply(Cells,GetParent)]),
         main=paste(i, "R=",round(cor(PosDefectCount[Cells],CCDefectCount[sapply(Cells,GetParent)],use="pairwise.complete.obs",method="spearman"),digits=4)),
         xlab="Position defects", ylab="Parent CC defects")
    

    plot(jitter(PosDefectCount[Cells]), jitter(CCDefectCount[sapply(Cells,GetParent)] + CCDefectCount[Cells]),
         main=paste(i, "R=",round(cor(PosDefectCount[Cells],CCDefectCount[sapply(Cells,GetParent)] + CCDefectCount[Cells],use="pairwise.complete.obs",method="spearman"),digits=4)),
         xlab="Position defects", ylab="This cell + parent CC defects")
    

    
    ExpPosDef=sum(ExpPeak[Cells]>expCutoff & PosDefectCount>0,na.rm=T)
    ExpPosOK=sum(ExpPeak[Cells]>expCutoff & PosDefectCount==0,na.rm=T)
    NotExpPosDef=sum(ExpPeak[Cells]<expCutoff & PosDefectCount>0,na.rm=T)
    NotExpPosOK=sum(ExpPeak[Cells]<expCutoff & PosDefectCount==0,na.rm=T)
    
    plot(ExpPeak[Cells],jitter(PosDefectCount[Cells]), xlab="Peak expression", ylab="Position defect count",
         main=paste("NN Cutoff: ", NN_sig, "max cutoff: ", max_sig,"mean cutoff: ", mean_sig, "micron cutoff:", microns,"\nPOSITION DEFECTS (defect / no defect / fractionDefect)",
                    "\nExpressing:", ExpPosDef, ExpPosOK, round(ExpPosDef/(ExpPosDef+ExpPosOK),3), 
                    "\nNot Expressing:", NotExpPosDef, NotExpPosOK, round(NotExpPosDef/(NotExpPosDef+NotExpPosOK),3), "\n"))   
    
    cat(paste("NN Cutoff: ", NN_sig, "max cutoff: ", max_sig,"mean cutoff: ", mean_sig, "micron cutoff:", microns,"\nPOSITION DEFECTS (defect / no defect / fractionDefect)",
              "\nExpressing:", ExpPosDef, ExpPosOK, round(ExpPosDef/(ExpPosDef+ExpPosOK),3), 
              "\nNot Expressing:", NotExpPosDef, NotExpPosOK, round(NotExpPosDef/(NotExpPosDef+NotExpPosOK),3), "\n"))



    BothDefectExp = sum(ExpPeak[Cells]>expCutoff & PosDefectCount[Cells] > 0 & (CCDefectCount[sapply(Cells,GetParent)] + CCDefectCount[Cells]) > 0,na.rm=T)
    PosDefectExp = sum(ExpPeak[Cells]>expCutoff & PosDefectCount[Cells] > 0 & (CCDefectCount[sapply(Cells,GetParent)] + CCDefectCount[Cells]) == 0,na.rm=T)
    CCDefectExp = sum(ExpPeak[Cells]>expCutoff & PosDefectCount[Cells] == 0 & (CCDefectCount[sapply(Cells,GetParent)] + CCDefectCount[Cells]) > 0,na.rm=T)
    NoDefectExp = sum(ExpPeak[Cells]>expCutoff & PosDefectCount[Cells] == 0 & (CCDefectCount[sapply(Cells,GetParent)] + CCDefectCount[Cells]) == 0,na.rm=T)

    BothDefectNotExp = sum(ExpPeak[Cells]<=expCutoff & PosDefectCount[Cells] > 0 & (CCDefectCount[sapply(Cells,GetParent)] + CCDefectCount[Cells]) > 0,na.rm=T)
    PosDefectNotExp = sum(ExpPeak[Cells]<=expCutoff & PosDefectCount[Cells] > 0 & (CCDefectCount[sapply(Cells,GetParent)] + CCDefectCount[Cells]) == 0,na.rm=T)
    CCDefectNotExp = sum(ExpPeak[Cells]<=expCutoff & PosDefectCount[Cells] == 0 & (CCDefectCount[sapply(Cells,GetParent)] + CCDefectCount[Cells]) > 0,na.rm=T)
    NoDefectNotExp = sum(ExpPeak[Cells]<=expCutoff & PosDefectCount[Cells] == 0 & (CCDefectCount[sapply(Cells,GetParent)] + CCDefectCount[Cells]) == 0,na.rm=T)


    cat(paste("POS+CC DEFECTS: Both CC and POS defects / only POS / only CC/ neither\n","Exp: Both:", BothDefectExp,"Pos:",PosDefectExp,"CC:",CCDefectExp,"No:",NoDefectExp,
              "\nNo Exp: Both:", BothDefectNotExp,"Pos:",PosDefectNotExp,"CC:",CCDefectNotExp,"No:",NoDefectNotExp),"\n")
    boxplot(
        ExpPeak[Cells[PosDefectCount[Cells] > 0 & (CCDefectCount[sapply(Cells,GetParent)] + CCDefectCount[Cells]) > 0]],
        ExpPeak[Cells[PosDefectCount[Cells] == 0 | 
                         (CCDefectCount[sapply(Cells,GetParent)] + CCDefectCount[Cells]) == 0 & 
                         (PosDefectCount[Cells] > 0 | (CCDefectCount[sapply(Cells,GetParent)] + CCDefectCount[Cells]) > 0)]],
        ExpPeak[Cells[PosDefectCount[Cells] == 0 & (CCDefectCount[sapply(Cells,GetParent)] + CCDefectCount[Cells]) == 0]]
        ,main=paste(i,"expression in cells with 2, 1 or no defect types"))






    ##Alternate approach - using absolute deviation thresholds as well

    PosCells = rownames(theseMeanDevs) 
    print(paste("Number of cells dividing in WT or mutant: ",length(PosCells),
                  " Expressing: ", sum(ExpPeak[PosCells]>=expCutoff,na.rm=T),
                  " Nonexpressing: ", sum(ExpPeak[PosCells]<expCutoff,na.rm=T)))

    CellDefects=theseMeanDevs[Cells,]>microns & theseMeanZ > mean_sig & theseNNs[Cells,] > NN_cutoff & theseNNZ > NN_sig    
    CellDefects[is.na(CellDefects)] <- 0
    rownames(CellDefects) <- Cells
    print(paste("Pos Defects, total: ", 
                  sum(CellDefects)))
    
    ##Defects in expressing cells
    print(paste("Pos Defects, expressing: ", 
                  sum(CellDefects[ExpPeak[Cells]>=expCutoff],na.rm=T)))

    ##Defects in nonexpressing cells
    print(paste("Pos Defects, nonexpressing: ", 
                  sum(CellDefects[ExpPeak[Cells]<expCutoff],na.rm=T)))

    ##Expressing cells defective in one or more embryos
    CellDefectCounts = rowSums(CellDefects,na.rm=T)
    print(paste("Cells with at least one Pos Defect, expressing: ", 
                  sum(CellDefectCounts[ExpPeak[Cells]>=expCutoff]>0,na.rm=T)))
    print(paste("Cells with at least one Pos Defect, nonexpressing: ", 
                  sum(CellDefectCounts[ExpPeak[Cells]<expCutoff]>0,na.rm=T)))

    ##Expressing cells defective in two or more embryos
    print(paste("Cells with at least two Pos Defect, expressing: ", 
                  sum(CellDefectCounts[ExpPeak[Cells]>=expCutoff]>1,na.rm=T)))
    print(paste("Cells with at least two Pos Defect, nonexpressing: ", 
                  sum(CellDefectCounts[ExpPeak[Cells]<expCutoff]>1,na.rm=T)))

    ##list cells defective in 2 or more embryo
    DoubleDefectsExp = Cells[CellDefectCounts>1 & ExpPeak[Cells]>=expCutoff]
    DoubleDefectsNon = Cells[CellDefectCounts>1 & ExpPeak[Cells]<expCutoff]

    
    write.csv(data.frame(
        Cell = Cells,
        Time = CellTimes[Cells],
        Defect_Count = CellDefectCounts[Cells],
        ExpLevel = ExpPeak[Cells],
        Expressing = ExpPeak[Cells]>expCutoff,
        Mean_Dev = theseMeanMean[Cells],
        Max_Dev = theseMaxMean[Cells],
        Mean_NN_score = theseNNMean[Cells],
        Max_NN_score = theseNNMax[Cells]
    ),file=paste(i,"PosDefectSummary.csv",sep="_"),row.names=F)



    ##but now propagate to terminal cells and see if earlier CC defects correlate
    

    #3)?Load orientation devs

    # AngleArray=NULL
    # Angles = read.table(paste(i,"/",i,".angles_aligned.txt",sep=""),header=T)

    # Angles[,1]<- gsub("-", ".", Angles[,1])


    # for(j in Embs){
    #     j2=substr(start=2,stop=nchar(j),x=j)
    #     print(j)
    #     print(j2)
    #     thisCol=as.vector(Angles[Angles[,1]==j2,17])
    #     names(thisCol)=as.vector(Angles[Angles[,1]==j2,2])
    #     AngleArray=cbind(AngleArray,as.numeric(thisCol[Cells]))
    #     }

    # colnames(AngleArray)=Embs
    # rownames(AngleArray)=Cells


    AngleArray=read.csv(paste(i,"_dots.csv",sep=""),header=T,row.names=1)
    AngleArray=AngleArray[Cells,7:length(AngleArray)]
    rownames(AngleArray)=Cells

    divObservedCounts=rowSums(!is.na(AngleArray[Cells,]),na.rm=T)

    
    sum(AngleArray<.5,na.rm=T)

    ##fixme allow more nuanced angle deviation. Maybe better to look at orientation relative to neighbors?
    AngCutoff=.5

    AngOutliers=abs(AngleArray)<AngCutoff

    AngOutliers[AngOutliers==T]<-1
    AngOutliers[AngOutliers==F]<-0

    AngDefectCount=rowSums(AngOutliers,na.rm=T)
    plot(ExpPeak[Cells],jitter(AngDefectCount), xlab="Peak expression", ylab="Angle defet count")
    
    plot(jitter(AngDefectCount),jitter(PosDefectCount),xlab="Angle defect count", ylab="Position defect count")

    plot(jitter(AngDefectCount),jitter(CCDefectCount),xlab="Angle defect count", ylab="CC defect count")
    
    ExpAngDef=sum(ExpPeak[Cells]>expCutoff & AngDefectCount>0,na.rm=T)
    ExpAngOK=sum(ExpPeak[Cells]>expCutoff & PosDefectCount==0,na.rm=T)
    NotExpAngDef=sum(ExpPeak[Cells]<expCutoff & AngDefectCount>0,na.rm=T)
    NotExpAngOK=sum(ExpPeak[Cells]<expCutoff & AngDefectCount==0,na.rm=T)

    cat(paste("ANGLES", "Ang cutoff: ", AngCutoff,
              "\nExpressing:", ExpAngDef, ExpAngOK, round(ExpAngDef/(ExpAngDef+ExpAngOK),3), 
              "\nNot Expressing:", NotExpAngDef, NotExpAngOK, round(NotExpAngDef/(NotExpAngDef+NotExpAngOK),3), "\n"))

    plot(AngDefectCount[Cells],CCDefectCount[Cells], xlab="Angle defect count", ylab="CC defect count",
         main=paste(i,"R=", round(cor(AngDefectCount[Cells],CCDefectCount[Cells],use="pairwise.complete.obs",method="spearman"),digits=4))
         )
    plot(AngDefectCount[Cells],CCDefectCount[sapply(Cells,GetParent)], xlab="Angle defect count", ylab="parent CC defect count",
         main=paste(i,"R=", round(cor(AngDefectCount[Cells],CCDefectCount[sapply(Cells,GetParent)],use="pairwise.complete.obs",method="spearman"), digits=4))
         )
    plot(AngDefectCount[Cells],CCDefectCount[Cells]+CCDefectCount[sapply(Cells,GetParent)],xlab="Angle defect count", ylab="cell plus parent CC defect count",
         main=paste(i,"R=", round(cor(AngDefectCount[Cells],CCDefectCount[Cells]+CCDefectCount[sapply(Cells,GetParent)],use="pairwise.complete.obs",method="spearman"), digits=4))
         )

    plot(AngDefectCount[Cells],PosDefectCount[Cells], xlab="Angle defect count", ylab="Position defect count",
         main=paste(i,"R=", round(cor(AngDefectCount[Cells],PosDefectCount[Cells],use="pairwise.complete.obs",method="spearman"), digits=4))
         )

    plot(AngDefectCount[sapply(Cells,GetParent)],PosDefectCount[Cells],xlab="Parent angle defect count", ylab="Position defect count",
         main=paste(i,"R=", round(cor(AngDefectCount[sapply(Cells,GetParent)],PosDefectCount[Cells],use="pairwise.complete.obs",method="spearman"),digits=4))
         )
    plot(AngDefectCount[Cells]+AngDefectCount[sapply(Cells,GetParent)],PosDefectCount[Cells] , xlab="Cell plus parent angle defect count", ylab="Position defect count",
         main=paste(i,"R=", round(cor(AngDefectCount[Cells]+AngDefectCount[sapply(Cells,GetParent)],PosDefectCount[Cells],use="pairwise.complete.obs",method="spearman"),digits=4))
         )

    #4)Integrate onto list of "terminal cells" matching expression tree

    PosOutliers2=PosOutliers
    AngOutliers2=AngOutliers
    CCDefectArray2=CCDefectArray

    PosOutliers2[is.na(PosOutliers2)]<-0
    AngOutliers2[is.na(AngOutliers2)]<-0
    CCDefectArray2[is.na(CCDefectArray)]<-0

    NumDefects = PosOutliers2+AngOutliers2+CCDefectArray2+CCDefectArray2[sapply(Cells,GetParent)]

    
    #Get # of observations for each defect table
    
    #Import list of cells in lineage order (350 minute sulston cells?)

    #Map "Cells" to these terminal cells
    GetTerminalCell <- function(X){
        #cell is a terminal cell - this works
        if(X %in% CellsLineageOrder[,1]){
            return(X)
            }

        #cell is a descendent of a terminal cell - this works
        tmpCell=X
        while(!is.na(GetParent(tmpCell)) && GetParent(tmpCell) != "P"){
            tmpCell=GetParent(tmpCell)
            if(tmpCell %in% CellsLineageOrder[,1]){
                return(tmpCell)
                }
            }

        #cell is an ancestor of a terminal cell - works but sends out a list of cells
        Descendants=NULL
        for(k in CellsLineageOrder[,1]){
            tmpCell=k
            while(!is.na(GetParent(tmpCell)) && GetParent(tmpCell) != "P"){
                tmpCell=GetParent(tmpCell)
                if(tmpCell == X){
                    Descendants=c(Descendants,k)
                    }
                }
            }
        if(length(Descendants)>0){
            return(Descendants)
            }
        return(NA)
        }
    #NewCells=sapply(Cells,GetTerminalCell)
    #Output plots of defect frequency
    
    AncestorDefects<-function(Cell,DefectArray){
        MyDefects=0
        if(Cell %in% rownames(DefectArray)){
            MyDefects=sum(DefectArray[Cell,],na.rm=T)
            }
        Parent=GetParent(Cell)
        if(Parent!="P" && !is.na(Parent)){
            MyDefects = MyDefects + AncestorDefects(Parent,DefectArray)
            }
        return(MyDefects)

        }

    SumDefects<-function(X,DefectArray){
        return(sum(DefectArray[grep(paste(X,".",sep=""),rownames(DefectArray)),],na.rm=T) +
              AncestorDefects(X,DefectArray))
        
        }
    CCTermDefects=sapply(as.vector(CellsLineageOrder[,1]),function(X){SumDefects(X,CCDefectArray)})
    names(CCTermDefects)=CellsLineageOrder[,1]
    
    PosTermDefects=sapply(as.vector(CellsLineageOrder[,1]),function(X){SumDefects(X,PosOutliers)})
    names(PosTermDefects)=CellsLineageOrder[,1]

    AngTermDefects=sapply(as.vector(CellsLineageOrder[,1]),function(X){SumDefects(X,AngOutliers)})
    names(AngTermDefects)=CellsLineageOrder[,1]


    plot(CCTermDefects,ylim=c(min(PosTermDefects,AngTermDefects,CCTermDefects),max(PosTermDefects,AngTermDefects,CCTermDefects)),xlab=CellsLineageOrder[,1],xaxt="n",main=paste(i,"black=CC red=Pos green=Ang"), cex=0.1)
    axis(1,at=1:length(CellsLineageOrder[,1]),labels=CellsLineageOrder[,1],las=2, ylab="Count, includes ancestors")
    points(PosTermDefects+.05,col=2, cex=0.1)
    points(AngTermDefects-.05,col=3, cex=0.1)

    plot(CCTermDefects,main=i, ylab="CC defect count including ancestors", xlab="Terminal cells, lineage order", cex=0.1)
    plot(PosTermDefects,main=i, ylab="Position defect count including ancestors", xlab="Terminal cells, lineage order", cex=0.1)
    plot(AngTermDefects,main=i, ylab="Angle defect count including ancestors", xlab="Terminal cells, lineage order", cex=0.1)
    plot(CCTermDefects+PosTermDefects+AngTermDefects,main=i , ylab="Total defect count including ancestors", xlab="Terminal cells, lineage order", cex=0.1)

#Fixme generalize
    for(Founder in c("ABalpa", "ABara", "ABplp", "ABprp", "MSaa", "MSpa")){
        plot(CCTermDefects[grep(Founder,names(CCTermDefects))],main=paste(i,Founder), ylab="CC defect count", xlab=paste(Founder, "cells lineage order"))
#             labels= names(CCTermDefects)[grep(Founder,names(CCTermDefects))])
        plot(PosTermDefects[grep(Founder,names(CCTermDefects))],main=paste(i,Founder), ylab="Position defect count", xlab=paste(Founder, "cells lineage order"))
#             labels= names(CCTermDefects)[grep(Founder,names(CCTermDefects))])
        plot(AngTermDefects[grep(Founder,names(CCTermDefects))],main=paste(i,Founder), ylab="Angle defect count", xlab=paste(Founder, "cells lineage order"))
#             labels= names(CCTermDefects)[grep(Founder,names(CCTermDefects))])

        plot(CCTermDefects[grep(Founder,names(CCTermDefects))],ylim=c(min(PosTermDefects,AngTermDefects,CCTermDefects),max(PosTermDefects,AngTermDefects,CCTermDefects)),xlab=CellsLineageOrder[,1],xaxt="n",main=paste(i, Founder, "black=CC red=Pos green=Ang"))
#             labels= names(CCTermDefects)[grep(Founder,names(CCTermDefects))])

        points(PosTermDefects[grep(Founder,names(CCTermDefects))],col=2)
        points(AngTermDefects[grep(Founder,names(CCTermDefects))],col=3)


    }
    Parents=sapply(Cells,GetParent)
    Sisters=sapply(Cells,GetSister)
    Aunts=sapply(Parents,GetSister)


    write.csv(data.frame(
        Cells=CellsLineageOrder[,1],
        Fates=CellNames[as.vector(CellsLineageOrder[,1])],
        CCTermDefects,
        PosTermDefects,
        AngTermDefects,
        Expression=ExpPeak[CellsLineageOrder[,1]],
        ceh36_peak[CellsLineageOrder[,1]],
        unc30_peak[CellsLineageOrder[,1]],
        mls2_peak[CellsLineageOrder[,1]],
        nob1_peak[CellsLineageOrder[,1]],
        ceh13_peak[CellsLineageOrder[,1]]
    ),file=paste(i,"Term_defects.csv",sep="_"),row.names=F)

    write.csv(
        data.frame(Cells,
                   Fates=CellNames[Cells],
                   Times=CellTimes[Cells],
                   CC_Defects=CCDefectCount[Cells],
                   Pos_Defects=PosDefectCount[Cells],
                   Angle_Defects=AngDefectCount[Cells],
                   Cell_Count=cellObservedCounts[Cells],
                   Division_Count=divObservedCounts[Cells],
                   Expression=ExpPeak[Cells],
                   ceh36_peak[Cells],
                   unc30_peak[Cells],
                   mls2_peak[Cells],
                   nob1_peak[Cells],
                   ceh13_peak[Cells],
                   MeanMaxPositionDev=theseMaxMean[Cells],
                   MeanMeanPositionDev=theseMeanMean[Cells]),
        file=paste(i,"_defect_counts.csv",sep="_"),row.names=F)


    write.csv(data.frame(Cells,
                         CellNames[Cells],
                         CellTimes[Cells],
                         CCDefectCount[Cells],
                         CCDefectCount[Parents],
                         CCDefectCount[Sisters],
                         CCDefectCount[Aunts],
                         PosDefectCount[Cells],
                         PosDefectCount[Parents],
                         PosDefectCount[Sisters],
                         PosDefectCount[Aunts],
                         AngDefectCount[Cells],
                         AngDefectCount[Parents],
                         AngDefectCount[Sisters],
                         AngDefectCount[Aunts],
                         cellObservedCounts[Cells],
                         cellObservedCounts[Parents],
                         cellObservedCounts[Sisters],
                         cellObservedCounts[Aunts],
                         divObservedCounts[Cells],
                         divObservedCounts[Parents],
                         divObservedCounts[Sisters],
                         divObservedCounts[Aunts],
                         ExpPeak[Cells],
                         ceh36_peak[Cells],
                         unc30_peak[Cells],
                         mls2_peak[Cells],
                         nob1_peak[Cells],
                         ceh13_peak[Cells]),
              file=paste(i,"defect_counts_with_relatives.csv",sep="_"),row.names=F)

    write.csv(data.frame(Cells,CellNames[Cells],CellTimes[Cells],
                         CCDefectArray,
                         AncestralDefects,
                         PosOutliers),
              file=paste(i,"defects_with_ancestral_CC.csv",sep="_"),row.names=F)
                         

    tmp=data.frame(Cells,CellNames[Cells],CellTimes[Cells],ExpPeak[Cells],CCDefectArray[Cells,],PosOutliers[Cells,],AngOutliers[Cells,])
    colnames(tmp)=c("Cells","Names", "Birth Times", "peak exp", paste("CC",colnames(PosOutliers)),paste("Position",colnames(PosOutliers)),paste("Angle",colnames(PosOutliers)))
    write.csv(tmp,file=paste(i,"defects.csv",sep="_"),row.names=F)
    

    dev.off()
    setwd(startDir)
}






# tmp=aggregate(. ~ cell,  data=read.csv("ACDunc-30_mean.csv",row.names=1), mean, na.rm=T)
# unc30_exp=data.frame(tmp[,"blot"])
# rownames(unc30_exp)=tmp[,"cell"]
# unc30_peak=sapply(Cells,function(X){
#         GetPeak(unc30_exp,X)
#         })  

# tmp=aggregate(. ~ cell,  data=read.csv("ACDceh-36_mean.csv",row.names=1), mean, na.rm=T)
# ceh36_exp=data.frame(tmp[,"blot"])
# rownames(ceh36_exp)=tmp[,"cell"]
# ceh36_peak=sapply(Cells,function(X){
#         GetPeak(ceh36_exp,X)
#         })  

# tmp=aggregate(. ~ cell,  data=read.csv("ACDmls-2_mean.csv",row.names=1), mean, na.rm=T)
# mls2_exp=data.frame(tmp[,"blot"])
# rownames(mls2_exp)=tmp[,"cell"]
# mls2_peak=sapply(Cells,function(X){
#         GetPeak(mls2_exp,X)
#         })  

# tmp=aggregate(. ~ cell,  data=read.csv("ACDnob-1GFP_mean.csv",row.names=1), mean, na.rm=T)
# nob1_exp=data.frame(tmp[,"blot"])
# rownames(nob1_exp)=tmp[,"cell"]
# nob1_peak=sapply(Cells,function(X){
#         GetPeak(nob1_exp,X)
#         })  


