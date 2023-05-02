source("functions.R")
library(rgl)
library(gdata)
library(plotrix)

WTDivTimes=read.table("Richard_et_al_plus_comma_WT/Richard_et_al_plus_comma_WTDivTimeNorm.tsv", header=T,row.names=1,stringsAsFactors=F)
    #Remove sulston which is included in these tables for historical reasons
WTDivTimes[,1] <- NULL
WTCellCycles=read.table("Richard_et_al_plus_comma_WT/Richard_et_al_plus_comma_WTCCLengthNorm.tsv", header=T,row.names=1,stringsAsFactors=F)
WTCellCycles[,1]<-NULL

Cells=rownames(WTDivTimes)

    #calculate mean and stdev of WT CC and div times


WTDivMeans=rowMeans(WTDivTimes,na.rm=T)
WTDivSDs=apply(WTDivTimes,1,sd,na.rm=T)
WTCCMeans=rowMeans(WTCellCycles,na.rm=T)
WTCCSDs=apply(WTCellCycles,1,sd,na.rm=T)

WTCCCounts=rowSums(!is.na(WTCellCycles))
WTDivCounts=rowSums(!is.na(WTDivTimes))

## Require at least 2 observations in WT
WTDivMeans[WTDivCounts<2]=NA
WTDivSDs[WTDivCounts<2]=NA

WTCCMeans[WTCCCounts<2]=NA
WTCCSDs[WTCCCounts<2]=NA

AnalyzeDivTimes <- function(Name){
    message(paste("Analyzing division times", Name))
    startDir=getwd()
    setwd(Name)    
    #calculate mean and stdev of 'mutant' CC and div times
    MutantDivTimes=read.table(paste(Name,"DivTimeNorm.tsv",sep=""), header=T,row.names=1,stringsAsFactors=F)
    MutantDivTimes[,1]<-NULL
    MutantCellCycles=read.table(paste(Name,"CCLengthNorm.tsv",sep=""), header=T,row.names=1,stringsAsFactors=F)
    MutantCellCycles[,1]<-NULL
    MutantTerminalLengths=read.table(paste(Name,"CCLengthMinTerminal.tsv",sep=""), header=T,row.names=1,stringsAsFactors=F)


    MutantDivMeans=rowMeans(MutantDivTimes,na.rm=T)
    MutantDivSDs=apply(MutantDivTimes,1,sd,na.rm=T)
    MutantCCMeans=rowMeans(MutantCellCycles,na.rm=T)
    MutantCCSDs=apply(MutantCellCycles,1,sd,na.rm=T)

    #Make some useful lists 
    Cells=union(rownames(MutantDivTimes),rownames(WTDivTimes))
    ABpxpCells=Cells[grep("^ABp[lr]p",Cells)]
    MSxaCells=Cells[grep("^MS[ap]a",Cells)]
    ABara_alpaCells=union(Cells[grep("^ABalpa",Cells)],Cells[grep("^ABara",Cells)])
    OtherCells=Cells[-grep("^ABp[lr]p",Cells)]
    WTEmbryos=colnames(WTCellCycles)
    MutantEmbryos=colnames(MutantCellCycles)
    ABalaCells=Cells[grep("^ABala",Cells)]
    ABalpCells=Cells[grep("^ABalp",Cells)]
    ABaraCells=Cells[grep("^ABara",Cells)]
    ABarpCells=Cells[grep("^ABarp",Cells)]
    ABplaCells=Cells[grep("^ABpla",Cells)]
    ABplpCells=Cells[grep("^ABplp",Cells)]
    ABpraCells=Cells[grep("^ABpra",Cells)]
    ABprpCells=Cells[grep("^ABprp",Cells)]
    MSCells=Cells[grep("^MS",Cells)]
    ECells=Cells[grep("^E",Cells)]
    CCells=Cells[grep("^C",Cells)]
    DCells=Cells[grep("^D",Cells)]
    PCells=Cells[grep("^P",Cells)]
    pdf(paste(Name,"CC_plots.pdf",sep=""))

    # THIS PLOTS DIV TIME FOR EACH MUTANT SERIES VS WT AVERAGE
    for(i in MutantEmbryos){
        plot(MutantDivTimes[OtherCells,i],WTDivMeans[OtherCells],xlab=i,main="Div Time ABpxp: red, MSxa: green, ABara/alpa: blue, Other: black")
        points(MutantDivTimes[ABpxpCells,i],WTDivMeans[ABpxpCells],col=2)
        points(MutantDivTimes[MSxaCells,i],WTDivMeans[MSxaCells],col=3)
        points(MutantDivTimes[ABara_alpaCells,i],WTDivMeans[ABara_alpaCells],col=4)
        }
    # THIS PLOTS CC LENGTH FOR EACH MUTANT SERIES VS WT AVERAGE
    for(i in MutantEmbryos){
        plot(MutantCellCycles[OtherCells,i],WTCCMeans[OtherCells],xlab=i,main="CC Length ABpxp: red, MSxa: green, ABara/alpa: blue, Other: black")
        points(MutantCellCycles[ABpxpCells,i],WTCCMeans[ABpxpCells],col=2)
        points(MutantCellCycles[MSxaCells,i],WTCCMeans[MSxaCells],col=3)
        points(MutantCellCycles[ABara_alpaCells,i],WTCCMeans[ABara_alpaCells],col=4)
        }
    # This calculates and plots p-value for whether mean of each cell cycle length is significantly different between WT and mutant
    CC_Dif_pValues=sapply(Cells, function(X){
            obj<-try(t.test(WTCellCycles[X,],MutantCellCycles[X,]),silent=T)
            if (is(obj, "try-error")) return(NA) else return(obj$p.value)})
    plot(WTDivMeans[OtherCells],log(CC_Dif_pValues[OtherCells]), xlab="WT DIV TIME, min", ylab = "log(10) p Value; mean Mutant CC == mean WT CC",
         main="CC Length ABpxp: red, MSxa: green, ABara/alpa: blue, Other: black")   
    points(WTDivMeans[ABpxpCells],log(CC_Dif_pValues[ABpxpCells]),col=2)
    points(WTDivMeans[MSxaCells],log(CC_Dif_pValues[MSxaCells]),col=3)
    points(WTDivMeans[ABara_alpaCells],log(CC_Dif_pValues[ABara_alpaCells]),col=4)


    AlphOrder=1:length(Cells)
    names(AlphOrder)=Cells

    plot(AlphOrder[Cells],log(CC_Dif_pValues[Cells]),xlab="cells alph order",main="red:ABaxa,green:ABaxp,blue:ABpxa,cyan:ABpxp,black:C,red:D,grey:E,yellow:MSxa,magenta:MSxp,green:P")
    points(AlphOrder[ABalaCells],log(CC_Dif_pValues[ABalaCells]),xlab="cells alph order",col=2)
    points(AlphOrder[ABaraCells],log(CC_Dif_pValues[ABaraCells]),xlab="cells alph order",col=2)
    points(AlphOrder[ABalpCells],log(CC_Dif_pValues[ABalpCells]),xlab="cells alph order",col=3)
    points(AlphOrder[ABarpCells],log(CC_Dif_pValues[ABarpCells]),xlab="cells alph order",col=3)
    points(AlphOrder[ABplaCells],log(CC_Dif_pValues[ABplaCells]),xlab="cells alph order",col=4)
    points(AlphOrder[ABpraCells],log(CC_Dif_pValues[ABpraCells]),xlab="cells alph order",col=4)
    points(AlphOrder[ABpxpCells],log(CC_Dif_pValues[ABpxpCells]),xlab="cells alph order",col=5)
    points(AlphOrder[MSCells],log(CC_Dif_pValues[MSCells]),xlab="cells alph order",col=6)
    points(AlphOrder[MSxaCells],log(CC_Dif_pValues[MSxaCells]),xlab="cells alph order",col=7)
    points(AlphOrder[ECells],log(CC_Dif_pValues[ECells]),xlab="cells alph order",col=8)
    points(AlphOrder[CCells],log(CC_Dif_pValues[CCells]),xlab="cells alph order",col=9)
    points(AlphOrder[DCells],log(CC_Dif_pValues[DCells]),xlab="cells alph order",col=10)
    points(AlphOrder[PCells],log(CC_Dif_pValues[PCells]),xlab="cells alph order",col=11)


    #Generate and output table of embryo/cell deviations
    CCdevs=data.frame()
    for(i in MutantEmbryos){
        CCdevs=rbind(CCdevs,data.frame(Cells,i,MutantDivTimes[Cells,i],WTDivMeans[Cells],WTDivSDs[Cells],MutantCellCycles[Cells,i],WTCCMeans[Cells],WTCCSDs[Cells],
                                       MutantCellCycles[Cells,i]-WTCCMeans[Cells],(MutantCellCycles[Cells,i]-WTCCMeans[Cells])/WTCCSDs[Cells],
                                       MutantTerminalLengths[Cells,i],MutantTerminalLengths[Cells,i]-WTCCMeans[Cells],(MutantTerminalLengths[Cells,i]-WTCCMeans[Cells])/WTCCSDs[Cells],
                                       row.names=paste(Cells,i,sep="_")))
        }
    colnames(CCdevs)<-c("Cell","Embryo","Div time", "WT Div time", "WT Div time SD", "CC", "WT CC", "WT CC SD", "CC Deviation", "CC Z score", 
                        "Terminal Cell Length", "Terminal cell delta", "Terminal Cell Z")

    #Output full list, and filtered list containing just potentially interesting deviant cells
    write.csv(CCdevs,file=paste(Name,"_ccDevs.csv",sep=""))
    write.csv(
        subset(CCdevs,
               #Filter for CC deviation > 5 and z score > 3
               ((abs(CCdevs[,9])>5&abs(CCdevs[,10])>3)|
               #Filter for cells that should but don't divide and have observed lengths > 5min more than WT CC and z score > 3
               (CCdevs[,12]>5&CCdevs[,13]>3))|
               #Filter for cells that divide inappropriately
               (CCdevs[,3]>0&is.nan(CCdevs[,4]))),
        file=paste(Name,"5min_3sd_Devs.csv",sep=""))

    #Plot CC deviation versus Z score
    plot(CCdevs[is.element(CCdevs[,1],OtherCells),9],CCdevs[is.element(CCdevs[,1],OtherCells),10],xlab="CC Deviation", ylab="CC z score")
    points(CCdevs[is.element(CCdevs[,1],ABpxpCells),9],CCdevs[is.element(CCdevs[,1],ABpxpCells),10],col=2)
    points(CCdevs[is.element(CCdevs[,1],MSxaCells),9],CCdevs[is.element(CCdevs[,1],MSxaCells),10],col=3)
    points(CCdevs[is.element(CCdevs[,1],ABara_alpaCells),9],CCdevs[is.element(CCdevs[,1],ABara_alpaCells),10],col=4)


    #plot CC Z score vs timeW
    plot(CCdevs[,3],rowMeans(data.frame(CCdevs[,10],CCdevs[,13]),na.rm=T),xlab="Division time", ylab="z score for CC deviation")

    #plot CC Z score vs lineage
    plot(AlphOrder[CCdevs[,1]],rowMeans(data.frame(CCdevs[,10],CCdevs[,13]),na.rm=T),xlab="Cells alph order",ylab="z score for CC deviation")
    
    dev.off()

    # OUTPUTS CLUSTERABLE DEVIATION AND CC LENGTH SPREADSHEETS
    write.table(data.frame(paste(Cells,round(CC_Dif_pValues,4)),round(data.frame(WTCellCycles[Cells,],MutantCellCycles[Cells,]))),file=paste(Name,"CC.txt",sep=""),sep="\t",na="",quote=F,row.names=F)
    write.table(data.frame(paste(Cells,round(CC_Dif_pValues,4)),round(data.frame(WTCellCycles[Cells,]-WTCCMeans[Cells],MutantCellCycles[Cells,]-WTCCMeans[Cells]))),file=paste(Name,"CCdev.txt",sep=""),sep="\t",na="",quote=F,row.names=F)

    #Count and output table of outliers per cell
    LateCalls=sapply(Cells, function(X) {
            test=CCdevs[CCdevs[,1]==X&CCdevs[,9]>5&CCdevs[,10]>3,]
            sum(!is.na(test[3]))})

    EarlyCalls=sapply(Cells, function(X) {
            test=CCdevs[CCdevs[,1]==X&CCdevs[,9]<(-5)&CCdevs[,10]<(-3),]
            sum(!is.na(test[3]))})
    
    MissedCalls=sapply(Cells, function(X) {
            test=CCdevs[CCdevs[,1]==X&CCdevs[,12]>5&CCdevs[,13]>3,]
            sum(!is.na(test[4]))})

    EctopicCalls=sapply(Cells, function(X) {
            test=CCdevs[CCdevs[,1]==X&is.na(CCdevs[,4])&!is.na(CCdevs[,3]),]
            sum(!is.na(test[3]))})
    TotalDefects=LateCalls+EarlyCalls+MissedCalls+EctopicCalls

    write.csv(data.frame(Cells,WTDivMeans[Cells],LateCalls,EarlyCalls,MissedCalls,EctopicCalls,TotalDefects,CC_Dif_pValues,WTCCMeans[Cells],MutantCCMeans[Cells]),file=paste(Name,"CellDefectSummary.csv",sep=""),quote=F,row.names=F)
    setwd(startDir)

    ##quantitative deviations
    Deviations=MutantCellCycles[Cells,]-WTCCMeans[Cells]

    ##Ectopic divisions - enter CC and set to zero if not observed (conservative on inclusion but may be too high level
    #set divisions not observed in WT as zero
#    Deviations[is.na(WTCCMeans[Cells]),]=0
    #set subset of those observed in mutants to mutant CC (consider instead making this a small # based on scale of CC devs ~10 min or so?)
#    Deviations[!is.na(MutantCellCycles[Cells,]) & is.na(WTCCMeans[Cells])] = MutantCellCycles[!is.na(MutantCellCycles[Cells,])& is.na(WTCCMeans[Cells])]
    Deviations[!is.na(MutantCellCycles[Cells,]) & is.na(WTCCMeans[Cells])] = -50

    ##non-observed divisions (conservative evidence for late)
    NonObservedMinDevs=MutantTerminalLengths[Cells,2:length(MutantTerminalLengths)]-WTCCMeans[Cells]
    Deviations[NonObservedMinDevs>0 & !is.na(NonObservedMinDevs)]=NonObservedMinDevs[NonObservedMinDevs>0 & !is.na(NonObservedMinDevs)]

    return(Deviations)
    
}
