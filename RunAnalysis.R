#Can comment this out if already in workspace
#print("Running ProcessError")
#source("ProcessError.R")
#print("Running AnalyzeDivTime")
#source("AnalyzeWTDivTime.R")
#print("Running ProcessAngles")
#source("ProcessAngles.R")
#print("Running AnalyzePosition")
#source("AnalyzePosition.R")

load("Model.RData")

#BaseFilename="test"
#BaseFilename="26degrees_all"
#BaseFilename="Over500WT"
#BaseFilename="ceh-36ok795.txt"
#BaseFilename="24_25degrees"
WTFilename="Over500WT"
#BaseFilename="Over500All"

source("GetList.R")

WTPositionFileName=paste(BaseFilename,"/",BaseFilename,"positions.txt",sep="")

PositionFileName=paste(BaseFilename,"/",BaseFilename,"positions.txt",sep="")

TestPositions=read.table(PositionFileName, header=T,stringsAsFactors=F,sep="\t")
TestPositions[,"concat"] = paste(TestPositions[,1],TestPositions[,2],sep=":")
rownames(TestPositions)=paste(TestPositions[,1],TestPositions[,2],sep=":")

#calculate mean positions and compare by cell?
#or try calculating mean deviation for each cell
 #align positions by concat
 #

compPositions=TestPositions[rownames(MeanPositionFrame),]

tbaseIndeces=1:((length(compPositions[1,])-2)/3)
tXindeces=3*tbaseIndeces
tYindeces=3*tbaseIndeces+1
tZindeces=3*tbaseIndeces+2

tXpositions=compPositions[,tXindeces]
tYpositions=compPositions[,tYindeces]
tZpositions=compPositions[,tZindeces]

tXdev=(tXpositions-MeanPositionFrame[,3])^2
tYdev=(tYpositions-MeanPositionFrame[,4])^2
tZdev=(tZpositions-MeanPositionFrame[,5])^2

#deviations by cell time point
tDeviations=sqrt(tXdev+tYdev+tZdev)


#mean deviations by cell
cellLookups=data.frame(rownames(MeanPositionFrame),MeanPositionFrame[,1],row.names=1)
tDeviations=data.frame(cellLookups,tDeviations)
colnames(tDeviations)[1]="Cell"
tCellDeviations=aggregate(tDeviations,tDeviations["Cell"],mean)
rownames(tCellDeviations)=tCellDeviations[,1]

pdf(paste(BaseFilename,"PostionDeviations.pdf"),width=30)

MeanDeviationsPerCell=(apply(tCellDeviations[,2:length(tCellDeviations)],1,mean,na.rm=T))
boxplot(MeanDeviationsPerCell,xlab=BaseFilename,ylab="Mean Deviation(microns)")
plot(CellTmeans[names(MeanDeviationsPerCell)],MeanDeviationsPerCell,xlab="Mean T(minutes)",ylab="Mean Deviation(microns)")

write.csv(tCellDeviations,paste(BaseFilename,"tCellDeviations.csv"))
for(i in colnames(tCellDeviations)){
    plot(tCellDeviations[,i],xlab="cells",ylab=paste("Deviation",i),ylim=c(0,30),xaxt="n")


    axis(1,1:length(tCellDeviations[,i]),labels=tCellDeviations[,1])
    plot(CellTmeans[tCellDeviations[,1]],tCellDeviations[,i],xlab="time",ylab=paste("Deviation",i),ylim=c(0,30))
}
dev.off()

#Same thing for just x deviations
tXDeviations=sqrt(tXdev)
tXDeviations=data.frame(cellLookups,tXDeviations)
colnames(tXDeviations)[1]="Cell"
tXCellDeviations=aggregate(tXDeviations,tXDeviations["Cell"],mean)
rownames(tXCellDeviations)=tXCellDeviations[,1]

pdf(paste(BaseFilename,"XPositionDeviations.pdf"),width=30)

MeanXDeviationsPerCell=(apply(tXCellDeviations[,2:length(tXCellDeviations)],1,mean,na.rm=T))
boxplot(MeanXDeviationsPerCell,xlab=BaseFilename,ylab="Mean X Deviation(microns)")
plot(CellTmeans[names(MeanXDeviationsPerCell)],MeanXDeviationsPerCell,xlab="Mean T(minutes)",ylab="Mean X Deviation(microns)")
plot(abs(dX[names(MeanXDeviationsPerCell)]),MeanXDeviationsPerCell)
cor(abs(dX[names(MeanXDeviationsPerCell)]),MeanXDeviationsPerCell)

MeanXDeviationsPerLateCell=MeanXDeviationsPerCell[CellTmeans[names(MeanXDeviationsPerCell)]>200]
plot(MeanXDeviationsPerLateCell,abs(dX[names(MeanXDeviationsPerLateCell)]))
cor(MeanXDeviationsPerLateCell,abs(dX[names(MeanXDeviationsPerLateCell)]))

#sort(MeanXDeviationsPerCell[CellTmeans[names(MeanXDeviationsPerCell)]>200])

write.csv(tXCellDeviations,paste(BaseFilename,"tXCellDeviations.csv"))
for(i in colnames(tXCellDeviations)){
    plot(tXCellDeviations[,i],xlab="cells",ylab=paste("X Deviation",i),ylim=c(0,30),xaxt="n")

    axis(1,1:length(tXCellDeviations[,i]),labels=tXCellDeviations[,1])
    plot(CellTmeans[tXCellDeviations[,1]],tXCellDeviations[,i],xlab="time",ylab=paste("X Deviation",i),ylim=c(0,30))
}
dev.off()




TimeFileName=paste(BaseFilename,"/",BaseFilename,"CCLengthNorm.tsv",sep="")
TestTimes=read.table(TimeFileName, header=T,stringsAsFactors=F,sep="\t",row.names=1)
#TestTimes=TestTimes[,2:length(TestTimes)]
compTimes=TestTimes[names(CellAverages),]


DivFileName=paste(BaseFilename,"/",BaseFilename,"DivTimeNorm.tsv",sep="")
DivTestTimes=read.table(DivFileName, header=T,stringsAsFactors=F,sep="\t",row.names=1)
#DivTestTimes=DivTestTimes[,2:length(DivTestTimes)]
theseDivTimes=DivTestTimes[names(CellAverages),]



pdf(paste(BaseFilename,"CC_Length_Deviations.pdf",sep=""))
for(i in colnames(compTimes)){
    plot(theseDivTimes[,i],DivTimeAverages[names(CellAverages)],xlab=i,ylab="WT AVERAGE DIV TIME")
    plot(compTimes[,i],CellAverages,xlab=i,ylab="WT AVERAGE CC")
    plot(compTimes[,i]/CellAverages,ylim=c(0,4),xaxt="n",ylab=paste("CC RATIO TO WT",i))
    axis(1,1:length(compTimes[,i]),labels=rownames(compTimes))
}
dev.off()

plot(theseDivTimes[,i],compTimes[,i])

#ADD ECTOPIC AND MISSED DIVISIONS.  TRY Z SCORES

write.csv(data.frame(CellAverages[rownames(compTimes)],DivTimeAverages[rownames(compTimes)],CC_SD[rownames(compTimes)],theseDivTimes,compTimes,compTimes-CellAverages[rownames(compTimes)]),paste(BaseFilename,"CC_Deviations.csv"))


MissedFileName=paste(BaseFilename,"/",BaseFilename,"CCLengthMinTerminal.tsv",sep="")
TestTerminalBranches=divTimes=read.table(MissedFileName, header=T,row.names=1,stringsAsFactors=F)
#TestTerminalBranches=TestTerminalBranches[,2:length(TestTerminalBranches)]

#count >3sigma failed to divide branches missedDivisionsCount=apply(TerminalBranches>CellAverages+3*CC_SD,2,sum,na.rm=T)
TestMissedDivisionsCount=apply(TestTerminalBranches>(CellAverages+3*CC_SD),2,sum,na.rm=T)
TestEarlyDivisionsCount=apply(TestTimes<(CellAverages-3*CC_SD),2,sum,na.rm=T)
TestLateDivisionsCount=apply(TestTimes>(CellAverages+3*CC_SD),2,sum,na.rm=T)

TestCC_Zscores=(TestTimes-CellAverages)/CC_SD
TestCC_Deltas=TestTimes-CellAverages

TestUnscoredDivisionMinDeltas=TestTerminalBranches-CellAverages
TestUnscoredDivisionMinZ=(TestTerminalBranches-CellAverages)/CC_SD

TestMissedDivisionsCountCells=apply(TestTerminalBranches>(CellAverages+3*CC_SD),1,sum,na.rm=T)
TestEarlyDivisionsCountCells=apply(TestTimes<(CellAverages-3*CC_SD),1,sum,na.rm=T)
TestLateDivisionsCountCells=apply(TestTimes>(CellAverages+3*CC_SD),1,sum,na.rm=T)
TestVeryEarlyDivisionsCountCells=apply(TestTimes<(CellAverages-3*CC_SD)&TestTimes<(CellAverages-10),1,sum,na.rm=T)
TestVeryLateDivisionsCountCells=apply(TestTimes>(CellAverages+3*CC_SD)&TestTimes>(CellAverages+10),1,sum,na.rm=T)

write.csv(data.frame(DivTimeAverages[names(TestMissedDivisionsCountCells)],TestMissedDivisionsCountCells,TestEarlyDivisionsCountCells,TestLateDivisionsCountCells,TestVeryEarlyDivisionsCountCells,TestVeryLateDivisionsCountCells),paste(BaseFilename,"DivisionOutliers.csv"))

pdf(paste(BaseFilename,'Outliers.pdf',sep=""),width=30)
boxplot(TestCC_Zscores,xlab="Embryos",ylab="Cell cycle Z scores",xaxt="n")
axis(1,1:length(TestCC_Zscores),labels=colnames(TestCC_Zscores))

boxplot(TestCC_Deltas,xlab="Embryos",ylab="Cell cycle offset (minutes)",xaxt="n")
axis(1,1:length(TestCC_Zscores),labels=colnames(TestCC_Zscores))

boxplot(TestUnscoredDivisionMinDeltas,xlab="Embryos",ylab="Failed to divide - observed length delta (minutes)",xaxt="n")
axis(1,1:length(TestCC_Zscores),labels=colnames(TestCC_Zscores))

boxplot(TestUnscoredDivisionMinZ,xlab="Embryos",ylab="Failed to divide - observed length Z score",xaxt="n")
axis(1,1:length(TestCC_Zscores),labels=colnames(TestCC_Zscores))

plot(TestMissedDivisionsCount)
plot(TestEarlyDivisionsCount)
plot(TestLateDivisionsCount)

plot(DivTimeAverages[names(TestLateDivisionsCountCells)],TestMissedDivisionsCountCells)
plot(DivTimeAverages[names(TestLateDivisionsCountCells)],TestEarlyDivisionsCountCells)
plot(DivTimeAverages[names(TestLateDivisionsCountCells)],TestLateDivisionsCountCells)
plot(DivTimeAverages[names(TestLateDivisionsCountCells)],TestVeryEarlyDivisionsCountCells)
plot(DivTimeAverages[names(TestLateDivisionsCountCells)],TestVeryLateDivisionsCountCells)

plot(TestMissedDivisionsCountCells,xaxt="n")
points(TestVeryEarlyDivisionsCountCells,col=2)
points(TestVeryLateDivisionsCountCells,col=3)
axis(1,1:length(TestMissedDivisionsCountCells),labels=names(TestMissedDivisionsCountCells))

Events=TestMissedDivisionsCountCells+TestVeryEarlyDivisionsCountCells+TestVeryLateDivisionsCountCells
Events
dev.off()

plot(tXCellDeviations[rownames(TestCC_Deltas),3],TestCC_Deltas[,1])




#Base stats, rate, lineage-specific rates

#align CC lengths by cell
#plot CC vs WT
#compute z scores, deviations, plot, rank


#align division angle by cell
#plot deviation angle by cell (and z scores)

#TestAngles=read.table("Over500ALL/Over500ALL.angles_unaligned.txt", header=T,stringsAsFactors=F,sep="\t")
