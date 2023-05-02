#setwd("Over500All")
#library("gplots")
source("functions.R")
source("PlotDefectSummaries.R")
source("PlotPositionDevs.R")
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


#source("GetPeak.R")
#This should be in functions.R now



ceh36_exp_orig=read.csv("data/CA20120714_JIM136_L3.csv",row.names=2)["blot"]
ceh36_peak_orig=sapply(Cells,function(X){
        GetPeak(ceh36_exp_orig,X)
        })  

nhr67_peak=ReadPeakExpression("data/CA20130527_nhr-67_JIM202_L2.csv")
ceh36_peak=ReadPeakExpression("data/CA20120714_JIM136_L3.csv")
mls2_peak=ReadPeakExpression("data/CA20110209_UP2051_mls-2_L2.csv")
ceh13_peak=ReadPeakExpression("data/CA20120305_ceh-13_L3.csv")
unc30_peak=ReadPeakExpression("data/CA20121115_unc-30_18_L2.csv")
nob1_peak=ReadPeakExpression("data/CA20140825_nob-1_JIM284_L3.csv")
#ref2_peak=ReadPeakExpression("data/CA20141119_ref-2_OP520_L2.csv")
ref2_peak=ReadPeakExpression("data/CA20160203_ref-2_L2.csv")
elt1_peak=ReadPeakExpression("data/CA20120130_elt-1_17_L2.csv")
unc130_peak=ReadPeakExpression("data/CA20110927_RW11144_L4.csv")
php3_peak=ReadPeakExpression("data/CA20201007_JIM644_L1.csv")
ceh32_peak=ReadPeakExpression("data/CA_CEH-32_SYS85.csv")
tbx8_peak=ReadPeakExpression("data/CA20230227_RW10265_tbx-8_L1.csv")

All_Exp=data.frame(nhr67_peak,ceh36_peak,mls2_peak,ceh13_peak,unc30_peak,nob1_peak,elt1_peak,ref2_peak,unc130_peak,ceh32_peak,tbx8_peak)

pdf("Expression/Between_Genes_ExpressionComparisons.pdf")
plot(All_Exp)
dev.off()
#Being extracted
#Needs to be fully edited
#ceh27_peak=ReadPeakExpression("data/CA20141028_ceh-27_JIM290_L4.csv"/)

tmpCellNames=read.csv("CellNames.csv")
CellNames=as.vector(tmpCellNames[,2])
names(CellNames)=tmpCellNames[,1]

write.csv(ceh36_peak,file="Expression/ceh36_peak.csv")
write.csv(All_Exp,file="Expression/Peak_Expression_All.csv")

source("AnalyzeDivTimes.R")

source("AnalyzePositions.R")

source("AnalyzeRotation.R")

source("PlotDefectSummaries.R")
    #FIXME - ADD ROTATION CORRECTED ANGLE ANALYSIS
source("PlotPositionDevs.R")




StartDir=getwd()
die1_cc_dev = AnalyzeDivTimes("die-1")
setwd(StartDir)
die1_dev=AnalyzePositions("die-1", Expression=ceh32_peak, CalculateNeighbors=TRUE)
#ceh32_dev=AnalyzePositions("ceh-32_mutant", Expression=ceh32_peak, CalculateNeighbors=FALSE)
setwd(StartDir)
AnalyzeRotation("die-1")
setwd(StartDir)
PlotDefectSummaries("die-1", ceh32_peak, minDivTime=50)
PlotDefectSummaries("die-1", ceh32_peak,sig=5,microns=5,expCutoff=500,minDivTime=70)
PlotDefectSummaries("die-1", ceh32_peak,sig=3,microns=5,expCutoff=2000,minDivTime=70)
setwd(StartDir)



StartDir=getwd()
tbx_cc_dev = AnalyzeDivTimes("tbx-8_9_double")
setwd(StartDir)
tbx_dev=AnalyzePositions("tbx-8_9_double", Expression=tbx8_peak, CalculateNeighbors=TRUE)
setwd(StartDir)
AnalyzeRotation("tbx-8_9_double")
setwd(StartDir)
PlotDefectSummaries("tbx-8_9_double", tbx8_peak, minDivTime=50)
PlotDefectSummaries("tbx-8_9_double", tbx8_peak,sig=5,microns=5,expCutoff=500,minDivTime=70)
PlotDefectSummaries("tbx-8_9_double", tbx8_peak,sig=3,microns=5,expCutoff=2000,minDivTime=70)
setwd(StartDir)



StartDir=getwd()
ceh32_cc_dev = AnalyzeDivTimes("ceh-32_mutant")
setwd(StartDir)
ceh32_dev=AnalyzePositions("ceh-32_mutant", Expression=ceh32_peak, CalculateNeighbors=TRUE)
#ceh32_dev=AnalyzePositions("ceh-32_mutant", Expression=ceh32_peak, CalculateNeighbors=FALSE)
setwd(StartDir)
AnalyzeRotation("ceh-32_mutant")
setwd(StartDir)
PlotDefectSummaries("ceh-32_mutant", ceh32_peak, minDivTime=50)
PlotDefectSummaries("ceh-32_mutant", ceh32_peak,sig=5,microns=5,expCutoff=500,minDivTime=70)
PlotDefectSummaries("ceh-32_mutant", ceh32_peak,sig=3,microns=5,expCutoff=2000,minDivTime=70)
setwd(StartDir)


ref2MeanPosDevs = PlotDeviationsList(Name="ceh-32_mutant",exp=ceh32_peak, t=200)
PlotExpVsDev("ceh-32_mutant",outfile="ceh-32_mutant_ExpVsDev.pdf",exp=ceh32_peak)
PlotExpVsDev("Richard_et_al_plus_comma_WT",outfile="WT_ExpVsDev_ceh-32_exp.pdf",exp=ceh32_peak)




StartDir=getwd()
ref2_cc_dev = AnalyzeDivTimes("ref2_mutant")
setwd(StartDir)
#ref2_dev=AnalyzePositions("ref2_mutant", Expression=ref2_peak, CalculateNeighbors=TRUE)
ref2_dev=AnalyzePositions("ref2_mutant", Expression=ref2_peak, CalculateNeighbors=FALSE)
setwd(StartDir)
AnalyzeRotation("ref2_mutant")
setwd(StartDir)
##PlotDefectSummaries("ref2_mutant", ref2_peak)
PlotDefectSummaries("ref2_mutant", ref2_peak, minDivTime=50)
PlotDefectSummaries("ref2_mutant", ref2_peak,sig=5,microns=5,expCutoff=100,minDivTime=70)

setwd(StartDir)
ref2MeanPosDevs = PlotDeviationsList(Name="ref2_mutant",exp=ref2_peak, t=200)
PlotExpVsDev("ref2_mutant",outfile="ref-2_mutant_ExpVsDev.pdf",exp=ref2_peak)
PlotExpVsDev("Richard_et_al_plus_comma_WT",outfile="WT_ExpVsDev_ref-2_exp.pdf",exp=ref2_peak)




StartDir=getwd()
unc130_cc_dev = AnalyzeDivTimes("unc-130mutant")
setwd(StartDir)
#unc130_dev=AnalyzePositions("unc-130mutant", Expression=unc130_peak, CalculateNeighbors=TRUE)
unc130_dev=AnalyzePositions("unc-130mutant", Expression=unc130_peak, CalculateNeighbors=FALSE)
setwd(StartDir)
AnalyzeRotation("unc-130mutant")
PlotDefectSummaries("unc-130mutant", ref2_peak)
setwd(StartDir)
PlotDefectSummaries("unc-130mutant", unc130_peak,sig=5,microns=5)
unc130MeanPosDevs = PlotDeviationsList(Name="unc-130mutant",exp=unc130_peak, t=250)



StartDir=getwd()
ceh51_cc_dev = AnalyzeDivTimes("20190715_ceh-51.txt")
setwd(StartDir)
ceh51_dev=AnalyzePositions("20190715_ceh-51.txt", Expression=ref2_peak, CalculateNeighbors=FALSE)
setwd(StartDir)
AnalyzeRotation("20190715_ceh-51.txt")
setwd(StartDir)
PlotDefectSummaries("20190715_ceh-51.txt", ref2_peak)
setwd(StartDir)




StartDir=getwd()
cec4_cc_dev = AnalyzeDivTimes("cec4_mutant.txt")
setwd(StartDir)
cec4_dev=AnalyzePositions("cec4_mutant.txt", Expression=elt1_peak, CalculateNeighbors=TRUE)
setwd(StartDir)
AnalyzeRotation("cec4_mutant.txt")
setwd(StartDir)



StartDir=getwd()
elt1_cc_dev = AnalyzeDivTimes("elt-1_mutant")
setwd(StartDir)
elt1_dev=AnalyzePositions("elt-1_mutant", Expression=elt1_peak)
#ceh13_dev=AnalyzePositions("elt-1_mutant", Expression=elt1_peak, CalculateNeighbors=TRUE)
setwd(StartDir)
AnalyzeRotation("elt-1_mutant")
setwd(StartDir)



StartDir=getwd()
AnalyzeDivTimes("ceh-13_mutant")
setwd(StartDir)
ceh13_dev=AnalyzePositions("ceh-13_mutant", Expression=ceh13_peak)
#ceh13_dev=AnalyzePositions("ceh-13_mutant", Expression=ceh13_peak, CalculateNeighbors=TRUE)
setwd(StartDir)
AnalyzeRotation("ceh-13_mutant")
setwd(StartDir)
PlotDefectSummaries("ceh-13_mutant", ceh13_peak,sig=5,microns=5)
ceh13MeanPosDevs = PlotDeviationsList(Name="ceh-13_mutant",exp=ceh13_peak, t=200)


StartDir=getwd()
AnalyzeDivTimes("nob-1_mutant")
setwd(StartDir)
#nob1_dev=AnalyzePositions("nob-1_mutant", Expression=nob1_peak)
nob1_dev=AnalyzePositions("nob-1_mutant", Expression=nob1_peak, CalculateNeighbors=TRUE)
setwd(StartDir)
AnalyzeRotation("nob-1_mutant")
setwd(StartDir)
PlotDefectSummaries("nob-1_mutant", nob1_peak,sig=5,microns=5)
nob1MeanPosDevs = PlotDeviationsList("nob-1_mutant",exp=nob1_peak, t=200)
plot(nob1_peak[nob1_dev$Cell],nob1_dev[,3])

PlotDefectSummaries("nob-1_mutant", php3_peak,sig=5,microns=5,expCutoff=100)
nob1MeanPosDevs = PlotDeviationsList("nob-1_mutant",exp=nob1_peak, t=200)
plot(nob1_peak[nob1_dev$Cell],nob1_dev[,3])




pdf("ceh13_nob1_Pos_comparisons.pdf")
xcor = round(cor(nob1MeanPosDevs$x,ceh13MeanPosDevs$x,use="pairwise.complete.obs"),3)
ycor = round(cor(nob1MeanPosDevs$y,ceh13MeanPosDevs$y,use="pairwise.complete.obs"),3)
zcor = round(cor(nob1MeanPosDevs$z,ceh13MeanPosDevs$z,use="pairwise.complete.obs"),3)
cor(unlist(nob1MeanPosDevs),unlist(ceh13MeanPosDevs),use="pairwise.complete.obs")
plot(nob1MeanPosDevs$x,ceh13MeanPosDevs$x,cex=0.1,xlab="nob-1 Xdev", ylab=paste("ceh-13 Xdev, cor:",xcor))
plot(nob1MeanPosDevs$y,ceh13MeanPosDevs$y,cex=0.1,col=2,xlab="nob-1 Ydev", ylab=paste("ceh-13 Ydev, cor:",ycor))
plot(nob1MeanPosDevs$z,ceh13MeanPosDevs$z,cex=0.1,col=3,xlab="nob-1 Zdev", ylab=paste("ceh-13 Zdev, cor:",zcor))
CellDots = rowSums(nob1MeanPosDevs*ceh13MeanPosDevs)
plot(sort(CellDots),xlab="rank",ylab="ceh-13 / nob-1 pos deviation dot product")

dev.off()
   




PlotExpVsDev("nob-1_mutant",outfile="nob-1_mutant_ExpVsDev.pdf",exp=nob1_peak)
PlotExpVsDev("ceh-13_mutant",outfile="ceh-13_mutant_ExpVsDev.pdf",exp=ceh13_peak)
PlotExpVsDev("Richard_et_al_plus_comma_WT",outfile="WT_ExpVsDev_nob-1_exp.pdf",exp=nob1_peak)
PlotExpVsDev("Richard_et_al_plus_comma_WT",outfile="WT_ExpVsDev_ceh-13_exp.pdf",exp=ceh13_peak)

PlotExpVsDev("ceh-36",outfile="ceh-36_mutant_ExpVsDev.pdf",exp=ceh36_peak)
ceh36MeanPosDevs = PlotDeviationsList("ceh-36",exp=ceh36_peak,t=200)



PlotDefectSummaries("Richard_et_al_plus_comma_WT", nob1_peak,sig=5,microns=5)
PlotDefectSummaries("Richard_et_al_plus_comma_WT", ceh13_peak,sig=5,microns=5)



PlotDefectSummaries("ceh-13_mutant", ceh13_peak)
setwd(StartDir)
PlotDefectSummaries("ceh-36", ceh36_peak)
setwd(StartDir)
PlotDefectSummaries("unc-30", unc30_peak)
setwd(StartDir)
PlotDefectSummaries("nob-1_mutant", nob1_peak)
setwd(StartDir)
PlotDefectSummaries("mls-2_mutants", mls2_peak)
setwd(StartDir)
PlotDefectSummaries("ceh-36_unc-30", ceh36_peak)
setwd(StartDir)
PlotDefectSummaries("nhr-67_mutant_edited", nhr67_peak)
setwd(StartDir)
PlotDefectSummaries("Richard_et_al_plus_comma_WT", nhr67_peak)
setwd(StartDir)
PlotDefectSummaries("mls2ceh36.txt", mls2_peak)
setwd(StartDir)



###load CC and pos dev tables for different mutants for clustering
##ceh-13, nob-1, ceh-36, unc-30, mls-2
###plot 3D vector maps of cell mispositioning?



## strains for embryo seq
# nob-1 - rescued with dup
# ceh-36 JIM
# unc-30
# ceh-36_unc-30 - rescued with extrachromosomal transgene
# ceh-13 - rescued with extrachromosomal transgene
# nhr-67 - rescued with extrachromosomal transgene

stop("Not an error - just reached the stop point!")




##CODE to load position devs for each mutant
#
##FIXMES - also get a single metric for 
##Instead report position devs as 3D vector. 
##Find way to distinguish delay from misdirection
##How to report expression for double mutants?
##Output summary of cell fates for progeny

WT_CCdev=AnalyzeDivTimes("Richard_et_al_plus_comma_WT")
setwd(StartDir)
ceh13_CCdev=AnalyzeDivTimes("ceh-13_mutant")
setwd(StartDir)
nob1_CCdev=AnalyzeDivTimes("nob-1_mutant") 
setwd(StartDir)
ceh36_CCdev=AnalyzeDivTimes("ceh-36") 
setwd(StartDir)
ceh36_unc30_CCdev=AnalyzeDivTimes("ceh-36_unc-30")
setwd(StartDir)
unc30_CCdev=AnalyzeDivTimes("unc-30")
setwd(StartDir)
mls2_CCdev=AnalyzeDivTimes("mls-2_mutants") 
setwd(StartDir)
mls2_ceh36_CCdev=AnalyzeDivTimes("mls2ceh36.txt")
setwd(StartDir)
nhr67_CCdev=AnalyzeDivTimes("nhr-67_mutant_edited")
setwd(StartDir)

CC_Combined=data.frame(WT_CCdev[Cells,],ceh13_CCdev[Cells,],nob1_CCdev[Cells,],ceh36_CCdev[Cells,],
                       ceh36_unc30_CCdev[Cells,],unc30_CCdev[Cells,],mls2_CCdev[Cells,],
                       mls2_ceh36_CCdev[Cells,],nhr67_CCdev[Cells,])
write.table(CC_Combined,file="CC_Combined.tsv",sep="\t",na="",quote=FALSE,col.names=NA)


AnalyzeRotation("Richard_et_al_plus_comma_WT")
setwd(StartDir)
AnalyzeRotation("ceh-13_mutant")
setwd(StartDir)
AnalyzeRotation("nob-1_mutant") 
setwd(StartDir)
AnalyzeRotation("ceh-36")
setwd(StartDir)
AnalyzeRotation("unc-30")
setwd(StartDir)
AnalyzeRotation("ceh-36_unc-30")
setwd(StartDir)
AnalyzeRotation("mls-2_mutants")
setwd(StartDir)
AnalyzeRotation("mls2ceh36.txt")
setwd(StartDir)
AnalyzeRotation("nhr-67_mutant_edited")
setwd(StartDir)




WT_dev=AnalyzePositions("Richard_et_al_plus_comma_WT")
setwd(StartDir)
ceh13_dev=AnalyzePositions("ceh-13_mutant", Expression=ceh13_peak)
setwd(StartDir)
ceh36_dev=AnalyzePositions("ceh-36",Expression=ceh36_peak)
setwd(StartDir)
ceh36_unc30_dev=AnalyzePositions("ceh-36_unc-30",Expression=unc30_peak)
setwd(StartDir)
unc30_dev=AnalyzePositions("unc-30",Expression=unc30_peak)
setwd(StartDir)
mls2_dev=AnalyzePositions("mls-2_mutants",Expression=mls2_peak)
setwd(StartDir)
mls2_ceh36_dev=AnalyzePositions("mls2ceh36.txt",Expression=mls2_peak)
setwd(StartDir)
nhr67_dev=AnalyzePositions("nhr-67_mutant_edited",Expression=mls2_peak)
setwd(StartDir)

nob1_dev=AnalyzePositions("nob-1_mutant",Expression=nob1_peak) ## needs NN
setwd(StartDir)

CellTimes=unique(c(rownames(WT_dev),rownames(ceh13_dev),rownames(nob1_dev),rownames(ceh36_dev),rownames(ceh36_unc30_dev),rownames(unc30_dev),rownames(mls2_dev),rownames(mls2_ceh36_dev),rownames(nhr67_dev)))

Position_Combined_Dev=data.frame(WT_dev[CellTimes,-1:-5],ceh13_dev[CellTimes,-1:-5],ceh36_dev[CellTimes,-1:-5],ceh36_unc30_dev[CellTimes,-1:-5],unc30_dev[CellTimes,-1:-5],
                             mls2_dev[CellTimes,-1:-5],mls2_ceh36_dev[CellTimes,-1:-5],nhr67_dev[CellTimes,-1:-5],nob1_dev[CellTimes,-1:-5])
write.table(Position_Combined_Dev,file="Pos_Dev_Combined.tsv",sep="\t",na="",quote=FALSE,col.names=NA)

CellLookups = WT_dev[CellTimes,1]


Position_Combined_CellMeans = aggregate(Position_Combined_Dev,by=list(CellLookups[,1]),mean,na.rm=T)
rownames(Position_Combined_CellMeans)=Position_Combined_CellMeans[,1]
write.table(Position_Combined_CellMeans,file="Pos_Dev_CellMeans.tsv",sep="\t",na="",quote=FALSE,col.names=NA)

Position_MeanDevs = data.frame(WT_dev[CellTimes,3],ceh13_dev[CellTimes,3],ceh36_dev[CellTimes,3],ceh36_unc30_dev[CellTimes,3],unc30_dev[CellTimes,3],
                             mls2_dev[CellTimes,3],mls2_ceh36_dev[CellTimes,3],nhr67_dev[CellTimes,3],nob1_dev[CellTimes,3])
write.table(Position_MeanDevs,file="Pos_Dev_Means.tsv",sep="\t",na="",quote=FALSE,col.names=NA)


Position_MeanDevs_CellMeans= aggregate(Position_MeanDevs,by=list(CellLookups[,1]),mean,na.rm=T)
rownames(Position_MeanDevs_CellMeans)=Position_MeanDevs_CellMeans[,1]
write.table(Position_MeanDevs_CellMeans,file="Pos_Dev_Means_CellMeans.tsv",sep="\t",na="",quote=FALSE,col.names=NA)

##Also zen-4? pop-1, sys-1, lag-1 RNAi etc? Stressed embryos?
#
#
#
###





PlotDefectSummaries("nhr-67_mutant_edited", nhr67_peak)


StartDir=getwd()
AnalyzeDivTimes("nob-1_mutant")
setwd(StartDir)
nob1_dev=AnalyzePositions("nob-1_mutant",Expression=nob1_peak)
#nob1_dev=AnalyzePositions("nob-1_mutant",Expression=nob1_peak,CalculateNeighbors=TRUE)
setwd(StartDir)
AnalyzeRotation("nob-1_mutant")
setwd(StartDir)





AnalyzeDivTimes("SM2373")
zen4_dev=AnalyzePositions("SM2373")
AnalyzeRotation("SM2373")



AnalyzeDivTimes("unc-30")
unc30_dev=AnalyzePositions("unc-30")
AnalyzeRotation("unc-30")

AnalyzeDivTimes("SM2373")
zen_dev=AnalyzePositions("SM2373")
AnalyzeRotation("SM2373")


#rescued_dev=AnalyzePositions("ceh-36_rescued")
#double_rescued_dev=AnalyzePositions("ceh-36_unc-30_rescued")

#ceh_unc_double_dev=AnalyzePositions("ceh-36_unc-30")

#RichardsCommaDev=AnalyzePositions("Richard_et_al_plus_comma_WT")

nob_dev=AnalyzePositions("nob-1")

ceh_dev=AnalyzePositions("ceh-36")

mls_dev=AnalyzePositions("mls-2_mutants")
mls_ceh_double_dev=AnalyzePositions("mls2ceh36.txt")
sys1_rnai_dev=AnalyzePositions("ALZ_sys-1RNAi.txt")

AnalyzeDivTimes("ALZ_sys-1RNAi.txt")
AnalyzeDivTimes("Richard_et_al_plus_comma_WT")
AnalyzeDivTimes("ceh-36_unc-30_rescued")
AnalyzeDivTimes("ceh-36")
AnalyzeDivTimes("nob-1")
AnalyzeDivTimes("mls-2_mutants")
AnalyzeDivTimes("mls2ceh36.txt")
AnalyzeDivTimes("ceh-36_rescued")

AnalyzeDivTimes("mutantsCombined")
AnalyzeDivTimes("ceh-36_unc-30")

AnalyzeRotation("ceh-36_unc-30")
AnalyzeRotation("ALZ_SYS-1RNAi.txt")
AnalyzeRotation("Richard_et_al_plus_comma_WT")
AnalyzeRotation("ceh-36_unc-30_rescued")
AnalyzeRotation("ceh-36")
AnalyzeRotation("nob-1")
AnalyzeRotation("mls-2_mutants")
AnalyzeRotation("mls2ceh36.txt")
AnalyzeRotation("ceh-36_rescued")




#AllMutantsDev=AnalyzePositions("mutantsCombined")

WTMeanDevs=rowMeans(RichardsCommaDev[,3:length(RichardsCommaDev)],na.rm=T)
WTdevSDs=apply(RichardsCommaDev[,3:length(RichardsCommaDev)],1,sd,na.rm=T)
WTMeanDevsperCell=tapply(WTMeanDevs,RichardsCommaDev[,1],mean,na.rm=T)
WTdevSDsperCell=tapply(WTdevSDs,RichardsCommaDev[,1],mean,na.rm=T)
#wtDevCellMeans=tapply(RichardsCommaDev[,3:length(RichardsCommaDev)],RichardsCommaDev[,1],mean,na.rm=T)
#wtDevTimeMeans=tapply(RichardsCommaDev[,3:length(RichardsCommaDev)],RichardsCommaDev[,1],mean,na.rm=T)




pdf("PositionDevComparisons.pdf")

test=data.frame(RichardsCommaDev[,1],RichardsCommaDev[,2],WTMeanDevs)
colnames(test)=c("Cell","Time","MeanDev")
test$Cell<-with(test,reorder(Cell,MeanDev,mean))
boxplot(test$MeanDev~test$Cell, main="WT dev vs cell")
boxplot(test$MeanDev~test$Time, main="WT dev vs time")



for(i in c("RichardsCommaDev","AllMutantsDev")){
    #),"ceh_unc_double_dev","ceh_dev")){

     thesedata=get(i)



    fudge= median(mutantMeanDevsperCell)-median(WTMeanDevsperCell)
    thesedata=thesedata+fudge



    mutantMeanDevs=rowMeans(thesedata[,3:length(thesedata)],na.rm=T)
    mutantdevSDs=apply(thesedata[,3:length(thesedata)],1,sd,na.rm=T)
    mutantMeanDevsperCell=tapply(mutantMeanDevs,thesedata[,1],mean,na.rm=T)
    mutantdevSDsperCell=tapply(mutantdevSDs,thesedata[,1],mean,na.rm=T)

#    mutantDevCellMeans=tapply(thesedata[,3:length(thesedata)],RichardsCommaDev[,1],mean,na.rm=T)
     #    mutantDevTimeMeans=tapply(thesedata[,3:length(thesedata)],RichardsCommaDev[,1],mean,na.rm=T)


    model=lm(mutantMeanDevsperCell~WTMeanDevsperCell[names(mutantMeanDevsperCell)])
    plot(WTMeanDevsperCell[names(mutantMeanDevsperCell)],mutantMeanDevsperCell, xlab="WT Devs",ylab = paste(i,"Devs"),
         main=paste("rsq:",summary(model)$adj.r.squared,
                    "intercept",summary(model)$coefficients[1,1],
                    "slope",summary(model)$coefficients[2,1])
         )
    abline(model)

    boxplot(mutantMeanDevs~thesedata[,1],ylim=c(0,10),main=i)


    

    plot(mutantMeanDevsperCell,mutantdevSDsperCell,main=i,xlim=c(0,10),ylim=c(0,10))
    points(WTMeanDevsperCell,WTdevSDsperCell,col=2)
    

    test=data.frame(thesedata[,1],thesedata[,2],mutantMeanDevs)
    colnames(test)=c("Cell","Time","MeanDev")
    test$Cell<-with(test,reorder(Cell,MeanDev,mean))
    boxplot(test$MeanDev~test$Cell, main=paste(i, " dev vs cell"))
    boxplot(test$MeanDev~test$Time, main=paste(i, " dev vs time"))

}

dev.off()



#This requires first running "CompareDivTime.pl" on lists with each name





# OLD STUFF
#                 #All times too crazy - try subset

#         for(j in c(50,100,150,200,250,300,350)){
#             Pos=data.frame(mutantX[MutantPositions[2]==j,i],mutantY[MutantPositions[2]==j,i],mutantZ[MutantPositions[2]==j,i],row.names=MutantPositions[MutantPositions[2]==j,1])
#             colnames(Pos)=c("X","Y","Z")
#             if(length(Pos[,1])<5){
#                 next()
#                 }
#             #find best fit time match
#             BestTime = j
#             WTMeanPos=data.frame(wtXMeans[WTPositions[2]==BestTime],wtYMeans[WTPositions[2]==BestTime],wtZMeans[WTPositions[2]==BestTime],row.names=WTPositions[WTPositions[2]==BestTime,1])
#             colnames(WTMeanPos)=c("X","Y","Z")
#             WTMeanPos=WTMeanPos[rownames(Pos),]
#             WTMutantDevs=WTMeanPos-Pos
#             WTMutantDists=sqrt(rowSums(WTMutantDevs^2))
#             MeanDist=mean(WTMutantDists,na.rm=T)

#             DistOpt=c(j,MeanDist)
#             for(k in (j-20):(j+20)){
#                 theseWTPos=data.frame(wtXMeans[WTPositions[2]==k],wtYMeans[WTPositions[2]==k],wtZMeans[WTPositions[2]==k],row.names=WTPositions[WTPositions[2]==k,1])
#                 colnames(theseWTPos)=c("X","Y","Z")
#                 theseWTPos=theseWTPos[rownames(Pos),]
#                 theseDevs=theseWTPos-Pos
#                 theseDists=sqrt(rowSums(theseDevs^2))
#                 thisMeanDist=mean(theseDists,na.rm=T)
#                 DistOpt=rbind(DistOpt,c(k,thisMeanDist))
#                 # if(!is.na(thisMeanDist)&&!is.na(MeanDist)){
#                 #         if(thisMeanDist<MeanDist){
#                 #             BestTime=k
#                 #             WTMeanPos=theseWTPos
#                 #             WTMutantDevs=theseDevs
#                 #             WTMutantDists=theseDists
#                 #             MeanDist=thisMeanDist
#                 #             }
#                 #         }
#                 }
#             #currently not updating time - more of a reality check
#             try(plot(DistOpt,xlab="WT time point", ylab="mean distance", main=paste(i,"Best Time Match for t=",j)),silent=T)
            
#  #Find optimum rotation around x axis
#             BestTheta=0
#             thetaOpt=NULL
#             for(k in (-20:20)/60){

#                 #y' = y cos f - z sin f
#                 #z' = z cos f + y sin f 
#                 thesePos=data.frame(Pos[,"X"],Pos[,"Y"]*cos(k)-Pos[,"Z"]*sin(k),Pos[,"Z"]*cos(k)+Pos[,"Y"]*sin(k))
#                 colnames(thesePos)=c("X","Y","Z")
#                 rownames(thesePos)=rownames(Pos)

#                 theseWTMutantDevs=WTMeanPos-thesePos
#                 theseWTMutantDists=sqrt(rowSums(theseWTMutantDevs^2))
#                 thisMeanDist=mean(theseWTMutantDists,na.rm=T)
#                 thetaOpt=rbind(thetaOpt,c(k,thisMeanDist))
#                 if(!is.na(thisMeanDist)&&!is.na(MeanDist)){
#                     if(thisMeanDist<MeanDist){
#                         BestTheta=k
#                         Pos=thesePos
#                         WTMutantDevs=theseWTMutantDevs
#                         WTMutantDists=theseWTMutantDists
#                         MeanDist=thisMeanDist
#                         }
                            
#                     }
#                 }

#             #TODO: Should also optimize x/y/z scaling after rotation! Check whether scaling and rotation can be optimized at a single time point.
                

#             try(plot(thetaOpt,xlab="theta", ylab="mean distance", main=paste(i,"Best theta match for yz rotation")),silent=T)




#             #plot xy space deviations
#             plot(WTMeanPos[,"X"],WTMeanPos[,"Y"],main=paste(i,"XY Deviations at t=",j,"wtT=",BestTime,"Theta=",BestTheta),xlim=c(-30,30),ylim=c(-30,30))
#             points(Pos[,"X"],Pos[,"Y"],col=2)
#             segments(WTMeanPos[,"X"],WTMeanPos[,"Y"],Pos[,"X"],Pos[,"Y"])

#             #plot xz space deviations
#             plot(WTMeanPos[,"X"],WTMeanPos[,"Z"],main=paste(i,"XZ Deviations at t=",j,"wtT=",BestTime,"Theta=",BestTheta),xlim=c(-30,30),ylim=c(-30,30))
#             points(Pos[,"X"],Pos[,"Z"],col=2)
#             segments(WTMeanPos[,"X"],WTMeanPos[,"Z"],Pos[,"X"],Pos[,"Z"])
            
# #           plot3d(WTMeanPos[,"X"],WTMeanPos[,"Y"],WTMeanPos[,"Z"],size=10,xlim=c(-30,30),ylim=c(-30,30),zlim=c(-30,30),main=i)
# #           points3d(Pos[,"X"],Pos[,"Y"],Pos[,"Z"],size=10,col=2)
# #            segments3d(interleave(as.matrix(data.frame(WTMeanPos[,"X"],WTMeanPos[,"Y"],WTMeanPos[,"Z"])),as.matrix(data.frame(Pos[,"X"],Pos[,"Y"],Pos[,"Z"]))))
#             #These didn't work - segments requires a single matrix of xyz coordinats where the start and end of each line are in alternating rows hence the interleave call above
#             #segments3d(data.frame(WTMeanPos[,"X"],WTMeanPos[,"Y"],WTMeanPos[,"Z"]),data.frame(Pos[,"X"],Pos[,"Y"],Pos[,"Z"]))
#             #abclines3d(WTMeanPos[,"X"],WTMeanPos[,"Y"],WTMeanPos[,"Z"],
#             #           Pos[,"X"]-WTMeanPos[,"X"],Pos[,"Y"]-WTMeanPos[,"Y"],Pos[,"Z"]-WTMeanPos[,"Z"])



#             dev.off()
#             #END LOOP THROUGH TIMES
#             }

