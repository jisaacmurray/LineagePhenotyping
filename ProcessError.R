
ERRORS=read.table("ERRORS.txt",header=T);

CELLS=ERRORS[,"tp"]+ERRORS[,"fn"]

FP_50=sum(ERRORS[CELLS<=50,"fp"])/sum(CELLS[CELLS<50])
FP_101=sum(ERRORS[CELLS>50&CELLS<=101,"fp"])/sum(CELLS[CELLS>50&CELLS<=101])
FP_180=sum(ERRORS[CELLS>101&CELLS<=180,"fp"])/sum(CELLS[CELLS>101&CELLS<=180])
FP_194=sum(ERRORS[CELLS>180&CELLS<=194,"fp"])/sum(CELLS[CELLS>180&CELLS<=194])
FP_350=sum(ERRORS[CELLS>194&CELLS<=350,"fp"])/sum(CELLS[CELLS>194&CELLS<=350])
FP_500=sum(ERRORS[CELLS>350&CELLS<=500,"fp"])/sum(CELLS[CELLS>350&CELLS<=500])
FP_600=sum(ERRORS[CELLS>350,"fp"])/sum(CELLS[CELLS>350])

FP=c(FP_50,FP_101,FP_180,FP_194,FP_350,FP_500,FP_600)

FN_50=sum(ERRORS[CELLS<=50,"fn"])/sum(CELLS[CELLS<50])
FN_101=sum(ERRORS[CELLS>50&CELLS<=101,"fn"])/sum(CELLS[CELLS>50&CELLS<=101])
FN_180=sum(ERRORS[CELLS>101&CELLS<=180,"fn"])/sum(CELLS[CELLS>101&CELLS<=180])
FN_194=sum(ERRORS[CELLS>180&CELLS<=194,"fn"])/sum(CELLS[CELLS>180&CELLS<=194])
FN_350=sum(ERRORS[CELLS>194&CELLS<=350,"fn"])/sum(CELLS[CELLS>194&CELLS<=350])
FN_500=sum(ERRORS[CELLS>350&CELLS<=500,"fn"])/sum(CELLS[CELLS>350&CELLS<=500])
FN_600=sum(ERRORS[CELLS>350,"fn"])/sum(CELLS[CELLS>350])
FN=c(FN_50,FN_101,FN_180,FN_194,FN_350,FN_500,FN_600)

M_50=sum(ERRORS[CELLS<=50,"move"])/sum(CELLS[CELLS<50])
M_101=sum(ERRORS[CELLS>50&CELLS<=101,"move"])/sum(CELLS[CELLS>50&CELLS<=101])
M_180=sum(ERRORS[CELLS>101&CELLS<=180,"move"])/sum(CELLS[CELLS>101&CELLS<=180])
M_194=sum(ERRORS[CELLS>180&CELLS<=194,"move"])/sum(CELLS[CELLS>180&CELLS<=194])
M_350=sum(ERRORS[CELLS>194&CELLS<=350,"move"])/sum(CELLS[CELLS>194&CELLS<=350])
M_500=sum(ERRORS[CELLS>350&CELLS<=500,"move"])/sum(CELLS[CELLS>350&CELLS<=500])
M_600=sum(ERRORS[CELLS>350,"move"])/sum(CELLS[CELLS>350])
M=c(M_50,M_101,M_180,M_194,M_350,M_500,M_600)



labels=c("1-50","51-101","102-180","181-194","195-350","351-500",">500")

pdf('ERROR.pdf')
barplot(FP,names.arg=labels,xlab="Cell stage", ylab="False Positive Rate",cex.names=0.8,ylim=c(0,0.1))

barplot(FN,names.arg=labels,xlab="Cell stage", ylab="False Negative Rate",cex.names=0.8,ylim=c(0,0.1))
barplot(M,names.arg=labels,xlab="Cell stage", ylab="Mispositioned",cex.names=0.8,ylim=c(0,0.1))

FNs=data.frame(stage=labels,mean=FN)
FPs=data.frame(stage=labels,mean=FP)
Ms=data.frame(stage=labels,mean=M)
for(seriesname in unique(ERRORS[,1])){

    ERRORS=read.table("ERRORS.txt",header=T);
    ERRORS=ERRORS[ERRORS[,1]==seriesname,]


    CELLS=ERRORS[,"tp"]+ERRORS[,"fn"]
    
    FP_50=sum(ERRORS[CELLS<=50,"fp"])/sum(CELLS[CELLS<50])
    FP_101=sum(ERRORS[CELLS>50&CELLS<=101,"fp"])/sum(CELLS[CELLS>50&CELLS<=101])
    FP_180=sum(ERRORS[CELLS>101&CELLS<=180,"fp"])/sum(CELLS[CELLS>101&CELLS<=180])
    FP_194=sum(ERRORS[CELLS>180&CELLS<=194,"fp"])/sum(CELLS[CELLS>180&CELLS<=194])
    FP_350=sum(ERRORS[CELLS>194&CELLS<=350,"fp"])/sum(CELLS[CELLS>194&CELLS<=350])
    FP_500=sum(ERRORS[CELLS>350&CELLS<=500,"fp"])/sum(CELLS[CELLS>350&CELLS<=500])
    FP_600=sum(ERRORS[CELLS>350,"fp"])/sum(CELLS[CELLS>350])
    
    FP=c(FP_50,FP_101,FP_180,FP_194,FP_350,FP_500,FP_600)
    
    FN_50=sum(ERRORS[CELLS<=50,"fn"])/sum(CELLS[CELLS<50])
    FN_101=sum(ERRORS[CELLS>50&CELLS<=101,"fn"])/sum(CELLS[CELLS>50&CELLS<=101])
    FN_180=sum(ERRORS[CELLS>101&CELLS<=180,"fn"])/sum(CELLS[CELLS>101&CELLS<=180])
    FN_194=sum(ERRORS[CELLS>180&CELLS<=194,"fn"])/sum(CELLS[CELLS>180&CELLS<=194])
    FN_350=sum(ERRORS[CELLS>194&CELLS<=350,"fn"])/sum(CELLS[CELLS>194&CELLS<=350])
    FN_500=sum(ERRORS[CELLS>350&CELLS<=500,"fn"])/sum(CELLS[CELLS>350&CELLS<=500])
    FN_600=sum(ERRORS[CELLS>350,"fn"])/sum(CELLS[CELLS>350])
    FN=c(FN_50,FN_101,FN_180,FN_194,FN_350,FN_500,FN_600)
    
    M_50=sum(ERRORS[CELLS<=50,"move"])/sum(CELLS[CELLS<50])
    M_101=sum(ERRORS[CELLS>50&CELLS<=101,"move"])/sum(CELLS[CELLS>50&CELLS<=101])
    M_180=sum(ERRORS[CELLS>101&CELLS<=180,"move"])/sum(CELLS[CELLS>101&CELLS<=180])
    M_194=sum(ERRORS[CELLS>180&CELLS<=194,"move"])/sum(CELLS[CELLS>180&CELLS<=194])
    M_350=sum(ERRORS[CELLS>194&CELLS<=350,"move"])/sum(CELLS[CELLS>194&CELLS<=350])
    M_500=sum(ERRORS[CELLS>350&CELLS<=500,"move"])/sum(CELLS[CELLS>350&CELLS<=500])
    M_600=sum(ERRORS[CELLS>350,"move"])/sum(CELLS[CELLS>350])
    M=c(M_50,M_101,M_180,M_194,M_350,M_500,M_600)
    

    FNs[seriesname] = FN
    FPs[seriesname] = FP
    Ms[seriesname]=M
    

    barplot(FP,names.arg=labels,xlab="Cell stage", ylab=paste("False Positive Rate -",seriesname),cex.names=0.8,ylim=c(0,0.1))
    barplot(FN,names.arg=labels,xlab="Cell stage", ylab=paste("False Negative Rate -",seriesname),cex.names=0.8,ylim=c(0,0.1))
    barplot(M,names.arg=labels,xlab="Cell stage", ylab=paste("Mispositioned -",seriesname),cex.names=0.8,ylim=c(0,0.1))


    }


plot(FPs[,2],col=2,pch=3,cex=2,ylim=c(0,0.1))
for(i in 3:8){
    points(FPs[,i])
}

plot(FNs[,2],col=2,pch=3,cex=2,ylim=c(0,0.1))
for(i in 3:8){
    points(FNs[,i])
}

plot(Ms[,2],col=2,pch=3,cex=2,ylim=c(0,0.1))
for(i in 3:8){
    points(Ms[,i])
}


dev.off()



#ERRORS_1_50=ERRORS[,ERRORS
#ERRORS_51_101
#ERRORS_102_180
#ERRORS_181_194
#ERRORS_194_350
#ERRORS_351plus
