
RULES=read.table("Over500WT/Over500WT.angles_unaligned.txt.average")

library("scatterplot3d")

x=c(c(RULES[,2]),c(-RULES[,2]))
y=c(c(RULES[,3]),c(-RULES[,3]))
z=c(c(RULES[,4]),c(-RULES[,4]))

pdf('angles.pdf')

scatterplot3d(x,y,z,highlight.3d=T,pch=16,cex.symbols=0.8,xlab="A-P",ylab="L-R",zlab="D-V",grid=F)
scatterplot3d(y,x,z,highlight.3d=T,pch=16,cex.symbols=0.8,xlab="L-R",ylab="A-P",zlab="D-V",grid=F)


#Compare cells for higher/lower variability in wt

ANGLES=read.table("Over500WT/Over500WT.angles_unaligned.txt",header=T,stringsAsFactors=F)


movies=unique(ANGLES[,1])
cells=unique(ANGLES[,2])


boxplot(abs(dot)~parent,data=ANGLES[substr(start=1,stop=5,x=ANGLES[,2])=="ABala",],cex.axis=0.3)
boxplot(abs(dot)~parent,data=ANGLES[substr(start=1,stop=5,x=ANGLES[,2])=="ABalp",],cex.axis=0.3)
boxplot(abs(dot)~parent,data=ANGLES[substr(start=1,stop=5,x=ANGLES[,2])=="ABara",],cex.axis=0.3)
boxplot(abs(dot)~parent,data=ANGLES[substr(start=1,stop=5,x=ANGLES[,2])=="ABarp",],cex.axis=0.3)
boxplot(abs(dot)~parent,data=ANGLES[substr(start=1,stop=5,x=ANGLES[,2])=="ABpla",],cex.axis=0.3)
boxplot(abs(dot)~parent,data=ANGLES[substr(start=1,stop=5,x=ANGLES[,2])=="ABplp",],cex.axis=0.3)
boxplot(abs(dot)~parent,data=ANGLES[substr(start=1,stop=5,x=ANGLES[,2])=="ABpra",],cex.axis=0.3)
boxplot(abs(dot)~parent,data=ANGLES[substr(start=1,stop=5,x=ANGLES[,2])=="ABprp",],cex.axis=0.3)
boxplot(abs(dot)~parent,data=ANGLES[substr(start=1,stop=3,x=ANGLES[,2])=="MSa",],cex.axis=0.3)
boxplot(abs(dot)~parent,data=ANGLES[substr(start=1,stop=3,x=ANGLES[,2])=="MSp",],cex.axis=0.3)
boxplot(abs(dot)~parent,data=ANGLES[substr(start=1,stop=1,x=ANGLES[,2])=="E",],cex.axis=0.3)
boxplot(abs(dot)~parent,data=ANGLES[substr(start=1,stop=1,x=ANGLES[,2])=="C",],cex.axis=0.3)
boxplot(abs(dot)~parent,data=ANGLES[substr(start=1,stop=1,x=ANGLES[,2])=="D",],cex.axis=0.3)
boxplot(abs(dot)~parent,data=ANGLES[substr(start=1,stop=1,x=ANGLES[,2])=="P",],cex.axis=0.3)
boxplot(abs(dot)~series,data=ANGLES,cex.axis=0.1)


means = c()
sds = c()
lengths = c()
DivTimesForAngles=c()
CC_CVsForAngles=c()

for(i in cells){
    dots=ANGLES[ANGLES[,2]==i,17]
    #    print(c(i,mean(dots),sd(dots),length(dots)))
    means=c(means,mean(dots))
    sds=c(sds,sd(dots))
    lengths=c(lengths,length(dots))
    DivTimesForAngles=c(DivTimesForAngles,DivTimeAverages[i])
    if(DivTimeAverages[i]>80){
        CC_CVsForAngles=c(CC_CVsForAngles,CV[i])
        }else{
        CC_CVsForAngles=c(CC_CVsForAngles,NA)
    
        }
}

means22=means
lengths22=lengths
stderrs=sds/lengths
sds22=sds
DOTs=data.frame(cell=cells,mean=means,sd=sds,length=lengths,se=stderrs,divtime=DivTimesForAngles)
write.table(DOTs[order(-means),],file="CellDotScores.tsv",sep="\t")


plot(CC_CVsForAngles,acos(means22)/pi*180)
dev.off()
    



#significance of cells with increased variability
#aov.anglesSeries=aov(dot~series,data=ANGLES)
#aov.anglesCells=aov(abs(dot)~parent,data=ANGLES)


#aov.angles=aov(dot~series*parent,data=ANGLES)

# compare variability at high temp
ANGLES_ALL=read.table("Over500All/Over500All.angles_unaligned.txt",header=T,stringsAsFactors=F)



for(i in c("_25C_","_26C_","_30C_")){
    THESEANGLES=ANGLES_ALL[(grep(i,ANGLES_ALL[,1],)),]

    if(i=="_25C_"){
        THESEANGLES=rbind(ANGLES_ALL[(grep("_24C_",ANGLES_ALL[,1],)),],ANGLES_ALL[(grep("_25C_",ANGLES_ALL[,1],)),])
        }
    if(i=="_26C_"){
        THESEANGLES=rbind(THESEANGLES,ANGLES_ALL[(grep("20100",ANGLES_ALL[,1],)),])
        }
    pdf(paste("ANGLES",i,".pdf",sep=""))
    
    boxplot(abs(dot)~parent,data=THESEANGLES[substr(start=1,stop=5,x=THESEANGLES[,2])=="ABala",],cex.axis=0.3)
    boxplot(abs(dot)~parent,data=THESEANGLES[substr(start=1,stop=5,x=THESEANGLES[,2])=="ABalp",],cex.axis=0.3)
    boxplot(abs(dot)~parent,data=THESEANGLES[substr(start=1,stop=5,x=THESEANGLES[,2])=="ABara",],cex.axis=0.3)
    boxplot(abs(dot)~parent,data=THESEANGLES[substr(start=1,stop=5,x=THESEANGLES[,2])=="ABarp",],cex.axis=0.3)
    boxplot(abs(dot)~parent,data=THESEANGLES[substr(start=1,stop=5,x=THESEANGLES[,2])=="ABpla",],cex.axis=0.3)
    boxplot(abs(dot)~parent,data=THESEANGLES[substr(start=1,stop=5,x=THESEANGLES[,2])=="ABplp",],cex.axis=0.3)
    boxplot(abs(dot)~parent,data=THESEANGLES[substr(start=1,stop=5,x=THESEANGLES[,2])=="ABpra",],cex.axis=0.3)
    boxplot(abs(dot)~parent,data=THESEANGLES[substr(start=1,stop=5,x=THESEANGLES[,2])=="ABprp",],cex.axis=0.3)
    boxplot(abs(dot)~parent,data=THESEANGLES[substr(start=1,stop=3,x=THESEANGLES[,2])=="MSa",],cex.axis=0.3)
    boxplot(abs(dot)~parent,data=THESEANGLES[substr(start=1,stop=3,x=THESEANGLES[,2])=="MSp",],cex.axis=0.3)
    boxplot(abs(dot)~parent,data=THESEANGLES[substr(start=1,stop=1,x=THESEANGLES[,2])=="E",],cex.axis=0.3)
    boxplot(abs(dot)~parent,data=THESEANGLES[substr(start=1,stop=1,x=THESEANGLES[,2])=="C",],cex.axis=0.3)
    boxplot(abs(dot)~parent,data=THESEANGLES[substr(start=1,stop=1,x=THESEANGLES[,2])=="D",],cex.axis=0.3)
    boxplot(abs(dot)~parent,data=THESEANGLES[substr(start=1,stop=1,x=THESEANGLES[,2])=="P",],cex.axis=0.3)
    boxplot(abs(dot)~series,data=THESEANGLES,cex.axis=0.05)
    dev.off()


    means = c()
    sds = c()
    lengths = c()
    DivTimesForAngles=c()
    for(j in cells){
        dots=THESEANGLES[THESEANGLES[,2]==j,17]
        #    print(c(i,mean(dots),sd(dots),length(dots)))
        means=c(means,mean(dots))
        sds=c(sds,sd(dots))
        lengths=c(lengths,length(dots))
        DivTimesForAngles=c(DivTimesForAngles,DivTimeAverages[j])
        }
    stderrs=sds/lengths
    assign(paste("means", i,sep=""),means)
    assign(paste("sds", i,sep=""),sds)
    assign(paste("lengths", i,sep=""),lengths)

    
    THESEDOTs=data.frame(cell=cells,mean=means,sd=sds,length=lengths,se=stderrs,divtime=DivTimesForAngles)
    write.table(THESEDOTs[order(-means),],file=paste(i,"CellDotScores.tsv",sep=""),sep="\t")
    
    }
pdf("summary_angles_by_temp.pdf")
dot22=ANGLES[,17]
dot25=rbind(ANGLES_ALL[(grep("_24C_",ANGLES_ALL[,1],)),],ANGLES_ALL[(grep("_25C_",ANGLES_ALL[,1],)),])
dot26=rbind(ANGLES_ALL[(grep("26C",ANGLES_ALL[,1])),],ANGLES_ALL[(grep("20100",ANGLES_ALL[,1],)),])
dot30=ANGLES_ALL[(grep("30C",ANGLES_ALL[,1],)),]


boxplot(abs(dot22),abs(dot25[,17]),abs(dot26[,17]),abs(dot30[,17]))
boxplot(acos(abs(dot22))/pi*180,acos(abs(dot25[,17]))/pi*180,acos(abs(dot26[,17]))/pi*180,acos(abs(dot30[,17]))/pi*180)

print("Hi2.2!")

#angleMeans22=aggregate(dot~parent,data=ANGLES,FUN=mean)
#angleMeans25=aggregate(dot~parent,data=dot25,FUN=mean)

#angleMeans26=aggregate(dot~parent,data=dot26,FUN=mean)
#angleMeans30=aggregate(dot~parent,data=dot30,FUN=mean)

angleMeans22=aggregate(ANGLES$dot,list(ANGLES$parent),FUN=mean)
angleMeans25=aggregate(dot25$dot,list(dot25$parent),FUN=mean)
angleMeans26=aggregate(dot26$dot,list(dot26$parent),FUN=mean)
angleMeans30=aggregate(dot30$dot,list(dot30$parent),FUN=mean)



print("Hi2.3!")

combined=merge(angleMeans22,angleMeans25,by="Group.1")

combined=merge(combined,angleMeans26,by="Group.1")

combined=merge(combined,angleMeans30,by="Group.1")
#boxplot(acos(abs(combined[2:5]))/pi*180)
print("Hi3!")
boxplot(acos(abs(means22))/pi*180,acos(abs(means_25C_))/pi*180,acos(abs(means_26C_))/pi*180,acos(abs(means_30C_))/pi*180)

t.test(acos(abs(means22))/pi*180,acos(abs(means_25C_))/pi*180)
t.test(acos(abs(means22))/pi*180,acos(abs(means_26C_))/pi*180)
t.test(acos(abs(means22))/pi*180,acos(abs(means_30C_))/pi*180)



#these compare variability of each division between wildtype and stress to test whether the high-temp variable divisions are the same as the low-temp variable divisions
plot(acos(means22)/pi*180,acos(means_25C_)/pi*180)
x=acos(means22)/pi*180
y=acos(means_25C_)/pi*180
model=lm(x~y)
abline(model,col=2)
text(50,25,round({summary(model)}$r.squared,3),adj=c(0,0))
text(50,20,round(coef(model)[1],2),adj=c(0,0))
text(50,15,round(coef(model)[2],3),adj=c(0,0))
text(50,10,{summary(model)}$df[2],adj=c(0,0))

text(47,25,"Rsquared:",adj=c(1,0))
text(47,20,"intercept:",adj=c(1,0))
text(47,15,"slope:",adj=c(1,0))
text(47,10,"deg freedom:",adj=c(1,0))

plot(acos(means22)/pi*180,acos(means_26C_)/pi*180)
x=acos(means22)/pi*180
y=acos(means_26C_)/pi*180
model=lm(x~y)
abline(model,col=2)
text(50,25,round({summary(model)}$r.squared,3),adj=c(0,0))
text(50,20,round(coef(model)[1],2),adj=c(0,0))
text(50,15,round(coef(model)[2],3),adj=c(0,0))
text(50,10,{summary(model)}$df[2],adj=c(0,0))

text(47,25,"Rsquared:",adj=c(1,0))
text(47,20,"intercept:",adj=c(1,0))
text(47,15,"slope:",adj=c(1,0))
text(47,10,"deg freedom:",adj=c(1,0))

print("Hi4!")

plot(acos(means22)/pi*180,acos(means_30C_)/pi*180)
x=acos(means22)/pi*180
y=acos(means_30C_)/pi*180
model=lm(x~y)
abline(model,col=2)
text(50,25,round({summary(model)}$r.squared,3),adj=c(0,0))
text(50,20,round(coef(model)[1],2),adj=c(0,0))
text(50,15,round(coef(model)[2],3),adj=c(0,0))
text(50,10,{summary(model)}$df[2],adj=c(0,0))

text(47,25,"Rsquared:",adj=c(1,0))
text(47,20,"intercept:",adj=c(1,0))
text(47,15,"slope:",adj=c(1,0))
text(47,10,"deg freedom:",adj=c(1,0))


dev.off()
