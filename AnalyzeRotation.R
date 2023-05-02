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


AnalyzeRotation <- function(Name, peak=ceh36_peak){
    message("Analyzing rotated division orientations",Name)
    
    WTRotatedX=read.csv("Richard_et_al_plus_comma_WT/Richard_et_al_plus_comma_WT_rotatedX.csv",row.names=1)
    WTRotatedY=read.csv("Richard_et_al_plus_comma_WT/Richard_et_al_plus_comma_WT_rotatedY.csv",row.names=1)
    WTRotatedZ=read.csv("Richard_et_al_plus_comma_WT/Richard_et_al_plus_comma_WT_rotatedZ.csv",row.names=1)
    
    startDir=getwd()
    setwd(Name)

    Times=data.frame(do.call('rbind', strsplit(as.character(rownames(WTRotatedX)),':',fixed=TRUE)))
    rownames(Times)=rownames(WTRotatedX)
    minTimes=tapply(as.numeric(as.character(Times[,2])),Times[,1],min)
    maxTimes=tapply(as.numeric(as.character(Times[,2])),Times[,1],max)


    WTXMeans=rowMeans(WTRotatedX,na.rm=T)
    WTYMeans=rowMeans(WTRotatedY,na.rm=T)
    WTZMeans=rowMeans(WTRotatedZ,na.rm=T)


    

    mutantRotatedX=read.csv(paste(Name,"_rotatedX.csv",sep=""),row.names=1)
    mutantRotatedY=read.csv(paste(Name,"_rotatedY.csv",sep=""),row.names=1)
    mutantRotatedZ=read.csv(paste(Name,"_rotatedZ.csv",sep=""),row.names=1)

    mutantXMeans=rowMeans(mutantRotatedX,na.rm=T)
    mutantYMeans=rowMeans(mutantRotatedY,na.rm=T)
    mutantZMeans=rowMeans(mutantRotatedZ,na.rm=T)

    #replace eventually with an sapply)
#    Mutant_Dot=data.frame(col.names=colnames(mutantRotatedX))
    Mutant_Dot=NULL
#    WT_Dot=data.frame(col.names=colnames(WTRotatedX))
    WT_Dot=NULL
    RowNames=NULL
    for(thisCell in Cells){
        thisSister=GetSister(thisCell)
        thisParent=GetParent(thisCell)

        Landmarks = c("AB", "E", "MS", "C" , "D")
        if(is.element(thisCell, Landmarks)){
            print(paste(thisCell,thisSister,thisParent))
            }
        birthTime=minTimes[thisCell]
        divTime=maxTimes[thisParent]
        if(!is.na(thisSister) && !is.element(thisParent,RowNames[,1])&& thisCell<thisSister && !is.na(WTXMeans[paste(thisCell,birthTime,sep=":")])  && !is.na(mutantXMeans[paste(thisCell,birthTime,sep=":")]) ){ 
            #only use a, l and d daughters to avoid redundancy
            
            WT_Vector=c(WTXMeans[paste(thisSister,birthTime,sep=":")]-WTXMeans[paste(thisCell,birthTime,sep=":")],
                        WTYMeans[paste(thisSister,birthTime,sep=":")]-WTYMeans[paste(thisCell,birthTime,sep=":")],
                        WTZMeans[paste(thisSister,birthTime,sep=":")]-WTZMeans[paste(thisCell,birthTime,sep=":")])
            WTUnitVector=WT_Vector/sqrt(sum(WT_Vector^2))
            
            #            rbind(Mutant_Dot,sapply(colnames(mutantRotatedX),function(X){
            #                    thisVector=c(mutantRotatedX[paste(thisSister,birthTime,sep=":"),X]-mutantRotatedX[paste(thisCell,birthTime,sep=":"),X],
            #                                 mutantRotatedY[paste(thisSister,birthTime,sep=":"),X]-mutantRotatedY[paste(thisCell,birthTime,sep=":"),X],
            #                                 mutantRotatedZ[paste(thisSister,birthTime,sep=":"),X]-mutantRotatedZ[paste(thisCell,birthTime,sep=":"),X])
            #                    thisUnitVector=thisVector/sqrt(sum(thisVector^2))
            #                    as.double(WTUnitVector %*% thisUnitVector)
            #                    }))
            
            
            mags=sqrt((mutantRotatedX[paste(thisSister,birthTime,sep=":"),]-mutantRotatedX[paste(thisCell,birthTime,sep=":"),])^2+
                      (mutantRotatedY[paste(thisSister,birthTime,sep=":"),]-mutantRotatedY[paste(thisCell,birthTime,sep=":"),])^2+
                      (mutantRotatedZ[paste(thisSister,birthTime,sep=":"),]-mutantRotatedZ[paste(thisCell,birthTime,sep=":"),])^2)
            
            
            Mutant_Dot=rbind(Mutant_Dot, 
                             (mutantRotatedX[paste(thisSister,birthTime,sep=":"),]-mutantRotatedX[paste(thisCell,birthTime,sep=":"),])/mags*WTUnitVector[1]+
                             (mutantRotatedY[paste(thisSister,birthTime,sep=":"),]-mutantRotatedY[paste(thisCell,birthTime,sep=":"),])/mags*WTUnitVector[2]+
                             (mutantRotatedZ[paste(thisSister,birthTime,sep=":"),]-mutantRotatedZ[paste(thisCell,birthTime,sep=":"),])/mags*WTUnitVector[3])
            
            WTmags=sqrt((WTRotatedX[paste(thisSister,birthTime,sep=":"),]-WTRotatedX[paste(thisCell,birthTime,sep=":"),])^2+
                        (WTRotatedY[paste(thisSister,birthTime,sep=":"),]-WTRotatedY[paste(thisCell,birthTime,sep=":"),])^2+
                        (WTRotatedZ[paste(thisSister,birthTime,sep=":"),]-WTRotatedZ[paste(thisCell,birthTime,sep=":"),])^2)
            
            WT_Dot=rbind(WT_Dot, 
                         (WTRotatedX[paste(thisSister,birthTime,sep=":"),]-WTRotatedX[paste(thisCell,birthTime,sep=":"),])/WTmags*WTUnitVector[1]+
                         (WTRotatedY[paste(thisSister,birthTime,sep=":"),]-WTRotatedY[paste(thisCell,birthTime,sep=":"),])/WTmags*WTUnitVector[2]+
                         (WTRotatedZ[paste(thisSister,birthTime,sep=":"),]-WTRotatedZ[paste(thisCell,birthTime,sep=":"),])/WTmags*WTUnitVector[3])
            
            #            rbind(WT_Dot,sapply(colnames(WTRotatedX),function(X){
                #                    thisVector=c(mutantRotatedX[paste(thisSister,birthTime,sep=":"),X]-mutantRotatedX[paste(thisCell,birthTime,sep=":"),X],
            #                                 mutantRotatedY[paste(thisSister,birthTime,sep=":"),X]-mutantRotatedY[paste(thisCell,birthTime,sep=":"),X],
            #                                 mutantRotatedZ[paste(thisSister,birthTime,sep=":"),X]-mutantRotatedZ[paste(thisCell,birthTime,sep=":"),X])
            #                    thisUnitVector=thisVector/sqrt(sum(thisVector^2))
            #                    as.double(WTUnitVector %*% thisUnitVector)
            #                    }))
            
            RowNames=rbind(thisParent,RowNames)
            
            
            
            #Get XYZ positions of this cell and sister in each movie in WT
            #Get XYZ positions of this cell and sister in each mutant movie
            #Calculate mean vector in WT and angles (mean and SD)
            #Calculate mean vector in mutant and angles relative to WT (dot scores)
            }
        }
        
    rownames(WT_Dot)=RowNames[,1]
    rownames(Mutant_Dot)=RowNames[,1]
     
    WT_Mean=rowMeans(WT_Dot,na.rm=T)
    WT_SD=apply(WT_Dot,1,sd,na.rm=T)

    Mutant_Mean=rowMeans(Mutant_Dot,na.rm=T)
    Mutant_SD=apply(Mutant_Dot,1,sd,na.rm=T)
    pVals=sapply(names(WT_Mean), function(X){
            obj<-try(t.test(WT_Dot[X,],Mutant_Dot[X,]),silent=T)
            if (is(obj, "try-error")) return(NA) else return(obj$p.value)})

    
    write.csv(data.frame(WT_Mean,WT_SD,Mutant_Mean,Mutant_SD,pVals,peak[names(WT_Mean)],Mutant_Dot),file=paste(Name,"_dots.csv",sep=""))


    setwd(startDir)

}
