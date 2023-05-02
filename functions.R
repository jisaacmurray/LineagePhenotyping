wideScreen <- function(howWide=Sys.getenv("COLUMNS")) {
  options(width=as.integer(howWide))
}
#wideScreen()

GetPeak<-function(Expression,X){
    Parent=GetParent(X)
    MyExp=Expression[X,]
    if(is.na(Parent)){
        if(is.na(MyExp)){
            return(0)
        }else{
            return(MyExp)
        }
    }
    
    ParentPeak=GetPeak(Expression,Parent)
    if(is.na(MyExp)){
        return(GetPeak(Expression,Parent))
        }
    if(is.na(ParentPeak)){
        return(MyExp)
        }
    return(max(MyExp,ParentPeak))
}

ReadPeakExpression <- function(FileName){
    exp = read.csv(FileName,row.names=2)["blot"]
    return(sapply(Cells,function(X){
        GetPeak(exp,X)
        }))
}


GetSister <- function(x){
    if(is.na(x)) return(NA)		
    if(substr(x,nchar(x),nchar(x))=="a"){
        return(paste(substr(x,1,nchar(x)-1),"p",sep=""))
        }
    if(substr(x,nchar(x),nchar(x))=="d"){
        return(paste(substr(x,1,nchar(x)-1),"v",sep=""))
        }
    if(substr(x,nchar(x),nchar(x))=="l"){
        return(paste(substr(x,1,nchar(x)-1),"r",sep=""))
        }	
    if(substr(x,nchar(x),nchar(x))=="p"){
        return(paste(substr(x,1,nchar(x)-1),"a",sep=""))
        }
    if(substr(x,nchar(x),nchar(x))=="v"){
        return(paste(substr(x,1,nchar(x)-1),"d",sep=""))
        }
    if(substr(x,nchar(x),nchar(x))=="r"){
        return(paste(substr(x,1,nchar(x)-1),"l",sep=""))
        }	



    if(x=="AB") return("P1")
    if(x=="EMS") return("P2")
    if(x=="C") return("P3")
    if(x=="D") return("P4")
    if(x=="Z2") return("Z3")
    if(x=="E") return("MS")

    if(x=="P1") return("AB")
    if(x=="P2") return("EMS")
    if(x=="P3") return("C")
    if(x=="P4") return("D")
    if(x=="Z3") return("Z2")
    if(x=="MS") return("E")

    return(NA)
}

GetParent <- function(x){
    if(is.na(x) || x=="P") return(NA)
    if(is.element(substr(x,nchar(x),nchar(x)),c("a","p","d","v","l","r"))) return (substr(x,1,nchar(x)-1))

    if(x=="AB" || x=="P1") return("P0")

    if(x=="EMS" || x=="P2") return("P1")
    if(x=="E" || x=="MS") return("EMS")
    if(x=="C" || x=="P3") return("P2")
    if(x=="D" || x=="P4") return("P3")
    if(x=="Z2" || x=="Z3") return("P4")
    return("P")

}

SeqArray <- function(x,y){
    myArray=c()
    for(n in names(x)){
        if(!is.na(x[n])&&!is.na(y[n])){
            myArray=c(myArray,seq(x[n],y[n])) 
            }
        else{
            myArray=c(myArray,c())
            }
        }
    return (myArray)
}
