source("functions.R")

GetPeak<-function(Expression,X){
    Parent=GetParent(X)
    MyExp=Expression[X,]
    if(is.na(Parent)){
        return(MyExp)
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

