source("functions.R")
source("AnalyzeDivTimes.R")
source("AnalyzePositions.R")
source("AnalyzeRotation.R")
source("PlotDefectSummaries.R")
source("PlotPositionDevs.R")

args <- commandArgs(trailingOnly = TRUE)
print(args)


moviesList = args[1]
print(moviesList)
print(length(args))
exp_file <- args[2]
peak=ReadPeakExpression(exp_file)
expCutoff = as.integer(args[3])
sig = as.integer(args[4])
microns=as.integer(args[5])
posDevTime = as.integer(args[6])

StartDir=getwd()
cc_dev = AnalyzeDivTimes(moviesList)

setwd(StartDir)
devs=AnalyzePositions(moviesList, Expression=peak, CalculateNeighbors=TRUE, expCutoff = expCutoff)
setwd(StartDir)
AnalyzeRotation(moviesList, peak=peak)
setwd(StartDir)
PlotDefectSummaries(moviesList, peak,sig=sig,microns=microns,expCutoff=expCutoff,minDivTime=70)
setwd(StartDir)

MeanPosDevs = PlotDeviationsList(Name=moviesList,exp=peak, t=posDevTime)
PlotExpVsDev(moviesList,outfile=paste(moviesList,"ExpVsDev.pdf",sep="."),exp=peak)

