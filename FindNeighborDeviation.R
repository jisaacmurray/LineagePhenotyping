

1)Find all cells within threshold distance T of each cell at each time point
dimensions(cellsXtimepointsXcells)

2)Make the same list for each individual WT embryo
-for each cell there are some % of neighbors are also within threshold
-find an inflated threshold that includes all neighbors (this could be either a global inflated threshold, or cell-specific)

3)Find cells in mutant that fail the inflated threshold.  AND find cells that have new neighbors that are inappropriately close.




Johns alternative

1)Calculate distance matrix of each cell to each other cell at each time point. (all time points?)

2)Repeat for each WT, mutant, reference

3)Convert mutant matrix or each WT matrix to p-values for each cell pair

4)Collapse to single value per cell

2)For each WT series, find distance of each predicted neighbor

Joshs alternative:

library("mlegp")
XRef=MeanPositionFrame[,3]
YRef=MeanPositionFrame[,4]
YRef=MeanPositionFrame[,5]
notNA=!is.na(tXpositions[,1])

plot(XRef[notNA],tXpositions[notNA,1])
mlegp(XRef[notNA],tXpositions[notNA,1])
]
dev=(tXpositions-)^2
tYdev=(tYpositions-MeanPositionFrame[,4])^2
tZdev=(tZpositions-MeanPositionFrame[,5])^2
