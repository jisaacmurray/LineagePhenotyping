# LineagePhenotyping
This document describes the procedure for identifying defects in embryos by lineaging

## Dependencies
R with rgl, gdata, plotrix, ggplot2, plotly, plot3D
Perl with Statistics::Descriptive, Math::Trig
 
 
Local copies of:
AnalyzeDivTimes.R
AnalyzePositions.R
AnalyzeRotation.R
PlotDefectSummaries.R
PlotPositionDevs.R
Functions.R
MakeDB.pm

## Example workflow (assuming murraylab embryoDB instance in /gpfs/fs0/l/murr/tools3/embryoDB/
 
1.     Edit, check, extract, reopen in acetree, resave, re-extract each movie

2.     Make list of movies in directory (text file with one movie name per line)
a.     Can at this point re-extract all movies using sageExtract.pl or sageNotextract.pl with the path to the list as the argument)
b.     E.g. in the tools3 directory: “sageExtract.pl /gpfs/fs0/l/murr/lists/myMutants”

3.     Use PrintTrees to make sure the trees/editing all look good (checking that you have your partial editing codes specified correctly in embryoDB etc).
a.     Make sure you are in the tools3 directory
b.     E.g. “PrintTrees.pl /gpfs/fsp0/l/murr/lists/myMutants -200 500 rainbow 5”

4.     GetACD.pl (with the list as argument)
a.     Make sure you are in the tools3 directory
b.     E.g. “GetACD.pl /gpfs/fs0/l/murr/lists/myMutants”

5.     Run RunAllAnalysis.pl
a.     This will download summary files of cell positions and division times to a local subdirectory with the same name as the list you use as an argument
b.     Make sure you are in the directory that has the phenotyping scripts
c.     E.g. “GetACD.pl /gpfs/fs0/l/murr/lists/myMutants”

6.     Optional: Copy the “CA” file for a movie with expression of the gene you are interested in (otherwise you can use one of the provided CA files in the “data” subdirectory – this is relevant for outputs related to whether defects are enriched in expressing cells. 

7.     Run LineagePhenotyping.R
a.     Make sure you are in the directory that has the phenotyping scripts. There should be a subdirectory with the same name as your list that was created by RunAllAnalysis.pl… Arguments are
                                               i.     Name of list
                                             ii.     Path to CA file with expression values
                                            iii.     Cutoff to use for expression calculation
                                            iv.     Significance cutoff for defects (z score)
                                              v.     Magnitude cutoff for position defects (microns)
                                            vi.     Time point to use for plotting expression vs deviation
b.     E.g. “Rscript LineagePhenotyping.R ceh-32_mutant data/CA_CEH-32_SYS85.csv 200 3 3 250 > ceh-32_mutant/ceh-32_mutant.stats.txt myMutants”

8.     Output files are described separately in this sheet; below is a suggested approach for prioritization since there are many possible downstream questions and analyses.



## Question 1: Is the rate of development comparable between
Quick: Look at the plots of divTime and ccLength in each mutant vs the WT mean in the “CC_plots.pdf” file
More involved: look in the “STATS.tsv” file for differences in rate within each lineage. Analyze the CCLength.tsv and DivTime.tsv files 

## Question 2: 
Plots in CC_plots.pdf may give some hints (how many cells with high Z scores or low p values are there, approximately where are they in the lineage?)
“CellDefectSummary.csv” file will give a quick view of whether there are any cells with large numbers of division defects (recommend sorting by “TotalDefects” column)
VERY IMPORTANT to iteratively QC this analysis to identify any residual errors:
Open the “CCdev.txt” fils
(Optional: just filter to focus on the largest defects e.g. increasing the threshold from the default Z>3 and |deltaCC|>5 minutes)
Sort by embryo
Open each embryo in AceTree and manually check the cells identified as having defects by 
