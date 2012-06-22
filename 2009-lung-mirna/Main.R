#Install Packages
library("AgiMicroRna")

#Pull the internal read.agiMicroRna function
SCRIPT_PATH <- "~/GitHub/2009-Lung-miRNA/2009-lung-mirna"
setwd(SCRIPT_PATH)
source("read.agiMicroRna.R")

#Read target file
TARGET_PATH <- "~/GitHub/2009-Lung-miRNA/data/sample_groups.txt"
targets <- readTargets(TARGET_PATH)

#Read the data files
DATA_PATH <- "M:/MRSTran/Lab Notebook/Array/Lung 2009/Raw Data/miRNA 12-2009/P.Tran_MicroRNA_012210"
setwd(DATA_PATH)

#Usual workflow would use readMicroRnaAFE(), which is a wrapper for
#read.maimages in the limma package, but because the files are in 
#compact rather than full format, I have to re-wrap the underlying function
#based on advice from the package author in
#https://www.stat.math.ethz.ch/pipermail/bioconductor/2010-August/035136.html
# and http://permalink.gmane.org/gmane.science.biology.informatics.conductor/28127
# and http://www.mentby.com/Group/bioconductor-list/agimicrorna-problems.html



dd <- read.agiMicroRna(targets,
                       columns=list(TGS="gTotalGeneSignal",
                                    TPS="gTotalProbeSignal",
                                    #meanS="gMeanSignal",
                                    procS="gProcessedSignal"),
                       other.columns=list(IsGeneDetected="gIsGeneDetected",
                                          IsSaturated="gIsSaturated",
                                          IsFeatNonUnifOF="gIsFeatNonUnifOL",
                                          IsFeatPopnOL="gIsFeatPopnOL",
                                          BGKmd="gBGMedianSignal",
                                          BGKus="gBGUsed"),
                       annotation = c( "ControlType", "ProbeName","GeneName"),
                       verbose=TRUE)