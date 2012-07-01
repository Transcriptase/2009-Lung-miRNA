#Install Packages
library("AgiMicroRna")

#Pull the internal read.agiMicroRna function
SCRIPT_PATH <- "~/GitHub/2009-Lung-miRNA/2009-lung-mirna"
setwd(SCRIPT_PATH)
source("read.agiMicroRna.R")
source("tgsMicroRna.R")
source("fixed_esetMicroRna.R")

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
                       annotation = c( "ControlType", "ProbeName","SystematicName"),
                       verbose=TRUE)

#Because we got the compact and not the full output files, the RMA
#normalization method is not available. Going with a process that
#mimics AFE's approach using Total Gene Signal instead.

ddTGS = tgsMicroRna(dd, half = TRUE, makePLOT = FALSE,verbose = FALSE)
ddNORM = tgsNormalization(ddTGS, "quantile", makePLOTpre = FALSE,
                          makePLOTpost = FALSE, targets,
                          verbose = TRUE)
ddPROC = filterMicroRna(ddNORM, 
                        dd, 
                        control = TRUE, 
                        IsGeneDetected = TRUE,
                        wellaboveNEG = FALSE, 
                        limIsGeneDetected = 50, 
                        limNEG = 25,
                        makePLOT = FALSE, 
                        targets, 
                        verbose = TRUE,
                        writeout = FALSE)
esetPROC = fixed_esetMicroRna(ddPROC, targets, makePLOT = FALSE, verbose = TRUE)

PVcut=0.05
Mcut=0.0    	 
MTestmethod = "BH"
DEmethod="separate"

levels.treatment=levels(factor(targets$Treatment))
treatment=factor(as.character(targets$Treatment),
                 levels=levels.treatment)
design=model.matrix(~ -1 + treatment)

contrast_matrix=cbind(CTRvsCR=c(0,-1, 0, 0, 0, 1, 0, 0))

fit2 = basicLimma(esetPROC, design, contrast_matrix, verbose = TRUE)
DE = getDecideTests(fit2, 
                    DEmethod = DEmethod, 
                    MTestmethod = MTestmethod, 
                    PVcut = PVcut, 
                    verbose = TRUE)
RESULTS_PATH = "G:/Project - Lung miRNA Expression Analysis/Results"
setwd(RESULTS_PATH)
write.fit(fit2, results = DE, "CTRvsCT.txt", adjust = "BH")

miRNAs_of_interest <- c("mmu-miR-9",
                        "mmu-miR-15b",
                        "mmu-miR-27",
                        "mmu-miR-29a",
                        "mmu-miR-30a",
                        "mmu-miR-103",
                        "mmu-miR-107",
                        "mmu-miR-155",
                        "mmu-miR-194",
                        "mmu-miR-200a",
                        "mmu-miR-200b",
                        "mmu-miR-200b*",
                        "mmu-miR-200c",
                        "mmu-miR-205",
                        "mmu-miR-204",
                        "mmu-miR-221",
                        "mmu-miR-222",
                        "mmu-miR-661",
                        "mmu-miR-7",
                        "mmu-miR-10a",
                        "mmu-miR-10b",
                        "mmu-miR-16",
                        "mmu-miR-17",
                        "mmu-miR-21",
                        "mmu-miR-22",
                        "mmu-miR-31",
                        "mmu-miR-122",
                        "mmu-miR-126",
                        "mmu-miR-146a",
                        "mmu-miR-146b",
                        "mmu-miR-194",
                        "mmu-miR-206",
                        "mmu-miR-214",
                        "mmu-miR-335",
                        "mmu-miR-378",
                        "mmu-miR-520c")
mat.DE <- matrix(DE)
rownames(mat.DE) <- rownames(DE)
miRNA_interest_results <- which(rownames(dat.DE) %in% miRNAs_of_interest)
sig_miRNA_of_interest <- rownames()