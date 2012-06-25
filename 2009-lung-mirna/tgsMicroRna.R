tgsMicroRna <- 
function (dd, offset = 0, half = TRUE, makePLOT = FALSE, verbose = FALSE) 
{
    if (!is(dd, "uRNAList")) {
        stop("'input' must be a uRNAList")
        if (is.null(dim(dd)[1])) {
            stop("'input' is empty")
        }
    }
    TT = dim(dd)[1]
    uniqueProbe = unique(dd$genes$ProbeName)
    LUP = length(uniqueProbe)
    uniqueGene = unique(dd$genes$SystematicName)
    LUG = length(uniqueGene)
    if (verbose) {
        cat("\n")
        cat("GETTING Agilent Feature Extraction TotalGeneSignal", 
            "\n")
        cat("\n")
        cat("\tTotal Probes:\t", TT, "\n")
        cat("\tUnique Probe: \t", LUP, "\n")
        cat("\tUnique Gene: \t", LUG, "\n")
        cat("\n")
    }
    ug = which(duplicated(dd$genes$SystematicName) == FALSE)
    nGEN = length(ug)
    ddTGS = dd[ug, ]
    if (half) {
        for (i in 1:dim(ddTGS)[2]) {
            index = which(ddTGS$TGS[, i] < 0.5)
            ddTGS$TGS[index, i] = 0.5
        }
    }
    else {
        min = min(ddTGS$TGS)
        for (i in 1:dim(ddTGS)[2]) {
            ddTGS$TGS[, i] = ddTGS$TGS[, i] + (abs(min) + offset)
        }
    }
    ddTGS$TGS = ddTGS$TGS
    ddTGS$TGS = ddTGS$TGS
    ddTGS$TGS = ddTGS$TGS
    ddTGS$TGS = ddTGS$TGS
    nARR = dim(ddTGS)[2]
    geneNames = list(c(dd$genes$SystematicName[ug]), c(1:nARR))
    if (!missing(makePLOT)) {
        if (makePLOT) {
            MMM = log2(ddTGS$TGS)
            maintitle = "TotalGeneSignal"
            colorfill = "green"
            dev.new()
            boxplotMicroRna(MMM, maintitle, colorfill)
            dev.new()
            plotDensityMicroRna(MMM, maintitle)
            dev.new()
            ddaux = ddTGS
            ddaux$meanS = MMM
            mvaMicroRna(ddaux, maintitle, TRUE)
            rm(ddaux)
        }
    }
    return(ddTGS)
}