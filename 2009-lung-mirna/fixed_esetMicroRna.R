fixed_esetMicroRna <- function (uRNAList, targets, makePLOT = FALSE, verbose = FALSE) 
{
    if (!is(uRNAList, "uRNAList")) {
        stop("'input' must be a uRNAList", call. = FALSE)
        if (is.null(dim(uRNAList)[1])) {
            stop("'input' is empty", call. = FALSE)
        }
    }
    if (missing(targets)) {
        stop("'targets' is missing ", call. = FALSE)
    }
    if ("GErep" %in% colnames(targets)) {
        GErep = targets$GErep
        nGE = sum(table(table(GErep)))
        g1 = targets$GErep
        g2 = rownames(targets)
    }
    else {
        stop("'targets' needs 'GErep' field")
    }
    goON = all(rownames(targets) == colnames(uRNAList$iTGS))
    if (!goON) {
        stop("rownames in pData(targets) different from colnames uRNAList$TGS", 
             call. = FALSE)
    }
    phenoData = new("AnnotatedDataFrame", data = targets)
    TMP = uRNAList$TGS
    rownames(TMP) = uRNAList$genes$SystematicName
    colnames(TMP) = g2
    nGEN = dim(TMP)[1]
    nARR = dim(TMP)[2]
    esetPROC = new("ExpressionSet", exprs = TMP, phenoData = phenoData)
    if (verbose) {
        cat("outPUT DATA: esetPROC", "\n")
        print(dim(esetPROC))
        cat("------------------------------------------------------", 
            "\n")
    }
    if (!missing(makePLOT)) {
        if (makePLOT) {
            dev.new()
            hierclusMicroRna(exprs(esetPROC), GErep, methdis = "euclidean", 
                             methclu = "complete", sel = FALSE, 100)
            dev.new()
            size = 100
            maintitle = "50 High Variance"
            HeatMapMicroRna(exprs(esetPROC), 50, maintitle)
            dev.new()
            PCAplotMicroRna(esetPROC, targets)
        }
    }
    return(esetPROC)
}