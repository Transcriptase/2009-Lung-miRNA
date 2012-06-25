`read.agiMicroRna` <-
function(targets,columns=NULL,other.columns = NULL,annotation=NULL,
          exp.names=NULL,verbose=FALSE){

## code adapted from read.maimages function from limma package 
## it creates the uRNAList class: see AgiMicroRNA-classes.R for 
## the definition of this class and associated methods  

filess=as.character(targets$FileName)

if (is.null(filess)) {
            stop("targets frame doesn't contain FileName column")
}

if (is.null(columns)){
  	columns=list(TGS="gTotalGeneSignal",
	             TPS="gTotalProbeSignal",
	  	     meanS="gMeanSignal",
		     procS="gProcessedSignal")
}

if (!is.list(columns)){
            stop("columns must be a list")
            }
 
if (is.null(other.columns)){
other.columns=list(IsGeneDetected="gIsGeneDetected",
				        IsSaturated="gIsSaturated",
				        IsFeatNonUnifOF="gIsFeatNonUnifOL",
				        IsFeatPopnOL="gIsFeatPopnOL",
				        BGKmd="gBGMedianSignal",
				        BGKus="gBGUsed")
} 

if (!is.list(other.columns)){
            stop("other.columns must be a list")
            }

if (is.null(exp.names)){
annotation = c( "ControlType", "ProbeName","SystematicName")
#changed "GeneName" to "SystematicName" to match file
#will need to track this change through future processing
}

if (is.null(exp.names)){
            exp.names = rownames(targets) 
}


cnames <- names(columns)
required.col = unique(c(annotation, unlist(columns), unlist(other.columns))) 


obj = read.columns(filess[1], required.col, text.to.search="",skip = 9, 
  sep = "\t",quote="\"", stringsAsFactors = FALSE,flush=TRUE)
narray = length(filess)
ngenes = dim(obj)[1]

# ngenes = length(scan(filess[1],skip=10,what="integer",flush=TRUE,quiet=TRUE))
Y = matrix(NA, nrow = ngenes, ncol = narray)
colnames(Y) =  exp.names 

# R, B, Rb, Gb 
Newagi = columns
for (a in cnames){
 Newagi[[a]] <- Y
}

# targets 
Newagi$targets = data.frame(targets$FileName)
rownames(Newagi$targets) = exp.names
colnames(Newagi$targets) = "FileName"

# $genes ("ControlType","ProbeName","SystematicName") 
j <- match(annotation, colnames(obj), 0)
if (any(j > 0)){
            Newagi$genes <- data.frame(obj[, j, drop = FALSE], check.names = FALSE)
}
# $other 
other.columns <- as.character(other.columns)
        j <- match(other.columns, colnames(obj), 0)
        if (any(j > 0)) {
            other.columns <- colnames(obj)[j]
            Newagi$other = list()
            for (j in other.columns) Newagi$other[[j]] <- Y
        }
        

for(n in 1:narray) {
 if(verbose){
  cat("reading file ",n," - ",filess[n],"\n")
  }
    if(n > 1){
    obj = read.columns(filess[n], required.col, text.to.search="",skip = 9, 
        sep = "\t",quote="\"", stringsAsFactors = FALSE,flush=TRUE)
    }
      
    for (a in cnames){
      Newagi[[a]][, n] <- obj[, columns[[a]]]    # R, G, Rb, Gb
    }
    
    for (j in other.columns) {
                Newagi$other[[j]][, n] <- obj[, j]
    }    

} ## for  

  # defined in AgiMicroRNA-classes.R USING classes.R limma FILE 
  new("uRNAList", Newagi)
   
} ## end function 
