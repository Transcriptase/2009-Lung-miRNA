#functions for data visualization
name_prep <- function(name){
    while (length(grep("-", name)) > 0 || length(grep("\\*", name)) > 0){
        name <- sub("-", ".", name)
        name <- sub("\\*", ".", name)
    }
    return(name)
}

miRNA_lookup <- function(miRNA, samples){
    data <- df.eset[samples, name_prep(miRNA)]
    return(data)
}

make_table <- function(miRNA_list, samples){
    names_list <- as.vector(miRNA_list)
    miRNA_data <- sapply(names_list, miRNA_lookup, samples)
    miRNA_mat <- t(data.matrix(miRNA_data))
    return(miRNA_mat)
}
