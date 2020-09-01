zz <- read.csv("single_gene_perturbations-v1.0.csv.tfs")

gsms <- c()
for (i in 1:nrow(zz)){
    index=grep('GSM', as.vector(as.matrix(zz[i,])))
    print(index)
    for (j in index) {
        gsms <- c(gsms, unlist(strsplit(as.character(zz[i,j]), split='\\|')))
    }
}

length(unique(gsms))
