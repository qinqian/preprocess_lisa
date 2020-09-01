z <- read.csv("single_gene_perturbations-v1.0.csv", stringsAsFactors=F)

mmeta <- read.delim("~/seqpos2/cistrome.txt", stringsAsFactors=F)

z.tfh <- z[z$hs_gene_symbol %in% mmeta$symbol, ]
z.tfm <- z[z$mm_gene_symbol %in% mmeta$symbol, ]

z.tf <- rbind(z.tfh, z.tfm)
z.tf <- unique(z.tf[order(z.tf[,1]),])


write.csv(z.tf, file="single_gene_perturbations-v1.0.csv.tfs", quote=F, row.names=F)

con <- file('single_gene_perturbations-v1.0.gmt', 'r')
line <- readLines(con)

i=0
while(i < length(line)) {
    i <- i+1
    newTxt <- strsplit(line[i], split = "\\t")[[1]]

    tf <- unlist(strsplit(newTxt[1], split= '-'))[1]
    istf <- match(newTxt[2], z.tf[,1])
    if (is.na(istf)){
        next
        }
    if (z.tf[istf,'organism'] == 'mouse')
        cat(newTxt[3:length(newTxt)], file=paste0(sub(":", "_", newTxt[2]), '_', newTxt[1], '.mm10'), sep='\n')
    else if (z.tf[istf,'organism'] == 'human')
        cat(newTxt[3:length(newTxt)], file=paste0(sub(":", "_", newTxt[2]), '_', newTxt[1], '.hg38'), sep='\n')
    else
        print('non')
}

close(con)
