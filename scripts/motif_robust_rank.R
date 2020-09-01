#install.packages('RobustRankAggreg')
library(RobustRankAggreg)
args <- commandArgs(T)

tf <<- args[1]
histone <<- args[2]

cutoffs <- 97:99
batch_delta_test <-function(cutoff) {
    df <- paste0(tf, '_delta_batch_', histone,'_withpromoter_dnase_', cutoff, '.txt')
    content <- read.table(df, skip=1, header=T)
    diff <- read.table(df, nrow=1)
    diff <- as.numeric(diff[1,])==1
    
    ks = apply(content, 1, function(x) ks.test(x[diff], x[!diff], alternative='less')$p.value)
    
    print(cutoff)
    fig = paste0(tf, '_', cutoff, '_motif_deltaRP_FDR.pdf')
    result = -log(sort(ks, decreasing=F), 10)
    pdf(fig, width=6, height=30)
    par(mar=c(3, 10, 1, 1), las=2, font=2)
    barplot(result[1:100], main="top 100 motif by KS test", col='blue', horiz=T)
    dev.off()
    return(ks)
}

glist <- apply(matrix(cutoffs), 1, batch_delta_test)
glist <- apply(glist, 1, min)

result = -log(sort(glist, decreasing=F), 10)
fig = paste0(tf, '_final_motif_deltaRP_FDR.pdf')
pdf(fig, width=6, height=30)
par(mar=c(3, 10, 1, 1), las=2, font=2)
barplot(result[1:100], main="top 100 motif by KS test", col='blue', horiz=T)
dev.off()

