args <- commandArgs(T)
tf <- args[1]
cutoff <- args[2]
histone <- args[3]

df <- paste0(tf, '_delta_batch_', histone, '_withpromoter_dnase_', cutoff, '.txt')
content <- read.table(df, skip=1, header=T)
diff <- read.table(df, nrow=1)

diff <- as.numeric(diff[1,])==1

ks = apply(content, 1, function(x) ks.test(x[diff], x[!diff], alternative='less')$p.value)

#print(ks[grep("AR", names(ks))])
print(ks)

fig = paste0(tf, '_', cutoff, '_motif_deltaRP_FDR.pdf')
pdf(fig, width=6, height=30)
par(mar=c(3, 10, 1, 1), las=2, font=2)
barplot(-log(sort(ks, decreasing=F)[1:200], 10), main="top 50 motif KS test FDR", col='blue', horiz=T)
dev.off()
