args <- commandArgs(T)
tf <- args[1]
histone <- args[2]

df <- paste0(tf, '_delta_batch_', histone, '_withpromoter_dnase_cistrome_dc.txt')
diff <- read.table(df, sep='\t', nrow=1)
content <- read.delim(df, skip=1, header=T, sep='\t')

diff <- as.numeric(diff[1,])==1

ks = apply(content, 1, function(x) ks.test(x[diff], x[!diff], alternative='less')$p.value)

#print(ks)

fig = paste0(tf, '_cistromedc_motif_deltaRP_FDR.pdf')
pdf(fig, width=18, height=20)
par(mar=c(3, 30, 1, 1), las=2, font=2, cex=1.5)
barplot(-log(sort(ks, decreasing=F)[1:80], 10), main="top 80 Cistrome TF ChIP-seq deletion KS test", col='black', horiz=T)
dev.off()

