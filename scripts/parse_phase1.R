library(gdata)


meta = read.xls('phase1_tables3.xls', skip=1, header=T)

test = subset(meta, Manipulated.TF=='Aes') # & Logratio.of.expression.change >= 1.589 & FDR.for.expression.change <=0.01)

for (i in unique(meta[,1])) {
    # fold change 1.5
    cat(as.character(subset(meta, Manipulated.TF==i & abs(Logratio.of.expression.change) >= 0.5 &FDR.for.expression.change <= 0.05)$Affected.genes), file=paste0(i, '.symbol'), sep='\n')
}

