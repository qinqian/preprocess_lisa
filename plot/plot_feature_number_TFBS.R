library(stringr)
args <- commandArgs(T)

tf <- args[1]
cell <- args[2]

clean.data <- function(x, tf)
{
    fs <- dir(x, pattern=paste0('^', tf, '.*top.*'))
print(fs)
    result <- lapply(fs, function(x) {
        read.table(x)[1,]
    })
    result <- Reduce(rbind, result)
    n <- as.numeric(sapply(fs, function(x) {substr(unlist(strsplit(x,"\\."))[4], 4,10)}))
    rownames(result) <- n
    result
}


result <- clean.data('.', tf)
result <- result[order(as.numeric(rownames(result))),]

pdf(paste0(tf, "_sample_size_effect.pdf"), width=5,height=4)
par(mar=c(6,5,5,1), font=2, cex=0.8, xpd=T)
plot(x=as.numeric(rownames(result)), y=result[,1], col='red', xlab="Sample Size",  yaxs='i', xaxs='i', lty=1, type='b', ylim=c(0,1), xlim=c(1,100), pch=16, axes=F, ylab="Performance")
points(x=as.numeric(rownames(result)), y=result[,2], col='darkblue', pch=16)
lines(x=as.numeric(rownames(result)), y=result[,2], col='darkblue', pch=16, lty=2)
axis(1, pos=0, at=seq(0, 110, 20), labels=seq(0, 110, 20), las=1, cex=0.9, font=2, tck=-0.02, line=1, lwd.ticks=0.5, padj=0)
axis(2, las=2, lwd.ticks=0.5, tck=-0.02, padj=0)
legend(2.8, 1.1, c("AUC", "PRAUC"), pch=c(16,16), col=c('red', 'darkblue'), lty=c(1,2),
       pt.cex=1, cex=0.6, ncol=2, bty='n')
title(paste0(tf, " Differential gene expression model\n directly predicts ", tf, " TF binding in ", cell), line=2.5)
dev.off()

