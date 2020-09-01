tfname=commandArgs(T)

fs <- list.files('.', pattern=paste0("^", tfname, ".*\\.bed\\.index.*motif.*"))

tf <- read.table(fs[1], row.names=1)
tfm <- read.table(fs[2], row.names=1)

pdf(paste0(tfname, "motif_chipseq_bar.pdf"), height=5,width=6)
par(las=2, cex=1, xpd=F)
n=1
## color =c(rgb(1, 1, 1, 0.25), rgb(1, 1, 0, 0.25), rgb(1, 0, 0, 0.25))
color = c("red", "blue", "black")
pchs=c(9, 10,16)
for (z in seq_along(rownames(tf))) {
d <- as.vector(as.matrix(tf[z,]))
dm <- as.matrix(tfm[z,])
names(d) <- colnames(tf)
names(dm) <- colnames(tfm)
print(d)
## hist(sort(d), n=100, main=rownames(tf)[z], col=)
## hist(sort(dm), n=100, main=rownames(tf)[z], col=rgb(1,0, 0, 0.2), add=T)
if (n==1)
    plot(d, dm, pch=pchs[n], col=color[n], xlim=c(0, 0.035), ylim=c(0,0.035), cex=0.45, xlab='Motif', ylab='Motif intersection with union DHS')
else
    points(d, dm, pch=pchs[n], col=color[n], cex=0.45)
text(d['AR']+0.001, dm['AR'], "AR", cex=0.4, col=color[n]) ## NR3C4 is AR
text(d['NR3C1']+0.0015, dm['NR3C1'], "NR3C1", cex=0.4, col=color[n])
text(d['Nr3c1']+0.001, dm['Nr3c1'], "Nr3c1",cex=0.4, col=color[n])
text(d['FOXA1']-0.001, dm['FOXA1']+0.0006, "FOXA1", cex=0.4, col=color[n])
text(d['Foxa1']-0.0023, dm['Foxa1']+0.001, "Foxa1", cex=0.4, col=color[n])
## text(d['Gata3']+0.0015, dm['Gata3'], "Gata3",cex=0.4, col=color[n])
## text(d['GATA3']+0.001, dm['GATA3'], "GATA3",cex=0.4, col=color[n])

## text(d['ESR1']+0.0015, dm['ESR1'], "ESR1",cex=0.4, col=color[n])
## text(d['Esr1']-0.0015, dm['Esr1'], "Esr1",cex=0.4, col=color[n])
## text(d['ESR2']+0.0015, dm['ESR2'], "ESR2",cex=0.4, col=color[n])
## text(d['Esr2']-0.0015, dm['Esr2'], "Esr2",cex=0.4, col=color[n])
grid(4,4)
n<-n+1
lines( par()$usr[1:2], par()$usr[3:4] )
}
legend("bottomright", paste0("Motif cutoff ", rownames(tf), " percentile"), col=color, pch=pchs, ncol=1, bty='n')
dev.off()
system(paste0("cp", " ", paste0(tfname, "motif_chipseq_bar.pdf"), " ~/public_html"))
