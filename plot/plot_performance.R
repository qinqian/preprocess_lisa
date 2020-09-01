fs <- c(list.files("../PhaseA_motif_deltarp/noDNase/30sample/", ".*perfo.*json", full.names=T),
       list.files("../PhaseC_newsystem/noDNase/30sample/", ".*perfo.*json", full.names=T))

fsn <- c(list.files("../PhaseA_motif_deltarp/noDNase/30sample/", ".*perfo.*json", full.names=F),
        list.files("../PhaseC_newsystem/noDNase/30sample/", ".*perfo.*json", full.names=F))

library(jsonlite)
fsn <- gsub("_performance.json", "", gsub("tables1_il4","IL4", fsn))
fsn <- gsub("_performance.json", "", gsub("tables1_stat6","STAT6_siRNA_IL4", fsn))
fsn <- gsub("_performance.json", "", gsub("Macrophage2Kdo26h","PU1", fsn))

system("cp 9cases_features.pdf ~/public_html")

AUC <- unlist(lapply(fs, function(x)
    fromJSON(x, flatten=T)[[2]][1]))
PRAUC <- unlist(lapply(fs, function(x)
    fromJSON(x, flatten=T)[[2]][2]))


names(AUC) <- fsn
names(PRAUC) <- fsn

o <- order(names(AUC))

AUC <- AUC[o]
PRAUC <- PRAUC[o]

pdf("9cases_exp_performance.pdf", width=8,height=5)
layout(mat=matrix(c(1,2,3), byrow=T), height=c(2.5,1.6,2))
par(cex=0.8, font=2, xpd=T,las=2, mar=c(0.1,2.5,3.5,0.1))
x <- barplot(AUC, names.arg=NULL, axisnames=F, axes=F, ylim=c(0,1), main="9 cases across 7 types of active epigenomics data", xpd=T, ylab="AUC")
axis(2, at=c(0,0.8,1), labels=c(0,0.8,1), cex.axis=0.7)
par(cex=0.8, font=2, xpd=T,las=2, mar=c(0,2.5,0,0.1))
plot(1, type='n', xlab="", ylab="",xlim=c(0,max(x)), ylim=c(0,0.3), axes=F)
text(x*0.995-0.02, 0.15, names(AUC), las=2, srt=90, adj=0.5, cex=0.5)
par(cex=0.8, font=2, xpd=T,las=2, mar=c(3,2.5,0,0.1))
barplot(-PRAUC, names.arg=NULL, axisnames=F, axes=F, ylim=c(-0.3,0), xpd=T, ylab="PRAUC")
axis(2, at=c(-0.3,0), labels=c(0.3,0), cex.axis=0.7)
dev.off()
system("cp 9cases_exp_performance.pdf ~/public_html")
