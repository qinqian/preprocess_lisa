tf <- commandArgs(T)

## tf <-  'AR'
##tf <-  'ESR1'
## tf <-  'PPARG'
## tf <-  'NOTCH'

#asig <- read.csv(paste0(tf, "_H3K27ac_H3K27ac_cluster_aver_signal.csv"))

pdf(paste0(tf, "_kl_rank.pdf"), width=10, height=8)
#for (hm in c("H3K27ac", "H3K4me1", "H3K4me3", "H3K4me1", "H3K9ac", "DNase", "ATAC-seq")) {
for (hm in c("H3K27ac", "H3K4me3", "H3K27me3", "DNase", "ATAC-seq")) {
kl.all <- read.csv(paste0(tf, "_", hm, "_", hm, "_cluster_stat_kl_divergence.csv_kl_rank.csv"))

klplot <- function(kl, motif=T) {
    kl[,3] <- rank(kl[,3])

    if (motif) {
        #labels = kl[grepl(paste0("^", tf, '$'), kl[,1], ignore.case=T),]
        # labels = rbind(labels, kl[grepl(paste0("^", "CEBPB", "$"), kl[,1], ignore.case=T),])
        labels = kl[grepl(paste0("^", "RBPJ", "$"), kl[,1], ignore.case=T),]
        # its cofactor
        # labels = rbind(kl[grepl("^ESR2$", kl[,1], ignore.case=T),], labels)
        ## labels = rbind(kl[grepl("^NR3C1$", kl[,1], ignore.case=T),], labels)
        ## labels = rbind(kl[grepl("^FOXA1$", kl[,1], ignore.case=T),], labels)
    } else {
        labels = kl[grepl("chipseq", kl[,1], ignore.case=T),]
        #labels = rbind(labels, kl[grepl("_ESR1_", kl[,1], ignore.case=T),])
        #labels = rbind(labels, kl[grepl("_ESR2_", kl[,1], ignore.case=T),])

        labels = rbind(labels,kl[grepl("_NOTCH1_", kl[,1], ignore.case=T),])
        labels = rbind(kl[grepl("_NOTCH_", kl[,1], ignore.case=T),], labels)
        labels = rbind(kl[grepl("_RBPJ_", kl[,1], ignore.case=T),], labels)

    #labels = rbind(labels, kl[grepl("_PPARG_", kl[,1], ignore.case=T),])
    #    labels = rbind(kl[grepl("_CEBPB_", kl[,1], ignore.case=T),], labels)

#labels = kl[grepl("_AR_", kl[,1], ignore.case=T),]
#labels = rbind(kl[grepl("_NR3C1_", kl[,1], ignore.case=T),], labels)
    }

    par(las=2, font=2, cex=1.2)
    ## hist(kl[,2], n=20, main=paste0("2556 CistromeDB ChIP-seq peak and 1061 Cistrome motif distribution KL divergence \n across ", "10 clusters of 10 selected H3K27ac samples (", tf, ")"), col='grey', xlab='KL Divergence', ylab='Frequency', axes=F)
    if(motif) {
        hist(kl[,2], n=20, main=paste0("1061 Cistrome motif distribution KL divergence \n across ", "10 clusters of 10 selected ",hm," samples (", tf, ")"), col='grey', xlab='KL Divergence', ylab='Frequency', axes=F)
    } else {
        hist(kl[,2], n=20, main=paste0("2556 CistromeDB ChIP-seq peak distribution KL divergence \n across ", "10 clusters of 10 selected ",hm," samples (", tf, ")"), col='grey', xlab='KL Divergence', ylab='Frequency', axes=F)
    }
axis(2, las=2)
axis(1, las=1)
for (i in 1:nrow(labels)) {
    abline(v=labels[i,2], col='blue', lty=2, lwd=2)
    if (i%%2==0)
        text(labels[i,2], 3, paste0(gsub("ChIP_seq_","",as.vector(labels[i,1])), " rank ",labels[i,3]), col="red", srt=90,cex=1.0, adj=0,font=2)
    else
        text(labels[i,2], 80, paste0(gsub("ChIP_seq_","",as.vector(labels[i,1])), " rank ",labels[i,3]), col='red', srt=90, cex=1.0, adj=0, font=2)
}
}

ci <- c(grep('ChIP_seq', kl.all[,1]), grep('chipseq', kl.all[,1]))
kl.chip <- kl.all[ci,]
kl.motif <- kl.all[-ci,]


klplot(kl.chip, F)
klplot(kl.motif,T)
}
    dev.off()
system(paste0("cp ", tf, "_kl_rank.pdf ~/public_html"))

#klplot(kl.all)
