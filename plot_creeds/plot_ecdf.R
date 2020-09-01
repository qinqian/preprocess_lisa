z=read.table('NIA_shRNA_LISA_p_val_matrix.txt',header=T,stringsAsFactors = F)

library(reshape2)
library(RobustRankAggreg)

N <- 2000  # some magic number, possibly an overestimate
DF <- data.frame(id=rep("", N), tf=rep("", N),
                hm=rep("", N), rank=rep(NA,N),
                stringsAsFactors=FALSE)          # you don't know levels yet

zz <- 0
for (i in unique(z$gene_inducted)) {
    zsub = subset(z, gene_inducted == i)
    target.m = dcast(zsub, motif_gene_symbol~factor, min, value.var="p_val")
    target.m = target.m[, c("motif_gene_symbol", "DNase", "H3K27ac")]
    symbols = target.m[,1]
    ranks = apply(target.m[,c(2,3)],2,rank)
    index = grep(paste0('^', i, '$'),symbols,ignore.case=T)
    tfs = symbols[index]
    print(length(tfs))
    if (length(tfs) > 0) {
        for (m in 1:ncol(ranks)) {
            zz <- zz+1
            DF[zz, 1] = zz
            DF[zz, 2] = tfs[1]
            DF[zz, 3] = colnames(ranks)[m]
            DF[zz, 4] = ranks[j, m]
        }
    }
}

DF <- na.omit(DF)

z=read.table('Harmonizome_mouse_lof_LISA_p_val_matrix.txt',header=T,stringsAsFactors = F)

library(reshape2)
library(RobustRankAggreg)

N <- 2000  # some magic number, possibly an overestimate
DF2 <- data.frame(id=rep("", N), tf=rep("", N),
                 hm=rep("", N), rank=rep(NA,N),
                 stringsAsFactors=FALSE)          # you don't know levels yet

zz <- 0
for (i in unique(z$gene_inducted)) {
    zsub = subset(z, gene_inducted == i)
    target.m = dcast(zsub, motif_gene_symbol~factor, min, value.var="p_val")
    ## print(head(target.m))
    target.m = target.m[, c("motif_gene_symbol", "DNase", "H3K27ac")]
    ## print(head(target.m))
    print(i)
    symbols = target.m[,1]
    ranks = apply(target.m[,c(2,3)],2,rank)
    print(head(ranks))
    index = grep(paste0('^', i, '$'),symbols,ignore.case=T)
    tfs = symbols[index]
    if (length(tfs) > 0){
        for (m in 1:ncol(ranks)) {
            zz <- zz+1
            DF2[zz, 1] = zz
            DF2[zz, 2] = tfs[1]
            DF2[zz, 3] = colnames(ranks)[m]
            DF2[zz, 4] = ranks[j, m]
        }
    }
}

DF2 <- na.omit(DF2)

DF$project <- 'NIA_shRNA'
DF2$project <- 'harmonizome'



DF <- rbind(DF, DF2)
DF$hm <- paste0(DF$hm, "_Motif")
DF$tf_id <- paste0(DF$tf, '_', DF$id)

write.table(DF, file="Harmonizome_NIAshRNA_motif_rank.txt", quote=F, sep='\t', row.names=F)


crh <- read.delim("CREEDS_hg38_chip_and_motif_rank.txt", row.names=NULL)
crm <- read.delim("CREEDS_mm10_chip_and_motif_rank.txt", row.names=NULL)

plot.rank <- function(x) {
    ts <- unique(x$hm)
    par(mfrow=c(length(ts)/2, 2), las=2, cex=1.2, font=2, xpd=T)
    ## par(mfrow=c(1,2), las=2, cex=1.2, font=2, xpd=T)
    for (i in ts){
        y <- subset(x, hm==i)
        if (length(grep("ChIP", i))>0) {
            ## y[,3] <- cut(y[,3], breaks=c(0,10,30,50,100,200,300, 500, 1061, 2556))
            ## levels(y[,3]) <- c("0~10", "0~30", "0~50", "0~100", "0~200", "0~300", "0~500", "0~1061", "0~2556")
            y[,3] <- cut(y[,3], breaks=c(0,10,30,50,100,200,300, 500, 1604))
            levels(y[,3]) <- c("0~10", "0~30", "0~50", "0~100", "0~200", "0~300", "0~500", "0~1604")
        }else {
            y[,3] <- cut(y[,3], breaks=c(0,10,30,50,100,200,300, 500, 1061))
            levels(y[,3]) <- c("0~10", "0~30", "0~50", "0~100", "0~200", "0~300", "0~500", "0~1061")
        }
        y <- y[order(y[,3]),]
        y.s <- table(y[,3])
        val <- cumsum(y.s)/sum(y.s)
        xbar <- barplot(val, main=i, ylab=paste0('Cumulative proportion of TF rank'), xlab="TF Rank Range", axisnames=F, axes=F)
        axis(2, las=2, cex=1, font=2)
        text(xbar, -0.08, labels=names(val), srt=45, cex=1, font=2, adj=1)
        text(xbar, val+0.04, labels=paste0(round(val,3)*100, "%"), cex=0.6, font=2, pos=3)
        lines(xbar, val, lty=6, lwd=1, col='blue')
        points(xbar, val, pch=16, col='blue')
    }
    dev.off()
}

crh.dn <- subset(crh, status=='Down-regulated genes')[,c(4,6,5)]
crh.up <- subset(crh, status=='UP-regulated genes')[,c(4,6,5)]


plot.rank(DF.nia)
plot.rank(DF.har)

plot.rank(crh.dn)
length(unique(gsub("_\\d*", "", crh.dn[,2])))
length(unique(crh.dn[,2]))

plot.rank(crh.up)
length(unique(gsub("_\\d*", "", crh.up[,2])))
length(unique(crh.up[,2]))

crm.dn <- subset(crm, status=='Down-regulated genes')[,c(4,6,5)]
crm.up <- subset(crm, status=='UP-regulated genes')[,c(4,6,5)]

plot.rank(crm.up)
length(unique(gsub("_\\d*", "", crm.up[,2])))
length(unique(crm.up[,2]))

plot.rank(crm.dn)
length(unique(gsub("_\\d*", "", crm.dn[,2])))
length(unique(crm.dn[,2]))

DF.nia <- subset(DF, project=='NIA_shRNA')[,c(3,6,4)]

DF.har <- subset(DF, project=='harmonizome')[,c(3,6,4)]

plot.rank(DF.nia)
plot.rank(DF.har)

## total tfs
length(unique(toupper(c(gsub("_\\d*", "", DF[,2]), as.vector(crm[,2]), as.vector(crh[,2])))))
