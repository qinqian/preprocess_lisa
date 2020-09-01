fs <- list.files("/data/home/qqin/public_html/creeds/mm10", pattern='.*cistrome.*ks', full=T)
meta <- read.delim("~/01_Projects/Programming/dc2/scripts/dc_samples.xls", sep='\t')

test <- lapply(fs, read.csv, header=F)

test <- Reduce(function(x,y,...) { return(merge(x,y, by='V1')) }, test)

ann <- as.vector(meta[match(test[,1], meta[,1]), 'factor'])

fi <- apply(test[,2:ncol(test)], 2, function(x) {return(rank(x,ties.method ='min', na.last=T))})

fs <- apply(matrix(fs), 1, function(x) {
    gsub("_delta_batch", "",gsub("gene_", "", gsub("_withpromoter_nodnase_cistrome_dc.txt.ks", "", basename(x))))
})

fn <- apply(matrix(fs), 1, function(x) {
    gsub("-.*", "", gsub("_H3K27ac", "", gsub("_DNase", "", gsub("\\d*_", "", x))))
})

N <- 2000  # some magic number, possibly an overestimate
DF <- data.frame(id=rep("", N), tf=rep("", N), status=rep("", N),  # as many cols as you need
                hm=rep("", N), rank=rep(NA,N),
                stringsAsFactors=FALSE)          # you don't know levels yet
z <- 1
for (i in 1:ncol(fi)) {
    val=fi[grep(fn[i],ann, ignore.case=T),i]
    if (length(val)>0) {
        z <- z+1
        DF$hm[i] <- unlist(strsplit(fs[i], split='_'))[3]
        DF$tf[i] <- unlist(strsplit(unlist(strsplit(fs[i], split='_'))[2], split='-'))[1]
        DF$status[i] <- unlist(strsplit(unlist(strsplit(fs[i], split='_'))[2], split='-'))[2]
        DF$id[i] <- unlist(strsplit(fs[i], split='_'))[1]
        DF$rank[i] <- min(as.numeric(val))
    }
}

DF <- na.omit(DF)
DF$hm <- paste0(DF$hm, '_ChIP-seq')

##write.table(DF, file="CREEDS_mm10_cistrome_dc_rank.txt", quote=F, sep='\t', row.names=F)


fs <- list.files("/data/home/qqin/public_html/creeds/mm10", pattern='.*ks', full=T)

fs <- fs[grep("cistrome", fs, invert=T)]

meta <- read.delim("~/seqpos2/cistrome.txt", sep='\t')

test <- lapply(fs, read.csv, header=F)

test <- Reduce(function(x,y,...) { return(merge(x,y, by='V1')) }, test)

ann <- as.vector(meta[match(test[,1], meta[,1]), 'symbol'])

fi <- apply(test[,2:ncol(test)], 2, function(x) {return(rank(x,ties.method ='min', na.last=T))})

fs <- apply(matrix(fs), 1, function(x) {
    gsub("_delta_batch", "",gsub("gene_", "", gsub("_withpromoter_nodnase_99.txt.ks", "", basename(x))))
})

fn <- apply(matrix(fs), 1, function(x) {
    gsub("-.*", "", gsub("_H3K27ac", "", gsub("_DNase", "", gsub("\\d*_", "", x))))
})

N <- 2000  # some magic number, possibly an overestimate
DF2 <- data.frame(id=rep("", N), tf=rep("", N), status=rep("", N),  # as many cols as you need
                hm=rep("", N), rank=rep(NA,N),
                stringsAsFactors=FALSE)          # you don't know levels yet

z <- 1
for (i in 1:ncol(fi)) {
    val=fi[grep(fn[i],ann, ignore.case=T),i]
    if (length(val)>0) {
        z <- z+1
        DF2$hm[i] <- unlist(strsplit(fs[i], split='_'))[3]
        DF2$tf[i] <- unlist(strsplit(unlist(strsplit(fs[i], split='_'))[2], split='-'))[1]
        DF2$status[i] <- unlist(strsplit(unlist(strsplit(fs[i], split='_'))[2], split='-'))[2]
        DF2$id[i] <- unlist(strsplit(fs[i], split='_'))[1]
        DF2$rank[i] <- min(as.numeric(val))
    }
}

DF <- na.omit(rbind(DF, DF2))


DF$tf_id = paste0(DF$tf, "_", DF$id)
DF <- DF[order(DF$tf),]

DF[DF$status=='dn', 3] = 'Down-regulated genes'
DF[DF$status=='up', 3] = 'UP-regulated genes'

DF$hm[grep('ChIP-seq', DF$hm, invert=T)] = paste0(DF$hm[grep('ChIP-seq', DF$hm, invert=T)], '_Motif')

write.table(DF, file="CREEDS_mm10_chip_and_motif_rank.txt", quote=F, sep='\t', row.names=F)

library(ggplot2)
p <- ggplot(DF) + geom_bar(aes(x=factor(hm), y=rank, fill=factor(status)), stat='identity', position='dodge')
p <- p+facet_wrap(~tf_id, scales = "free", ncol=8)## + ylim(c(0,100))
p <- p+theme(legend.position = "top",
          legend.text = element_text(size=8, family="sans", face="bold", colour="black"),
          strip.text = element_text(size=10, family="sans", face="bold", colour="blue"),
          legend.title = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size=8, angle=90, family="sans", face="bold", colour="black"),
          axis.text = element_text(size=8, family="sans", face="bold", colour="black"),
          title = element_text(size=8, family="sans", face="bold", colour="black")) +
     ggtitle("TF Rank from CREEDS human differential expression gene list")+xlab("Epigenomics feature") + ylab("TF Rank")
ggsave("mm10_creeds_dnase_h3k27ac_chipseq_motif_rank.pdf", height=100, limitsize=F, width=12)


## pdf("test_ecdf.pdf")
## plot(ecdf(subset(DF, status=='Down-regulated genes' & hm =='DNase')$rank), col='purple')
## lines(ecdf(subset(DF, status=='UP-regulated genes' & hm =='H3K27ac')$rank), col='red')
## lines(ecdf(subset(DF, status=='Down-regulated genes' & hm =='H3K27ac')$rank), col='blue')
## lines(ecdf(subset(DF, status=='UP-regulated genes' & hm =='DNase')$rank), col='black')
## dev.off()
