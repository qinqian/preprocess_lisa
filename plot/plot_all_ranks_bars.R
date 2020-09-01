## fs <- c(list.files("../PhaseA_motif_deltarp/noDNase/30sample/", ".*perfo.*json", full.names=T),
##        list.files("../PhaseC_newsystem/noDNase/30sample/", ".*perfo.*json", full.names=T))
## fsn <- c(list.files("../PhaseA_motif_deltarp/noDNase/30sample/", ".*perfo.*json", full.names=F),
##         list.files("../PhaseC_newsystem/noDNase/30sample/", ".*perfo.*json", full.names=F))

args = commandArgs(T)
tf = args[1]

library(gplots)
library(RobustRankAggreg)

plot_rra <- function(x, folder, motif=True, sp='hg38', active=T) {
print('ok')
    if (motif){
        fs <- list.files(folder, paste0("^", tf, ".*ks"), full.names=T)
        fs <- fs[grep("cistrome", fs, invert=T)]
        meta = read.delim('/data/home/qqin/seqpos2/cistrome.txt', stringsAsFactors=F)
        fsn <- list.files(folder, paste0("^", tf, ".*ks"), full.names=F)
        fsn <- fsn[grep("cistrome", fsn, invert=T)]
        marks = gsub(".txt.ks", "", gsub("_withpromoter_nodnase", "", gsub(paste0(tf, "_delta_batch_"), "", fsn)))
    } else {
        fs <- list.files(folder, paste0("^", tf, ".*dc.*ks"), full.names=T)
        if (sp=='hg38') {
            meta = read.delim('/data/home/qqin/01_Projects/Programming/dc2/scripts/best_dc_tfcr_basedon_frip_peak_dhs.xls', stringsAsFactors=F)
        } else {
            meta = read.delim('/data/home/qqin/01_Projects/Programming/dc2/scripts/best_dc_tfcr_basedon_frip_peak_dhs_mm.xls', stringsAsFactors=F)
            }
        fsn <- list.files(folder, paste0("^", tf, ".*dc.*ks"), full.names=F)
        marks = gsub("_withpromoter_nodnase_cistrome_dc.txt.ks", "", gsub(paste0(tf, "_delta_batch_"), "", fsn))
}

dc.ks = lapply(fs, function(x)
{
    y=read.csv(x, header=F,row.names=1)
    rownames(y)[order(y[,1])]
})

names(dc.ks) = marks
dc.rra = aggregateRanks(dc.ks)

    if (motif) {
dc.rra$factor = meta[match(dc.rra$Name, meta[,1]), "symbol"]
cl = gsub(" Family", "", meta[match(dc.rra$Name, meta[,1]), "dbd"])
dc.rra$cell <- cl
    }
    else {
dc.rra$factor = meta[match(dc.rra$Name, meta[,1]), "FactorName"]
cl = meta[match(dc.rra$Name, meta[,1]), "cell_line"]
cl[is.na(cl)] = ""
dc.rra$cell <- paste0(cl, " ", meta[match(dc.rra$Name, meta[,1]), "cell_type"])
    }

dc.ks = lapply(fs, function(x)
{
    read.csv(x, header=F,row.names=1)
})

dc.p = Reduce(cbind,dc.ks)

rown = apply(cbind(rownames(dc.p), dc.rra[match(rownames(dc.p), rownames(dc.rra)),3]), 1, function(x) {paste0(x[1], '_', x[2])})
dc.final=dc.p
rownames(dc.final) = rown
colnames(dc.final) = marks

n=1061

dc.plot.rank = apply(dc.final,2,function(x) {rank(x)})

dc.plot = dc.plot.rank[grepl('PPARG', rownames(dc.plot.rank), ignore.case =T), ]
dc.plot = rbind(dc.plot, dc.plot.rank[grepl('CEBP', rownames(dc.plot.rank), ignore.case =T),])

color = colorRampPalette(c("white", "red"))

    library(gplots)

    if (motif) {
        motifl="motif"
    } else {
        motifl = ""
    }
pdf(paste0(tf, "_", motifl, "_rra_dc_targetone_barplot.pdf"), width=20, height=12)
dc.plot = t(dc.plot)
par(xpd=T, mar=c(1,12,8,3), cex=1.5, font=2)
image(x=1:nrow(dc.plot), y=1:ncol(dc.plot), z=dc.plot/1061, col=color(n),
      xlab="", ylab="", axes=F)
for (i in 1:nrow(dc.plot)) {
for (j in 1:ncol(dc.plot)) {
    text(i, j,  labels=round(dc.plot[i,j],1), cex=0.8, col='darkblue')
}

}
axis(3, at=1:nrow(dc.plot), labels=rownames(dc.plot), tick=F, las=2, line=-0.3)
axis(2, at=1:ncol(dc.plot), labels = colnames(dc.plot), las=2, tick=F)
dev.off()
system(paste0("cp ", tf, "_", motifl, "_rra_dc_targetone_barplot.pdf ~/public_html"))
}

# l1
plot_rra(tf, "../PhaseA_motif_deltarp/noDNase/10sample_wider_lambda/", T, active=F)
## plot_rra(tf, "../PhaseA_motif_deltarp/noDNase/10sample_l2/", T, active=T) # l2 not work for AR???

## l1
## plot_rra(tf, "../PhaseC_newsystem//noDNase/10sample_widerlambda/", T, 'mm10', active=F)
##plot_rra(tf, "../PhaseC_newsystem//noDNase/10sample_l2/", F, 'mm10', active=F)

