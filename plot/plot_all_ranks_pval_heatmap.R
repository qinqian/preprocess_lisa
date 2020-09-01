## fs <- c(list.files("../PhaseA_motif_deltarp/noDNase/30sample/", ".*perfo.*json", full.names=T),
##        list.files("../PhaseC_newsystem/noDNase/30sample/", ".*perfo.*json", full.names=T))
## fsn <- c(list.files("../PhaseA_motif_deltarp/noDNase/30sample/", ".*perfo.*json", full.names=F),
##         list.files("../PhaseC_newsystem/noDNase/30sample/", ".*perfo.*json", full.names=F))

args = commandArgs(T)
tf = args[1]

library(gplots)
library(RobustRankAggreg)

plot_rra <- function(x, folder, motif=True, sp='hg38', active=T) {
    if (motif){
        fs <- list.files(folder, paste0("^", tf, ".*ks"), full.names=T)
        fs <- fs[grep("cistrome", fs, invert=T)]
        print(fs)
        meta = read.delim('/data/home/qqin/seqpos2/cistrome.txt', stringsAsFactors=F)
        fsn <- list.files(folder, paste0("^", tf, ".*ks"), full.names=F)
        fsn <- fsn[grep("cistrome", fsn, invert=T)]
        marks = gsub(".txt.ks", "", gsub("_withpromoter_nodnase", "", gsub(paste0(tf, "_delta_batch_"), "", fsn)))
        print(marks)
    } else {
        fs <- list.files(folder, paste0("^", tf, ".*dc.*ks"), full.names=T)
        print(fs)
        if (sp=='hg38') {
            meta = read.delim('/data/home/qqin/01_Projects/Programming/dc2/scripts/best_dc_tfcr_basedon_frip_peak_dhs.xls', stringsAsFactors=F)
        } else {
            meta = read.delim('/data/home/qqin/01_Projects/Programming/dc2/scripts/best_dc_tfcr_basedon_frip_peak_dhs_mm.xls', stringsAsFactors=F)
            }
        fsn <- list.files(folder, paste0("^", tf, ".*dc.*ks"), full.names=F)
        marks = gsub("_withpromoter_nodnase_cistrome_dc.txt.ks", "", gsub(paste0(tf, "_delta_batch_"), "", fsn))
        print(marks)
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

dc.final=cbind(dc.p, dc.rra[match(rownames(dc.p), rownames(dc.rra)),c(2,3,4)])
dc.final = dc.final[order(dc.final$Score),]
colnames(dc.final) = marks

n=30
        dc.plot = dc.final[1:n,]
        if (!motif) {
            yl1 = dc.plot[,12]
            yl2 = dc.plot[,13]
        } else {
            yl1 = dc.plot[,32] ## three cutoffs
            yl2 = dc.plot[,33]
        }
    if (active){
        active.index=c(1:10)[-apply(matrix(c("K27me3", "K36me3", "K9me3")), 1, function(x) {grep(x, marks)})]
        print(active.index)
        z=7
    } else {
        if (motif){
            z=30
            active.index=1:z
        }else{
            active.index=1:10
            z=10
        }
    }

dc.plot = t(-log10(dc.plot[,active.index]))

color = colorRampPalette(c("white", "red"))

    library(gplots)

    if (motif) {
        motifl="motif"
    } else {
        motifl = ""
    }

    if (motif) {
        width=14
    }else {
        width=7
        }
print(head(dc.plot))
    if (n<=30)
        pdf(paste0(tf, "_", motifl, "_rra_dc.pdf"), width=width, height=8)
    else
        pdf(paste0(tf, "_", motifl, "_rra_dc.pdf"), width=width, height=15)

    #if (motif)
    #    par(xpd=T, mar=c(1,6,10,15))
    #else
    #    par(xpd=T, mar=c(1,6,6,15))
image(x=1:nrow(dc.plot), y=1:ncol(dc.plot), z=dc.plot, col=color(n),
      xlab="", ylab="", axes=F)
axis(3, at=1:nrow(dc.plot), labels=rownames(dc.plot), tick=F, las=2, line=-0.3)
axis(2, at=1:ncol(dc.plot), labels = yl1, las=2, tick=F)
axis(4, at=1:ncol(dc.plot), labels = yl2, las=2, tick=F)
box()

for (i in 1:n-0.5)
{
    for (j in 1:z-0.5) {
        rect(j,i,j+1,i+1, col=NA,border='grey')
    }
}


    if(motif)
        pos=32
    else
        pos=12


    if (n<=20) {
    rect(seq(pos,pos+4,length.out = n), 22, seq(pos,pos+4,length.out = n)+0.2,22.3, col = color(n),border='grey')
text(pos, 21.8, round(min(dc.plot),1), cex=.6)
text(pos+4, 21.8, round(max(dc.plot),1), cex=.6)
text(pos+2, 22.6, "-log10(p-value)",cex=0.8)
    } else
    {
    rect(seq(pos,pos+4,length.out = n), n+2, seq(pos,pos+4,length.out = n)+0.2, n+3, col = color(n),border=NA)
text(pos, n+1.8, round(min(dc.plot),1), cex=.6)
text(pos+4, n+1.8, round(max(dc.plot),1), cex=.6)
text(pos+2, n+3.6, "-log10(p-value)",cex=0.8)
        }
dev.off()
#grid(6,n, col='grey')

system(paste0("cp ", tf, "_", motifl, "_rra_dc.pdf ~/public_html"))
}

# l1
#plot_rra(tf, "../PhaseA_motif_deltarp/noDNase/10sample_wider_lambda/", T, active=F)
## plot_rra(tf, "../PhaseA_motif_deltarp/noDNase/10sample_l2/", T, active=T) # l2 not work for AR???


## l1
plot_rra(tf, "../PhaseC_newsystem//noDNase/10sample_widerlambda/", T, 'mm10', active=F)
##plot_rra(tf, "../PhaseC_newsystem//noDNase/10sample_l2/", F, 'mm10', active=F)
