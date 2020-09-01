fs <- c(list.files("../PhaseA_motif_deltarp/noDNase/30sample/", ".*coef", full.names=T),
       list.files("../PhaseC_newsystem/noDNase/30sample/", ".*coef", full.names=T))

fsn <- c(list.files("../PhaseA_motif_deltarp/noDNase/30sample/", ".*coef", full.names=F),
        list.files("../PhaseC_newsystem/noDNase/30sample/", ".*coef", full.names=F))

fsn <- gsub("tables1_il4","IL4", fsn)
fsn <- gsub("tables1_stat6","STAT6_siRNA_IL4", fsn)
fsn <- gsub("Macrophage2Kdo26h","PU1", fsn)

n <<-1
result <- lapply(fs, function(x) {
    yt <- unlist(strsplit(fsn[n], "_"))[1]
    xt <- unlist(strsplit(unlist(strsplit(fsn[n], "_"))[2], '.',fixed=T))[1]
    print(yt)
    print(xt)
    y=read.table(x,sep='\t', row.names=NULL, header=T, stringsAsFactors=F)
    y <- y[order(abs(y[,2]), decreasing=T),,drop=F]
    y <- cbind(cbind(rep(xt, nrow(y)), rep(yt, nrow(y))), y)
    colnames(y) <- c("TF", "Epigenomics", "Sample", "Coefficient")
    n<<- n+1
    return(y)
})

result <- Reduce(rbind,result)
result[,3] <- sapply(result[,3], function(y) {return(unlist(strsplit(unlist(strsplit(y,"_"))[2], "-"))[1])})
result[,3] <- gsub("MCF", "MCF-7", result[,3])
result[,3] <- gsub("LP", "LP-1", result[,3])

pdf("9cases_features.pdf", width=6, height=10)
par(mfrow=c(9,7), font=2, cex=0.7, las=2)
combs <- unique(result[,1:2])
tfs <- unique(combs[,1])
hms <- unique(combs[,2])
n <- 1
for (tf in tfs) {
    for (hm in hms) {
    sub.mat <- subset(result, TF==tf & Epigenomics==hm)
    sub.mat <- sub.mat[order(sub.mat[,4], decreasing=T),]
    z <- sub.mat[1:2,4]
    names(z) <- sub.mat[1:2,3]
    if (n%%7==1){
        yl=tf
        lmar=1
        xp=T
        }
    else{
        yl=""
        lmar=0.1
        xp=F
        }
    if (n>7){
        par(mar=c(0,lmar,0,0.1), xpd=xp)
        mi=""
        }
    else {
        par(mar=c(0,lmar,3,0.1), xpd=xp)
        mi=hm
        }
    pos <- barplot(z, col='grey', main=mi, ylab="", xlab="", width=0.3, names.arg=NULL,axes=F, axisnames = F)
    bd <- par('usr')
    text(pos, bd[3]*(-1.1), labels=names(z), las=2, adj=0,srt=90, cex=0.75,font=2)
    text(0, (bd[3]+bd[4])/2, labels=yl, las=2, adj=0,srt=90, cex=0.8)
#    axis(2, las=2, cex.axis=0.6)
    n <- n+1
    }
    }
dev.off()
system("cp 9cases_features.pdf ~/public_html")
