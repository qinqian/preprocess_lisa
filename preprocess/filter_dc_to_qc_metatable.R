#!/usr/bin/env Rscript

meta <- read.delim('DC_samples_id.xls')

qc <- read.delim('DC_samples_quality.xls', check.names=F)

qc <- cbind(meta[match(qc$SampleIDs, meta$id), c("species_id", "unique_id", "series_id", "cell_type_id", "cell_line_id", "tissue_type_id", "IF", "platform")], qc)

rownames(qc) <- qc$SampleIDs
qc <- qc[,-9]
colnames(qc) <- gsub("_id", "", colnames(qc))

write.csv(qc, file="dc_meta_qc.csv", quote=F)

margeFactor <- rbind(subset(qc, species=='Homo sapiens'&FactorName=='H3K4me3'),
subset(qc, species=='Homo sapiens'&FactorName=='DNase'),
subset(qc, species=='Homo sapiens'&FactorName=='ATAC-seq'),
subset(qc, species=='Homo sapiens'&FactorName=='H3K4me1'),
subset(qc, species=='Homo sapiens'&FactorName=='H3K4me2'),
subset(qc, species=='Homo sapiens'&FactorName=='H3K9ac'),
subset(qc, species=='Homo sapiens'&FactorName=='H3K9me3'),
#subset(qc, species=='Homo sapiens'&FactorName=='H3'),
#subset(qc, species=='Homo sapiens'&FactorName=='H3K79me2'),
#subset(qc, species=='Homo sapiens'&FactorName=='H4K20me1'),
#subset(qc, species=='Homo sapiens'&FactorName=='H3K9me2'),
subset(qc, species=='Homo sapiens'&FactorName=='H3K27me3'),
subset(qc, species=='Homo sapiens'&FactorName=='H3K36me3'),
#subset(qc, species=='Homo sapiens'&FactorName=='H3'),
subset(qc, species=='Homo sapiens'&FactorName=='H3K27ac'))

print(table(as.vector(margeFactor$FactorName)))
print('a')
#print(head(sort(table(as.vector(subset(qc, FactorType=='hm'&species=='Homo sapiens')$FactorName)), decreasing=T), 15))
write.csv(margeFactor, quote=F, file='margeFactor.csv')
