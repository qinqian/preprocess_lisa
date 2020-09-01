#!/bin/bash

#mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
#        "select chrom, size from hg38.chromInfo" > hg38.genome

# human
# 1000bp windows
#bedtools makewindows -g hg38.genome -w 1000 | awk '{OFS="\t"; print $1,$2,$3,NR}'> hg38_window1kb.bed

# 100bp windows
#bedtools makewindows -g hg38.genome -w 100 | awk '{OFS="\t"; print $1,$2,$3,NR}'> hg38_window100bp.bed


## NOTE:: intersect  hg38_window1kb.bed and hg38_window100bp.bed to get the 100bp to 1kb mapping

# to consider motif boundary cases for 100bp windows, +/- 10bp
#bedtools makewindows -g hg38.genome -w 100 | bedtools slop -b 10 -i - -g hg38.genome | awk '{OFS="\t"; print $1,$2,$3,NR}'> hg38_window100bp_both10bp.bed

mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
        "select chrom, size from mm10.chromInfo" > mm10.genome
# 1000bp windows
bedtools makewindows -g mm10.genome -w 1000 | awk '{OFS="\t"; print $1,$2,$3,NR}'> mm10_window1kb.bed

# 100bp windows
bedtools makewindows -g mm10.genome -w 100 | awk '{OFS="\t"; print $1,$2,$3,NR}'> mm10_window100bp.bed

# to consider motif boundary cases for 100bp windows, +/- 10bp
bedtools makewindows -g mm10.genome -w 100 | bedtools slop -b 10 -i - -g mm10.genome | awk '{OFS="\t"; print $1,$2,$3,NR}'> mm10_window100bp_both10bp.bed
