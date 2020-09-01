#!/bin/bash


cp ../MARGEData_mm/margeFactor_mm.csv .
#python run_uDHS_count_h5_genome_window1kb.py H3K27ac
#python run_uDHS_count_h5_genome_window1kb.py DNase

#python run_uDHS_count_h5_genome_window1kb.py H3K4me3

# not done
#python run_uDHS_count_h5_genome_window1kb.py H3K4me1
#python run_uDHS_count_h5_genome_window1kb.py ATAC-seq
#python run_uDHS_count_h5_genome_window1kb.py H3K4me2
#python run_uDHS_count_h5_genome_window1kb.py H3K27me3
#python run_uDHS_count_h5_genome_window1kb.py H3K9me3
#python run_uDHS_count_h5_genome_window1kb.py H3K9ac
python run_uDHS_count_h5_genome_window1kb.py H3K36me3

