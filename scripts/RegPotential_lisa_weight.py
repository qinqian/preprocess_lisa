#!/usr/bin/env python2.7
"""
LISA weight with all peaks fold change >= 5, input macs xls files
"""
import sys, os, time
import math
from optparse import OptionParser

# original Score calc function
#Sg = lambda ldx: sum([math.exp(-0.5-4*t) for t in ldx])
# lisa weight function
Sg = lambda ldx: sum([2*math.exp(-10.0*t*math.log(3))/(1+math.exp(-10.0*t*math.log(3))) for t in ldx])

# print current time and information on screen
def Info(infoStr):
    print("[%s] %s" %(time.strftime('%H:%M:%S'), infoStr))

def prepare_optparser():
    """
    Prepare optparser object and validation.
    New options will be added in this function first.
    """
    usage = "usage: %prog ......"
    description = "Input a peak file, and ...... "
    optparser = OptionParser(version="%prog v1.00", description=description, usage=usage, add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    optparser.add_option("-t","--treat",dest="peakfile",type="string",
                         help="Input the MACS's result peak file(.xls), it will recognize it by extension.")
    optparser.add_option("-n","--name",dest="name",type="string",
                         help="this argument is used to name the result file.")
    optparser.add_option("-d","--distance", dest="distance", type="int",
                         help="Set a number which unit is 'base'. It will get peaks within this distance from gene TSS. default:100000 (100kb)", default=100000)
    optparser.add_option("-g","--genome",dest="genome",type="string",
                         help="Select a genome file (sqlite3 file) to search refGenes.")
    optparser.add_option("--top",dest="top",type="float",
                         help="input a number between 0-1, so that the script will only output a percentage of top genes.\
                               input a number bigger than 1, for example, 2000. so that the script will only output top 2000 genes. Default: output All ~30000 genes", default = -1)

    (options,args) = optparser.parse_args()
    if not options.peakfile and not options.genome and not options.name:
        optparser.print_help()
        sys.exit(1)
    if not os.path.isfile(options.peakfile):
        Info('ERROR: Cannot find peak file, a tab-peak file must be given through -t (--treat).')
        sys.exit(1)

    if not os.path.isfile(options.genome):
        Info("ERROR: Genome file not found! A annottion file must be given through -g (--genome).")
        sys.exit(1)

    if not options.name:
        options.name = os.path.splitext(options.peakfile)[0] + "_result"

    if options.top > 1:
        options.output_flag = "number"
    elif 0 < options.top < 1:
        options.output_flag = "percentage"
    elif options.top == -1:
        options.output_flag = "all"
    else:
        Info("--top options error, please set a number 0-1 or bigger than 1")
        sys.exit(1)

    # print arguments
    Info("Argument List: ")
    Info("Name = " + options.name)
    Info("peak file = " + options.peakfile)
    Info("distance = %d bp" %options.distance)
    Info("genome = %s" %options.genome)
    Info("top = %f" %options.top)
    return options

class PScore:
    # connect to sqlite, select genome.
    def __init__(self, options):
        self.peakfile = options.peakfile
        self.genome = options.genome
        self.db = open(self.genome)
        self.c = self.db.readlines()
        self.db.close()
        self.geneInfo = []
        self.output_flag = options.output_flag
        self.top = options.top
        self.opts_string = "# Argument List:\n" +\
                           "# Name = %s\n" %options.name +\
                           "# peak file = %s\n" %options.peakfile +\
                           "# distance = %d bp\n" %options.distance +\
                           "# top = %f\n" %options.top
        self.peaklist = {}

    def readfile(self): #reads the file and returns a vector: each element is a bed_row.
        peakf = open(self.peakfile)
        from collections import defaultdict
        self.peaklist = defaultdict(list)

        import os
        bed = open("%s_5fold.bed" % os.path.basename(self.peakfile).split('.')[0], 'w')
        for line in peakf:
            if line.startswith("#") or (not line.strip()): #skip "#" lines and empty lines
                continue
            # .xls
            # chr	start	end	length	abs_summit	pileup	-log10(pvalue)	fold_enrichment	-log10(qvalue)	name
            line = line.strip().split()
            if line[0] == 'chr':
                continue

            if float(line[7]) < 5:
                continue
            bed.write("%s\t%s\t%s\n"%(line[0],int(line[4]), int(line[4])+1))
            # chr -> list with peak summit bp position
            self.peaklist[line[0]].append(int(line[4]))

        peakf.close()
        for i in self.peaklist:
            self.peaklist[i].sort()

    def ScoreCalc(self, distance):
        # use lisa histonerp .tss file
        # a line: chr1	30365	30366	NR_036267:MIR1302-10	0	+
        for line in self.c:
            if line.startswith("#"):
                continue
            line = line.strip().split()
            self.geneInfo.append([line[0], line[3].split(":")[0], int(line[1]), int(line[2]), line[5], line[3].split(":")[-1]]) # (0:chrom, 1:name, 2:txStart, 3:txEnd, 4:strand, 5:symbol)

        self.geneInfo.sort()
        self.geneInfo = [list(t) for t in self.geneInfo]

        count = 0
        for igene in self.geneInfo:
            if igene[4] == '+':
                gTSS = igene[2]
            elif igene[4] == '-':
                gTSS = igene[3]
            try:
                peaks = self.peaklist[igene[0]]
            except KeyError:
                peaks = []
            #peaksInDistance = [t+[abs((t[1]+t[2])/2-gTSS)] for t in peaks if t[4]>5 and abs((t[1]+t[2])/2-gTSS) < distance ]
            #peaksInDistance.sort(key=lambda x:x[-1])
            # peaks-> 0: summit, 1: fold change
            peaksInDistance = [abs(t-gTSS)*1.0/distance for t in peaks if abs(t-gTSS) <= distance]
            peaksInDistance.sort()
            igene.append(Sg(peaksInDistance))
            count += 1
            if not count % 2000:
                Info('Process <%d> genes'%count)
        self.geneInfo.sort(key=lambda x:x[-1], reverse=True)

    def Output2File(self, name):
        outf = open(name, "w")
        outf.write(self.opts_string)
        outf.write('#chrom\ttxStart\ttxEnd\trefseq\tscore\tstrand\tsymbol\n')
        for line in self.geneInfo:
            outf.write('%s\t%d\t%d\t%s\t%.3f\t%s\t%s\n'%(
                       line[0], line[2], line[3], line[1], line[6], line[4], line[5]))
        outf.close()
        Info("Finished! result output to <%s>"%name)

def main():
    opts=prepare_optparser()
    g = PScore(opts)
    g.readfile()
    g.ScoreCalc(opts.distance)
    g.Output2File(opts.name)

if __name__ == "__main__":
    main()
