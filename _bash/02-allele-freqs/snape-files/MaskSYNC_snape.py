
"""
From 
https://github.com/DEST-bio/DEST_freeze1/blob/main/mappingPipeline/scripts/MaskSYNC_snape_complete.py
with edits by Lucas A. Nell (lines 129 and 137, plus comments above each)
"""

import sys
from collections import defaultdict as d
import gzip
from optparse import OptionParser, OptionGroup
import pickle

# Author: Martin Kapun
# Modified by Maria Bogaerts

#########################################################   HELP   #########################################################################
usage="""python %prog \
      --sync input.sync \
      --coverage input.cov \
      --indel input.indel \
      --te TE.gff \
      --mincov 10 \
      --maxcov 0.9 \
      --output output """
parser = OptionParser(usage=usage)
group=OptionGroup(parser,
"""
H E L P:
____________

This script masks sites with too low/high coverage or which are located in proximity (within 5bp) to InDels or within TE's. In addition, it generates a BED file, with coordinates of all masked sites. It requires the indel/cov output files from Mpileup2Sync.py and a GFF from Repeatmaser as input files. Modified to parse the fifth column with information from SNAPE-pooled.
""")

parser.add_option("--sync", dest="sync", help="The original sync file")
parser.add_option("--output", dest="out", help="Output file prefix")
parser.add_option("--te", dest="te", help="A GFF file containing the coordinates of TE's (optional) ",default="NA")
parser.add_option("--indel", dest="indel", help="Python Object file with indel positions in the input.sync")
parser.add_option("--coverage", dest="cov", help="Python Object file with the coverage distribution of the input.sync")
parser.add_option("--maxcov", dest="maxcov", help="The maximum coverage threshold percentile, e.g. 0.9 ")
parser.add_option("--mincov", dest="mincov", help="The minimum coverage threshold: e.g. 10")
parser.add_option("--SNAPE", action="store_true", default=False, dest="snape", help="Parsing SNAPE sync file")
parser.add_option("--maxsnape", dest="maxsnape", help="e.g. 0.9", type="float")


parser.add_option_group(group)
(options, args) = parser.parse_args()

def load_data(x):
    if x=="-":
        y=sys.stdin
    elif x.endswith(".gz"):
        y=gzip.open(x,"rt")
    else:
        y=open(x,"rt")
    return y

def sync2string(x):
    ''' convert sync format to string of nucleotides  where x is a string in sync format '''
    string=""
    if x==".:.:.:.:.:." or x=="0:0:0:0:0:0":
        return ""
    alleles=["A","T","C","G"]
    ah={ k:v for (k,v) in zip(alleles,[int(X) for X in x.split(":")[:4]])}
    #ah=dict(zip(alleles,[int(X) for X in x.split(":")[:4]]))
    for k,v in ah.items():
        string+=v*k
    return string

## make dict of positions to mask
exclude=d(lambda:d(str))

## calculate max cov TH
with open(options.cov, 'rb') as coverage_file:
    allcoverage = pickle.load(coverage_file)

maximumcov=d(int)

for chrom,covh in allcoverage.items():
    tot=0
    if list(set(covh))==[0]:
        maximumcov[chrom]=0
        continue
    t=0
    for c,val in sorted(covh.items()):
        if c==-1:
            continue
        elif tot/float(covh[-1])>float(options.maxcov):
            t=1
            maximumcov[chrom]=c
            break
        else:
            tot+=val
    if t==0:
        maximumcov[chrom]=c

if options.indel!="NA":
    with open(options.indel, 'rb') as indel_file:
        indel = pickle.load(indel_file)
    for chrom,indel1 in indel.items():
        for pos,tup in indel1.items():
            for x,y in tup:
                for i in range(int(pos)+x,int(pos)+y+2):
                    exclude[chrom][i]

if options.te!="NA":
    for l in load_data(options.te):
        if l.startswith("#"):
            continue
        C=l.split()[0]
        S,E=map(int,l.split()[3:5])
        for i in range(S,E+1):
            exclude[C][i]

Start=""
RR=0
CR=""
# BED=gzip.open(options.out+".bed.gz","wt")
# SO=gzip.open(options.out+"_masked.sync.gz","wt")
BED=open(options.out+".bed","wt")
SO=open(options.out+"_masked.sync","wt")

if options.snape:
    for l in load_data(options.sync):
        C,P,R,S,I=l.split()
        # From SNAPE-pooled docs, `p_n0` is 
        # "$1 - p (0)$ where $p (f)$ is the probability distribution function 
        # for the minor allele freqeuncy"
        p_n0 = float(I.split(":")[4])
        COV=len(sync2string(S))
        """
        I changed the following line bc it was originally returning mostly
        monomorphic SNPs.
        The paper this filtering was based on did it the way I did below (see 
        https://github.com/GonzalezLab/SNP_caller_benchmarking/blob/ccfe4135b79b2d24bde48de46189bd790cd15ddb/simulations/snape2_simulations.sh#L83).
        """
        if int(P) in exclude[C] or COV < int(options.mincov) or COV > maximumcov[C] or p_n0 < float(options.maxsnape):
            b = SO.write("\t".join([C,P,R,".:.:.:.:.:."])+"\n")
            if Start=="":
                Start=int(P)
                RR=int(P)
                CR=C
            elif RR < int(P)-1 or CR!=C:
                b = BED.write("\t".join([CR,str(Start-1),str(RR)])+"\n")
                Start=int(P)
                RR=int(P)
                CR=C
            else:
                RR=int(P)
        else:
            b = SO.write("\t".join([C,P,R,S,I])+"\n")
    if Start!=RR:
        b = BED.write("\t".join([CR,str(Start-1),str(RR)])+"\n")
else:
    for l in load_data(options.sync):
        C,P,R,S=l.split()
        COV=len(sync2string(S))
        if int(P) in exclude[C] or COV < int(options.mincov) or COV > maximumcov[C]:
            b = SO.write("\t".join([C,P,R,".:.:.:.:.:."])+"\n")
            if Start=="":
                Start=int(P)
                RR=int(P)
                CR=C
            elif RR < int(P)-1 or CR!=C:
                b = BED.write("\t".join([CR,str(Start-1),str(RR)])+"\n")
                Start=int(P)
                RR=int(P)
                CR=C
            else:
                RR=int(P)
        else:
            b = SO.write("\t".join([C,P,R,S])+"\n")
    if Start!=RR:
        b = BED.write("\t".join([CR,str(Start-1),str(RR)])+"\n")
