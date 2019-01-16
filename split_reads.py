import pysam
from math import log10
from random import shuffle
import timeit
from sortedcontainers import SortedDict, SortedList, SortedSet
import math

from variant_candidate import VariantCandidate
from variant_candidate_collection import VariantCandidateCollection
from read_pair import ReadPair

DEL_TYPE = 0
INS_TYPE = 1
DUP_TYPE = 2
INV_TYPE = 3
TRA_TYPE = 4

var_types = {
    DEL_TYPE: "DEL",
    INS_TYPE: "INS",
    DUP_TYPE: "DUP",
    INV_TYPE: "INV",
    TRA_TYPE: "TRANS"
    }
  
                
def to_phred(x):
    return -10*(log10(x))

def from_phred(x):
    return pow(10, -x/10)

def process_reads(my_reads):
    for read in my_reads:
        read_is_good = not (read.is_duplicate or 
                            read.is_qcfail or 
                            read.mapping_quality == 0 or 
                            read.is_secondary or 
                            read.is_supplementary or 
                            read.is_unmapped or
                            read.reference_id == None) and (read.is_paired)
        if read_is_good:
            read_id = read.query_name
            print read.cigartuples
            #print read

            

def call_variants(input_file, 
                  reference_file, 
                  output_file,
                  my_start,
                  my_stop,
                  my_chr, 
                  print_results=True):
    
    ref =  pysam.FastaFile(reference_file).fetch(region="{}:{}-{}".format(my_chr, my_start + 1, my_stop + 1))
    samfile = pysam.AlignmentFile(input_file, "rb")
    reads = samfile.fetch(region="{}:{}-{}".format(my_chr, my_start + 1, my_stop + 1))
    my_reads = [a for a in reads]
    #shuffle(my_reads, lambda: 0.4)
    #shuffle(my_reads)
    
    process_reads(my_reads)
   

def call_test():
    call_variants(my_start = 59000,
                  my_stop = 65000,
                  my_chr = "20",
                  input_file = "/Users/siakhnin/data/giab/mnist_na12878_chrom20.bam",
                  reference_file = "/Users/siakhnin/data/reference/genome.fa",
                  output_file = None)
    
#print timeit.timeit(call_test,number=1)
call_test()
#in_vcf = pysam.VariantFile("/Users/siakhnin/data/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.freebayes.vcf", "r")
#out_vcf = pysam.VariantFile("/Users/siakhnin/data/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.rheos.vcf", "w", header=in_vcf.header)