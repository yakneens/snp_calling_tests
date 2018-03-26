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

inserts = {
        "20FUK.1":{"ReadSize":101,"Median":386,"MAD":19,"UniqueDiscordantPairs":17},
        "20FUK.2":{"ReadSize":101,"Median":414,"MAD":18,"UniqueDiscordantPairs":37},
        "20FUK.3":{"ReadSize":101,"Median":386,"MAD":19,"UniqueDiscordantPairs":8},
        "20FUK.4":{"ReadSize":101,"Median":413,"MAD":19,"UniqueDiscordantPairs":29},
        "20FUK.6":{"ReadSize":101,"Median":413,"MAD":18,"UniqueDiscordantPairs":29},
        "20FUK.5":{"ReadSize":101,"Median":386,"MAD":19,"UniqueDiscordantPairs":14},
        "20GAV.5":{"ReadSize":101,"Median":385,"MAD":19,"UniqueDiscordantPairs":21},
        "20FUK.7":{"ReadSize":101,"Median":386,"MAD":19,"UniqueDiscordantPairs":17},
        "20FUK.8":{"ReadSize":101,"Median":413,"MAD":18,"UniqueDiscordantPairs":26},
        "20GAV.4":{"ReadSize":101,"Median":413,"MAD":18,"UniqueDiscordantPairs":36},
        "20GAV.1":{"ReadSize":101,"Median":386,"MAD":19,"UniqueDiscordantPairs":15},
        "20GAV.2":{"ReadSize":101,"Median":414,"MAD":18,"UniqueDiscordantPairs":18},
        "20GAV.3":{"ReadSize":101,"Median":385,"MAD":19,"UniqueDiscordantPairs":24},
        "20GAV.6":{"ReadSize":101,"Median":413,"MAD":18,"UniqueDiscordantPairs":23},
        "20GAV.7":{"ReadSize":101,"Median":386,"MAD":19,"UniqueDiscordantPairs":25},
        "20GAV.8":{"ReadSize":101,"Median":413,"MAD":18,"UniqueDiscordantPairs":23}
    }

          
                
def to_phred(x):
    return -10*(log10(x))

def from_phred(x):
    return pow(10, -x/10)


def test_deletion_insert_size(my_read_pair, all_variants, rg_name):
    insert_dist = inserts[rg_name]
    insert_median = insert_dist["Median"]
    insert_mad = insert_dist["MAD"]

    read_insert = abs(my_read_pair.get_insert_size())
    
    if read_insert > insert_median + 9 * insert_mad:
        #print "{} {}".format(read_insert, insert_median + 9 * insert_mad)
        all_variants.add_read_pair(my_read_pair, DEL_TYPE)
        
def test_insertion_insert_size(my_read_pair, all_variants, rg_name):
    insert_dist = inserts[rg_name]
    insert_median = insert_dist["Median"]
    insert_mad = insert_dist["MAD"]

    read_insert = abs(my_read_pair.get_insert_size())
  
    if read_insert < insert_median - 10 * insert_mad and not my_read_pair.reads_overlap():
        print "{} {}".format(read_insert, insert_median + 5 * insert_mad)
        all_variants.add_read_pair(my_read_pair, INS_TYPE)

def test_split_reads(my_read_pair, all_variants):
    first = my_read_pair.first
    second = my_read_pair.second
    
    SOFT_CLIP = 4
    HARD_CLIP = 5
    
    SOFT_CLIPE_SIZE_THRESHOLD = 20
    
    found = False
    
    for el in first.cigartuples:
        if el[0] == SOFT_CLIP and el[1] >= SOFT_CLIPE_SIZE_THRESHOLD:
            #print first.cigartuples
            found = True
        elif el[0] == HARD_CLIP:
            print "Hard clip"
            
    for el in second.cigartuples:
        if el[0] == SOFT_CLIP and el[1] >= SOFT_CLIPE_SIZE_THRESHOLD:
            #print second.cigartuples
            found = True
        elif el[0] == HARD_CLIP:
            print "Hard clip"
            
    return found


def process_new_read_pair(my_read_pair, all_variants):
    first = my_read_pair.first
    second = my_read_pair.second
    
    first_rg_name = first.get_tag("RG")
    second_rg_name = second.get_tag("RG")
    
    if first_rg_name != second_rg_name:
        print "Bad read group"
        exit(1)
    else:
        test_deletion_insert_size(my_read_pair, all_variants, first_rg_name)
        #test_insertion_insert_size(my_read_pair, all_variants, first_rg_name)
        test_split_reads(my_read_pair, all_variants)
    
#    print my_read_pair

def process_reads(my_reads, waiting_reads, all_variants):
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
            
            my_other_read = waiting_reads.get(read_id)
            
            if my_other_read != None:
                new_read_pair = None             
                if read.reference_start <= my_other_read.reference_start:
                    new_read_pair = ReadPair(read_id, read, my_other_read)
                elif read.reference_start > my_other_read.reference_start:
                    new_read_pair = ReadPair(read_id, my_other_read, read)
                    new_read_pair.was_reversed = True
                else:
                    print "Bad state " + read
                    
                process_new_read_pair(new_read_pair, all_variants)
            else:
                waiting_reads[read_id] = read

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
    shuffle(my_reads)
    waiting_reads = {}
    all_variants = VariantCandidateCollection()
    my_snp = []
    
    process_reads(my_reads, waiting_reads, all_variants)
    print(["type:{} depth:{} key:{} left:{} right:{}\n".format(var_types[all_variants.variants[variant].var_type], 
                                all_variants.variants[variant].depth, 
                                variant, 
                                all_variants.variants[variant].left_pos,
                                all_variants.variants[variant].right_pos) 
            for variant in all_variants.variants])
    print "Full " + str(len(all_variants.variants))
    all_variants.merge()
    print "Post-merge " + str(len(all_variants.variants))
    all_variants.purge()
    print "Post-purge " + str(len(all_variants.variants))
    for key in all_variants.variants:
        my_variant = all_variants.variants[key]
        print my_variant.inner_span_to_str()

def call_test():
    call_variants(my_start = 9999900,
                  my_stop = 10250000,
                  my_chr = "20",
                  input_file = "/Users/siakhnin/data/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam",
                  reference_file = "/Users/siakhnin/data/human_g1k_v37.20.21.fasta",
                  output_file = None)
    
#print timeit.timeit(call_test,number=1)
call_test()
in_vcf = pysam.VariantFile("/Users/siakhnin/data/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.freebayes.vcf", "r")
out_vcf = pysam.VariantFile("/Users/siakhnin/data/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.rheos.vcf", "w", header=in_vcf.header)