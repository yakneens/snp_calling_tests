import pysam
from math import log10
from random import shuffle
import timeit
import math

from variant_candidate import VariantCandidate
from variant_candidate_collection import VariantCandidateCollection
from model.read_pair import ReadPair

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
        "U1b_A5_L1":{"ReadSize":148,"Median":558,"MAD":89,"UniqueDiscordantPairs":41},
        "U0a_A_L2":{"ReadSize":148,"Median":561,"MAD":101,"UniqueDiscordantPairs":52},
        "U5a_A2_L1":{"ReadSize":148,"Median":557,"MAD":113,"UniqueDiscordantPairs":71},
        "U0a_A2_L1":{"ReadSize":148,"Median":564,"MAD":104,"UniqueDiscordantPairs":42},
        "U0c_B_L1":{"ReadSize":148,"Median":563,"MAD":99,"UniqueDiscordantPairs":50},
        "U5a_A_L2":{"ReadSize":148,"Median":552,"MAD":112,"UniqueDiscordantPairs":58},
        "U5b_A_L2":{"ReadSize":148,"Median":498,"MAD":97,"UniqueDiscordantPairs":40},
        "U0a_A2_L2":{"ReadSize":148,"Median":566,"MAD":102,"UniqueDiscordantPairs":55},
        "U0a_A2_L2":{"ReadSize":148,"Median":566,"MAD":102,"UniqueDiscordantPairs":55},
        "U4b_B4_L1":{"ReadSize":148,"Median":575,"MAD":99,"UniqueDiscordantPairs":39},
        "U2a_B4_L1":{"ReadSize":148,"Median":561,"MAD":92,"UniqueDiscordantPairs":39},
        "U0a_A_L1":{"ReadSize":148,"Median":563,"MAD":102,"UniqueDiscordantPairs":48},
        "U4b_B4_L2":{"ReadSize":148,"Median":579,"MAD":99,"UniqueDiscordantPairs":34},
        "U0b_B_L2":{"ReadSize":148,"Median":531,"MAD":118,"UniqueDiscordantPairs":35},
        "U0a_B_L1":{"ReadSize":148,"Median":562,"MAD":102,"UniqueDiscordantPairs":34},
        "U0a_B2_L1":{"ReadSize":148,"Median":563,"MAD":102,"UniqueDiscordantPairs":42},
        "U0c_A2_L1":{"ReadSize":148,"Median":563,"MAD":100,"UniqueDiscordantPairs":48},
        "U5c_B2_L2":{"ReadSize":148,"Median":550,"MAD":92,"UniqueDiscordantPairs":61},
        "U2b_B3_L2":{"ReadSize":148,"Median":552,"MAD":90,"UniqueDiscordantPairs":38},
        "U1a_A5_L1":{"ReadSize":148,"Median":573,"MAD":93,"UniqueDiscordantPairs":37},
        "U0a_B2_L2":{"ReadSize":148,"Median":563,"MAD":103,"UniqueDiscordantPairs":38},
        "U1a_A4_L1":{"ReadSize":148,"Median":563,"MAD":90,"UniqueDiscordantPairs":38},
        "U5c_A2_L2":{"ReadSize":148,"Median":549,"MAD":92,"UniqueDiscordantPairs":61},
        "U0b_A_L1":{"ReadSize":148,"Median":530,"MAD":117,"UniqueDiscordantPairs":37},
        "U0a_B_L2":{"ReadSize":148,"Median":561,"MAD":102,"UniqueDiscordantPairs":39},
        "U2b_A4_L1":{"ReadSize":148,"Median":555,"MAD":89,"UniqueDiscordantPairs":40},
        "U0b_A2_L1":{"ReadSize":148,"Median":531,"MAD":118,"UniqueDiscordantPairs":53},
        "U2b_B3_L1":{"ReadSize":148,"Median":552,"MAD":89,"UniqueDiscordantPairs":40},
        "U0b_A2_L2":{"ReadSize":148,"Median":531,"MAD":118,"UniqueDiscordantPairs":43},
        "U5b_A_L1":{"ReadSize":148,"Median":500,"MAD":96,"UniqueDiscordantPairs":50},
        "U0c_A_L1":{"ReadSize":148,"Median":562,"MAD":99,"UniqueDiscordantPairs":48},
        "U0b_A_L2":{"ReadSize":148,"Median":529,"MAD":117,"UniqueDiscordantPairs":65},
        "U5b_B_L1":{"ReadSize":148,"Median":501,"MAD":97,"UniqueDiscordantPairs":44},
        "U3a_B3_L2":{"ReadSize":148,"Median":554,"MAD":90,"UniqueDiscordantPairs":39},
        "U0b_B_L1":{"ReadSize":148,"Median":529,"MAD":118,"UniqueDiscordantPairs":50},
        "U3b_A4_L2":{"ReadSize":148,"Median":552,"MAD":89,"UniqueDiscordantPairs":35},
        "U0b_B2_L1":{"ReadSize":148,"Median":531,"MAD":117,"UniqueDiscordantPairs":46},
        "U5b_A2_L2":{"ReadSize":148,"Median":500,"MAD":97,"UniqueDiscordantPairs":40},
        "U0b_B2_L2":{"ReadSize":148,"Median":531,"MAD":118,"UniqueDiscordantPairs":44},
        "U0c_B_L2":{"ReadSize":148,"Median":562,"MAD":99,"UniqueDiscordantPairs":36},
        "U1a_A4_L2":{"ReadSize":148,"Median":561,"MAD":91,"UniqueDiscordantPairs":44},
        "U2b_A5_L1":{"ReadSize":148,"Median":563,"MAD":91,"UniqueDiscordantPairs":30},
        "U0c_A2_L2":{"ReadSize":148,"Median":562,"MAD":99,"UniqueDiscordantPairs":50},
        "U3b_B4_L1":{"ReadSize":148,"Median":554,"MAD":90,"UniqueDiscordantPairs":45},
        "U4b_A4_L1":{"ReadSize":148,"Median":576,"MAD":98,"UniqueDiscordantPairs":26},
        "U1a_B3_L2":{"ReadSize":148,"Median":562,"MAD":90,"UniqueDiscordantPairs":46},
        "U1a_B4_L2":{"ReadSize":148,"Median":565,"MAD":92,"UniqueDiscordantPairs":28},
        "U3b_A5_L1":{"ReadSize":148,"Median":564,"MAD":92,"UniqueDiscordantPairs":43},
        "U2a_A5_L2":{"ReadSize":148,"Median":569,"MAD":92,"UniqueDiscordantPairs":37},
        "U4a_B3_L1":{"ReadSize":148,"Median":543,"MAD":90,"UniqueDiscordantPairs":37},
        "U0c_B2_L1":{"ReadSize":148,"Median":561,"MAD":100,"UniqueDiscordantPairs":41},
        "U0c_B2_L2":{"ReadSize":148,"Median":562,"MAD":100,"UniqueDiscordantPairs":47},
        "U1a_B3_L1":{"ReadSize":148,"Median":565,"MAD":91,"UniqueDiscordantPairs":47},
        "U1a_B4_L1":{"ReadSize":148,"Median":563,"MAD":91,"UniqueDiscordantPairs":45},
        "U1b_A4_L1":{"ReadSize":148,"Median":548,"MAD":88,"UniqueDiscordantPairs":43},
        "U3a_A5_L1":{"ReadSize":148,"Median":565,"MAD":92,"UniqueDiscordantPairs":30},
        "U1b_A4_L2":{"ReadSize":148,"Median":547,"MAD":88,"UniqueDiscordantPairs":33},
        "U2a_B3_L1":{"ReadSize":148,"Median":561,"MAD":90,"UniqueDiscordantPairs":50},
        "U5c_B_L1":{"ReadSize":148,"Median":549,"MAD":91,"UniqueDiscordantPairs":51},
        "U1b_A5_L2":{"ReadSize":148,"Median":556,"MAD":89,"UniqueDiscordantPairs":49},
        "U3a_A4_L1":{"ReadSize":148,"Median":556,"MAD":89,"UniqueDiscordantPairs":40},
        "U1b_B3_L1":{"ReadSize":148,"Median":548,"MAD":88,"UniqueDiscordantPairs":39},
        "U2a_B4_L2":{"ReadSize":148,"Median":561,"MAD":92,"UniqueDiscordantPairs":40},
        "U1b_B3_L2":{"ReadSize":148,"Median":549,"MAD":88,"UniqueDiscordantPairs":45},
        "U3a_B4_L1":{"ReadSize":148,"Median":555,"MAD":90,"UniqueDiscordantPairs":42},
        "U2b_B4_L2":{"ReadSize":148,"Median":554,"MAD":90,"UniqueDiscordantPairs":38},
        "U1b_B4_L1":{"ReadSize":148,"Median":548,"MAD":88,"UniqueDiscordantPairs":47},
        "U2b_B4_L1":{"ReadSize":148,"Median":553,"MAD":90,"UniqueDiscordantPairs":43},
        "U3a_A4_L2":{"ReadSize":148,"Median":554,"MAD":90,"UniqueDiscordantPairs":40},
        "U2a_A4_L2":{"ReadSize":148,"Median":559,"MAD":91,"UniqueDiscordantPairs":29},
        "U5c_B_L2":{"ReadSize":148,"Median":550,"MAD":91,"UniqueDiscordantPairs":68},
        "U4b_B3_L1":{"ReadSize":148,"Median":576,"MAD":99,"UniqueDiscordantPairs":52},
        "U2a_A5_L1":{"ReadSize":148,"Median":571,"MAD":91,"UniqueDiscordantPairs":40},
        "U2b_A4_L2":{"ReadSize":148,"Median":552,"MAD":89,"UniqueDiscordantPairs":50},
        "U2b_A5_L2":{"ReadSize":148,"Median":563,"MAD":92,"UniqueDiscordantPairs":46},
        "U3a_A5_L2":{"ReadSize":148,"Median":563,"MAD":91,"UniqueDiscordantPairs":33},
        "U3a_B3_L1":{"ReadSize":148,"Median":555,"MAD":91,"UniqueDiscordantPairs":38},
        "U3a_B4_L2":{"ReadSize":148,"Median":555,"MAD":91,"UniqueDiscordantPairs":38},
        "U3b_A5_L2":{"ReadSize":148,"Median":562,"MAD":90,"UniqueDiscordantPairs":39},
        "U1a_A5_L2":{"ReadSize":148,"Median":572,"MAD":93,"UniqueDiscordantPairs":45},
        "U3b_A4_L1":{"ReadSize":148,"Median":554,"MAD":90,"UniqueDiscordantPairs":46},
        "U1b_B4_L2":{"ReadSize":148,"Median":550,"MAD":88,"UniqueDiscordantPairs":41},
        "U3b_B3_L1":{"ReadSize":148,"Median":553,"MAD":90,"UniqueDiscordantPairs":37},
        "U0c_A_L2":{"ReadSize":148,"Median":561,"MAD":99,"UniqueDiscordantPairs":45},
        "U4a_A4_L2":{"ReadSize":148,"Median":543,"MAD":89,"UniqueDiscordantPairs":54},
        "U4a_A5_L2":{"ReadSize":148,"Median":552,"MAD":90,"UniqueDiscordantPairs":28},
        "U5c_A_L1":{"ReadSize":148,"Median":548,"MAD":91,"UniqueDiscordantPairs":64},
        "U4a_B3_L2":{"ReadSize":148,"Median":543,"MAD":89,"UniqueDiscordantPairs":45},
        "U4a_B4_L1":{"ReadSize":148,"Median":545,"MAD":89,"UniqueDiscordantPairs":49},
        "U4b_A5_L1":{"ReadSize":148,"Median":589,"MAD":101,"UniqueDiscordantPairs":39},
        "U5c_A_L2":{"ReadSize":148,"Median":548,"MAD":91,"UniqueDiscordantPairs":53},
        "U4a_B4_L2":{"ReadSize":148,"Median":545,"MAD":90,"UniqueDiscordantPairs":43},
        "U2a_B3_L2":{"ReadSize":148,"Median":560,"MAD":90,"UniqueDiscordantPairs":44},
        "U5a_A2_L2":{"ReadSize":148,"Median":557,"MAD":114,"UniqueDiscordantPairs":55},
        "U4b_A5_L2":{"ReadSize":148,"Median":587,"MAD":100,"UniqueDiscordantPairs":46},
        "U4b_B3_L2":{"ReadSize":148,"Median":577,"MAD":100,"UniqueDiscordantPairs":58},
        "U4a_A4_L1":{"ReadSize":148,"Median":544,"MAD":89,"UniqueDiscordantPairs":22},
        "U4a_A5_L1":{"ReadSize":148,"Median":554,"MAD":90,"UniqueDiscordantPairs":27},
        "U5b_B_L2":{"ReadSize":148,"Median":501,"MAD":97,"UniqueDiscordantPairs":41},
        "U5a_A_L1":{"ReadSize":148,"Median":554,"MAD":112,"UniqueDiscordantPairs":54},
        "U2a_A4_L1":{"ReadSize":148,"Median":562,"MAD":90,"UniqueDiscordantPairs":51},
        "U4b_A4_L2":{"ReadSize":148,"Median":575,"MAD":98,"UniqueDiscordantPairs":38},
        "U5a_B2_L1":{"ReadSize":148,"Median":553,"MAD":112,"UniqueDiscordantPairs":58},
        "U5c_A2_L1":{"ReadSize":148,"Median":549,"MAD":92,"UniqueDiscordantPairs":50},
        "U5a_B2_L2":{"ReadSize":148,"Median":556,"MAD":113,"UniqueDiscordantPairs":65},
        "U5a_B_L1":{"ReadSize":148,"Median":554,"MAD":113,"UniqueDiscordantPairs":56},
        "U5b_A2_L1":{"ReadSize":148,"Median":501,"MAD":99,"UniqueDiscordantPairs":50},
        "U3b_B4_L2":{"ReadSize":148,"Median":555,"MAD":89,"UniqueDiscordantPairs":37},
        "U5b_B2_L2":{"ReadSize":148,"Median":501,"MAD":98,"UniqueDiscordantPairs":41},
        "U5c_B2_L1":{"ReadSize":148,"Median":549,"MAD":91,"UniqueDiscordantPairs":46},
        "U3b_B3_L2":{"ReadSize":148,"Median":554,"MAD":89,"UniqueDiscordantPairs":41},
        "U5a_B_L2":{"ReadSize":148,"Median":555,"MAD":113,"UniqueDiscordantPairs":48},
        "U5b_B2_L1":{"ReadSize":148,"Median":498,"MAD":97,"UniqueDiscordantPairs":37}
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
    call_variants(my_start = 0,
                  my_stop = 15000000,
                  my_chr = "20",
                  input_file = "/Users/siakhnin/data/giab/mnist_na12878_chrom20.bam",
                  reference_file = "/Users/siakhnin/data/reference/genome.fa",
                  output_file = None)
    
#print timeit.timeit(call_test,number=1)
call_test()
#in_vcf = pysam.VariantFile("/Users/siakhnin/data/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.freebayes.vcf", "r")
#out_vcf = pysam.VariantFile("/Users/siakhnin/data/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.rheos.vcf", "w", header=in_vcf.header)