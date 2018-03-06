import pysam
from math import log10
from random import shuffle
import timeit

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

class Locus:
    def __init__(self, pos, ref, alt, gl_ref, gl_het, gl_hom, history):
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.gl_ref = gl_ref
        self.gl_het = gl_het
        self.gl_hom = gl_hom
        self.history = history
        #Read Depth
        self.dp = 0
        #Reference observation count
        self.ro = 0
        #Sum of reference observation qualities
        self.qr = 0
        #Alternate observation count
        self.ao = 0
        #Sum of alternate observation qualities
        self.qa = 0
        
    def add_to_history(self, new_val):
        self.history.append(new_val)

    def get_gls(self):
        return (self.gl_ref, self.gl_het, self.gl_hom)
    
    def get_norm_gls(self):
        norm_const = self.gl_ref + self.gl_het + self.gl_hom
        return (self.gl_ref / norm_const, self.gl_het / norm_const, self.gl_hom / norm_const)
  
class ReadPair:
    def __init__(self, read_id, first, second):
        self.first = first
        self.second = second
        self.read_id = read_id
        
    def __str__(self):
        return "{} {} {}".format(self.read_id, self.first, self.second)
    
def to_phred(x):
    return -10*(log10(x))

def from_phred(x):
    return pow(10, -x/10)

def get_priors():
    heterozygosity = 0.8 * 0.001
    prior_1 = heterozygosity
    prior_0 = pow(heterozygosity, 2)
    prior_2 = 1 - prior_1 - prior_0
    return (prior_0, prior_1, prior_2)

(prior_0, prior_1, prior_2) = get_priors()

def process_new_read_pair(my_read_pair, all_variants):
    print my_read_pair

def process_reads(my_reads, waiting_read_pairs, all_variants):
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
            read_is_first = True if read.is_read1 else False
            
            my_read_pair = waiting_read_pairs.get(read_id)
            
            if my_read_pair != None:
                if read_is_first and my_read_pair.first == None:
                    my_read_pair.first = read
                elif not read_is_first and my_read_pair.second == None:
                    my_read_pair.second = read
                else:
                    print "Bad state " + read
                    
                process_new_read_pair(my_read_pair, all_variants)
            else:
                new_read_pair = None
                
                if read_is_first:
                    new_read_pair = ReadPair(read_id, read, None)                
                else:
                    new_read_pair = ReadPair(read_id, None, read)
            
                waiting_read_pairs[read_id] = new_read_pair

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
    waiting_read_pairs = {}
    all_variants = {}
    my_snp = []
    
    process_reads(my_reads, waiting_read_pairs, all_variants)
    print "here"
#     for cur_pos in sorted(all_loci):
#         cur_locus = all_loci[cur_pos]
#         (post_ref, post_het, post_hom)  = cur_locus.get_norm_gls()
#         
#         vqual = 0
#         if post_ref != 0:
#             vqual = to_phred(post_ref)
#         else:
#             vqual = 5000
#             
#         if vqual > 20:   
#             res = (cur_pos + 1, 0 if post_hom > post_het else 1, cur_locus.ref, cur_locus.alt, vqual, post_ref, post_het, post_hom)
#             my_snp.append(cur_locus)
#             probs = cur_locus.history
#             if print_results:
#                 print res
                #print probs

def call_test():
    call_variants(my_start = 10000000,
                  my_stop = 10001000,
                  my_chr = "20",
                  input_file = "/Users/siakhnin/data/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam",
                  reference_file = "/Users/siakhnin/data/human_g1k_v37.20.21.fasta",
                  output_file = None)
    
#print timeit.timeit(call_test,number=1)
call_test()
in_vcf = pysam.VariantFile("/Users/siakhnin/data/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.freebayes.vcf", "r")
out_vcf = pysam.VariantFile("/Users/siakhnin/data/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.rheos.vcf", "w", header=in_vcf.header)