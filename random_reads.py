import pysam
from math import log10
from random import shuffle
import timeit

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


def handle_match(read, read_idx, match_len, ref_offset, all_loci):
    
    seq = read.query_alignment_sequence[read_idx:read_idx+match_len]
    qual = read.query_alignment_qualities[read_idx:read_idx+match_len]
    read_ref = read.get_reference_sequence()[ref_offset:ref_offset+match_len]
    
    for i in range(len(seq)):
        cur = read.pos + read_idx + i
        
        my_locus = all_loci.get(cur)
        if not my_locus:
            my_locus = Locus(pos=cur+1, ref=read_ref[i], alt=None, gl_ref=prior_2, gl_het=prior_1, gl_hom=prior_0, history=[])
        
        p_error = from_phred(qual[i])
        
        my_locus.dp += 1
        
        my_locus.gl_het *= 0.5
        
        if(seq[i] != read_ref[i]):
            my_locus.alt = seq[i]
            my_locus.gl_ref *= p_error
            my_locus.gl_hom *= (1 - p_error)
            my_locus.ao += 1
            my_locus.qa += qual[i]

        else:
            my_locus.gl_ref *= (1 - p_error)
            my_locus.gl_hom *= p_error
            my_locus.ro += 1
            my_locus.qr += qual[i]

        my_locus.add_to_history(my_locus.gl_ref)        
        
        all_loci[cur] = my_locus



def process_reads(my_reads, all_loci):
    for read in my_reads:
        if not read.is_duplicate and not read.is_qcfail and not read.mapping_quality == 0:
            
            cigar_els = read.cigartuples
            
            if not cigar_els:
                print("Read {} has bad CIGAR {}, skipping".format(read.qname, read.cigar))
                continue
                
            ref_offset = 0
            
            read_idx = 0
            for el in cigar_els:
                el_type = el[0]
                el_len = el[1]
                
                if el_type == 0 or el_type == 2:
                    if el_type == 0:
                        handle_match(read, read_idx, el_len, ref_offset, all_loci)
                    ref_offset += el_len
                if el_type != 2 and el_type != 4:
                    read_idx += el_len

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
    all_loci = {}
    my_snp = []
    
    process_reads(my_reads, all_loci)
    
    for cur_pos in sorted(all_loci):
        cur_locus = all_loci[cur_pos]
        (post_ref, post_het, post_hom)  = cur_locus.get_norm_gls()
        
        vqual = 0
        if post_ref != 0:
            vqual = to_phred(post_ref)
        else:
            vqual = 5000
            
        if vqual > 20:   
            res = (cur_pos + 1, 0 if post_hom > post_het else 1, cur_locus.ref, cur_locus.alt, vqual, post_ref, post_het, post_hom)
            my_snp.append(cur_locus)
            probs = cur_locus.history
            if print_results:
                print res
                #print probs

def call_test():
    call_variants(my_start = 10000000,
                  my_stop = 10001000,
                  my_chr = "20",
                  input_file = "/Users/siakhnin/data/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam",
                  reference_file = "/Users/siakhnin/data/human_g1k_v37.20.21.fasta",
                  output_file = None)
    
print timeit.timeit(call_test,number=1)

in_vcf = pysam.VariantFile("/Users/siakhnin/data/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.freebayes.vcf", "r")
out_vcf = pysam.VariantFile("/Users/siakhnin/data/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.rheos.vcf", "w", header=in_vcf.header)