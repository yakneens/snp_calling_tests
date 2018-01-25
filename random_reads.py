import pysam
from math import log10
from random import shuffle



def to_phred(x):
    return -10*(log10(x))

def from_phred(x):
    return pow(10, -x/10)

def get_priors():
    heterozygosity = 0.8 * 0.001
    prior_1 = heterozygosity
    prior_0 = pow(heterozygosity, 2) #should this be theta over 2?
    prior_2 = 1 - prior_1 - prior_0
    return (prior_0, prior_1, prior_2)

(prior_0, prior_1, prior_2) = get_priors()

my_start = 10000000
my_stop = 10001000
my_chr = "20"
infile = "/Users/siakhnin/data/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam"
ref =  pysam.FastaFile("/Users/siakhnin/data/human_g1k_v37.20.21.fasta").fetch(region="{}:{}-{}".format(my_chr, my_start + 1, my_stop + 1))
samfile = pysam.AlignmentFile(infile, "rb")
reads = samfile.fetch(region="{}:{}-{}".format(my_chr, my_start + 1, my_stop + 1))
my_reads = [a for a in reads]
shuffle(my_reads, lambda: 0.4)
sites = {}
my_snp = []
#rec - (buf, gl_0, gl_1, gl_2, ref, alt)

def handle_match(read, read_idx, match_len, ref_offset):
#    print(read.query_alignment_sequence[read_idx:read_idx+match_len])
    
    seq = read.query_alignment_sequence[read_idx:read_idx+match_len]
    qual = read.query_alignment_qualities[read_idx:read_idx+match_len]
    ref = read.get_reference_sequence()[ref_offset:ref_offset+match_len]
    
    flag = True
    
    for i in range(len(seq)):
        cur = read.pos + read_idx + i
        
        my_rec = sites.get(cur)
        if not my_rec:
            my_rec = [[], [prior_0, prior_1, prior_2], ref[i], None, cur+1]
        
        p_error = from_phred(qual[i])
        
        my_rec[1][1] = my_rec[1][1] * 0.5
        my_rec[0].append((seq[i], ref[i], p_error, read.qname))
        if(seq[i] != ref[i]):
            if flag:
                print(read)
                flag = False    
            print (seq[i], qual[i], ref[i])
            
            my_rec[3] = seq[i]
            my_rec[1][2] = my_rec[1][2] * p_error
            my_rec[1][0] = my_rec[1][0] * (1 - p_error)
        else:
            my_rec[1][2] = my_rec[1][2] * (1 - p_error)
            my_rec[1][0] = my_rec[1][0] * p_error
        
        sites[cur] = my_rec

others = []    
a = ['20GAVAAXX100126:7:64:9392:84449']
for read in my_reads:
    read_name = read.qname

    if read_name in a:
        print("Found")
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
                    handle_match(read, read_idx, el_len, ref_offset)
                ref_offset = ref_offset + el_len
            if el_type != 2 and el_type != 4:
                read_idx = read_idx + el_len

print_results = True
for pos in sorted(sites):
    rec = sites[pos]
    gls = rec[1]
    ref = rec[2]
    alt = rec[3]
    #pos = rec[4]
    reads = rec[0]
    norm = sum(gls)
    gls = [gl/norm for gl in gls]
    (post_0, post_1, post_2) = gls
    
    vqual = 0
    if post_2 != 0:
        vqual = to_phred(post_2)
    else:
        vqual = 5000
        
    if vqual > 20:   
        res = (pos+1, 0 if post_0 > post_1 else 1, ref, alt, vqual, post_2, post_1, post_0)
        my_snp.append(rec)
        #print pos
        if print_results:
            print res
            
        n_gl_0 = 1
        n_gl_1 = 1
        n_gl_2 = 1
        read_ids = []   
        for read in reads:
            n_gl_1 = n_gl_1 * 0.5
            if(read[0] == read[1]):
                n_gl_0 = n_gl_0 * read[2]
                n_gl_2 = n_gl_2 * (1 - read[2])            
            else:
                n_gl_2 = n_gl_2 * read[2]
                n_gl_0 = n_gl_0 * (1 - read[2])
            
            read_ids.append(read[3])
        
        n_gl_0 = n_gl_0 * prior_0
        n_gl_1 = n_gl_1 * prior_1
        n_gl_2 = n_gl_2 * prior_2
        
        norm = n_gl_0 + n_gl_1 + n_gl_2
             
        n_gl_0 = n_gl_0 / norm
        n_gl_1 = n_gl_1 / norm
        n_gl_2 = n_gl_2 / norm
        
        #print "{} {} {}".format(n_gl_2, n_gl_1, n_gl_0)
        print read_ids
              
#print "Other CIGARs {}.".format(others)
#a = sites[10000438]
#a_sum = a[1][0] + a[1][1] + a[1][2]    
#print "{} {} {} {}".format(a, a[1][0] / a_sum, a[1][1] / a_sum, a[1][2] / a_sum)
