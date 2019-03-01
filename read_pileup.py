import pysam
from math import pow, log10
import timeit

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
#print "Priors: ref - {} het - {} hom-var - {}".format(prior_2, prior_1, prior_0)

def get_pileup(filename, my_start, my_stop, my_chr):
    samfile = pysam.AlignmentFile(filename, "rb")
    my_iter = samfile.pileup(contig=my_chr, start=my_start, end=my_stop, truncate=True)
    return my_iter




def call_variants(input_file, 
                  reference_file, 
                  output_file,
                  my_start,
                  my_stop,
                  my_chr, 
                  use_log=False, 
                  print_results=True):

    my_snp = []
    
    #This is 0-based so is off by 1
    my_iter = get_pileup(input_file, my_start, my_stop, my_chr)
    
    ref =  pysam.FastaFile(reference_file).fetch(region="{}:{}-{}".format(my_chr, my_start + 1, my_stop + 1))
    
    f = open(output_file, 'w')
    
    
    for el in my_iter:
        #This is 0-based so is off by 1
        pos = el.reference_pos + 1
        
        if print_results and pos % 1000000 == 0:
            print(pos)
        
        ref_ind = el.reference_pos - my_start
        
        
        #print pos
        gl_0 = 1
        gl_1 = 1
        gl_2 = 1
        
        if use_log:           
            lgl_0 = 0.0
            lgl_1 = 0.0
            lgl_2 = 0.0
        
        alt = None
        read_ids = []
        for cur_read in el.pileups:
            if not cur_read.is_del:
                aln = cur_read.alignment
                
                if not aln.is_duplicate and not aln.is_qcfail and aln.mapping_quality > 0:
                    read_ids.append(aln.qname)
                    loc = cur_read.query_position
                    #print loc
                    base = aln.seq[loc]
                    qual = aln.query_qualities[loc]
                    
                    p_error = from_phred(qual)
              
                    gl_1 = gl_1 * 0.5
                    
                    if use_log: 
                        log_10_of_half = -0.301029995663981
                        lgl_1 = lgl_1 + log_10_of_half
                           
                    if(base == ref[ref_ind]):
                        gl_0 = gl_0 * p_error
                        gl_2 = gl_2 * (1 - p_error)
                        
                        if use_log:
                            lgl_0 = lgl_0 + (-qual / 10)
                            lgl_2 = lgl_2 + (-to_phred(1-p_error)/10)
                    else:
                        alt = base
                        gl_2 = gl_2 * p_error
                        gl_0 = gl_0 * (1 - p_error)
                        
                        if use_log:
                            lgl_2 = lgl_2 + (-qual / 10)
                            lgl_0 = lgl_0 + (-to_phred(1-p_error)/10)
                    
                    #my_snp.append((aln.seq[loc], aln.query_qualities[loc]))
            
            
        norm_const = prior_0 * gl_0 + prior_1 * gl_1 + prior_2 * gl_2
        
        post_0 = prior_0 * gl_0 / norm_const
        post_1 = prior_1 * gl_1 / norm_const
        post_2 = prior_2 * gl_2 / norm_const
        
        if use_log:         
            log_prior_0 = log10(prior_0)
            log_prior_1 = log10(prior_1)
            log_prior_2 = log10(prior_2)
            
            #log_norm_const =  (log_prior_0 + lgl_0) * (log_prior_1 + lgl_1) * (log_prior_2 + lgl_2)       
            
            log_post_0 = log_prior_0 + lgl_0 
            log_post_1 = log_prior_1 + lgl_1
            log_post_2 = log_prior_2 + lgl_2
            
            de_log_post_0 = pow(10, log_post_0)
            de_log_post_1 = pow(10, log_post_1)
            de_log_post_2 = pow(10, log_post_2)
            
            de_log_norm_const = de_log_post_0 + de_log_post_1 + de_log_post_2
            
            norm_de_log_post0 = de_log_post_0 / de_log_norm_const
            norm_de_log_post1 = de_log_post_1 / de_log_norm_const
            norm_de_log_post2 = de_log_post_2 / de_log_norm_const
        
            post_0 = norm_de_log_post0
            post_1 = norm_de_log_post1
            post_2 = norm_de_log_post2
        
        vqual = 0
        if post_2 != 0:
            vqual = to_phred(post_2)
        else:
            vqual = 5000
            
        if vqual > 20:   
            rec = (pos, 0 if post_0 > post_1 else 1, ref[ref_ind], alt, vqual, post_2, post_1, post_0)
            my_snp.append(rec)
            #print pos
            if print_results:
                print(rec)
                #print read_ids
                
            f.write(str(rec))
            

            
    f.close()
    print(len(my_snp))
    return my_snp

#my_start = 1
#my_stop = 10140159
my_start = 10000000
my_stop = 10001000
my_chr = "20"

# call_variants("/Users/siakhnin/data/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam", 
#               "/Users/siakhnin/data/human_g1k_v37.20.21.fasta", 
#               '/Users/siakhnin/data/blah.txt', 
#               my_start,
#               my_stop,
#               my_chr, 
#               use_log=False,
#               print_results=True)

def call_with_log():
    print(call_variants("/Users/siakhnin/data/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam", 
              "/Users/siakhnin/data/human_g1k_v37.20.21.fasta", 
              '/Users/siakhnin/data/blah.txt',
              my_start,
              my_stop,
              my_chr, 
              use_log=True,
              print_results=True))

def call_no_log():
    print(call_variants("/Users/siakhnin/data/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam", 
              "/Users/siakhnin/data/human_g1k_v37.20.21.fasta", 
              '/Users/siakhnin/data/blah.txt', 
              my_start,
              my_stop,
              my_chr, 
              use_log=False,
              print_results=True))
    
print(timeit.timeit(call_no_log, number=1))
print(timeit.timeit(call_with_log, number=1))
