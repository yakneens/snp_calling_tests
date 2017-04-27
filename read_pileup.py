import pysam

my_snp = []
samfile = pysam.AlignmentFile("/Users/siakhnin/data/snp_calling_tests/inputs/CEUTrio.HiSeq.WGS.b37.NA12878.20.10000758.10.bam", "rb")
my_iter = samfile.pileup(reference="20", start=10000757, end=10000758)
for el in my_iter:
    #This is 0-based so is off by 1
    if el.reference_pos == 10000757:
        for cur_read in el.pileups:
            aln = cur_read.alignment
            loc = cur_read.query_position
            my_snp.append((aln.query_alignment_sequence[loc], 
                           aln.query_qualities[loc]))
        
        
print my_snp
ref =  pysam.FastaFile("/Users/siakhnin/data/human_g1k_v37.20.21.fasta").fetch(region="20:10000758-10000758")
print "Ref is : " + str(ref)
