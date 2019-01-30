import redis
from math import log10
import timeit
import pysam
import pickle


def to_phred(x):
    return -10*(log10(x))

def from_phred(x):
    return pow(10, -x/10)

def ascii_to_qual(qual_str):
    return map(lambda x: ord(x) - 33, qual_str)

def get_priors():
    heterozygosity = 0.8 * 0.001
    prior_1 = heterozygosity
    prior_0 = pow(heterozygosity, 2)
    prior_2 = 1 - prior_1 - prior_0
    return (prior_0, prior_1, prior_2)

(prior_0, prior_1, prior_2) = get_priors()

redis_cli = redis.Redis(host='localhost',port=6379)

def call_variants(input_file,
                  reference_file,
                  output_file,
                  my_start,
                  my_stop,
                  my_chr,
                  out_vcf,
                  base_rec,
                  print_results=False, ):
    hit_counter = 0
    my_snp = []
    for cur_pos in range(my_start, my_stop):
        my_locus_redis = redis_cli.get(cur_pos)

        if my_locus_redis:
            cur_locus = pickle.loads(my_locus_redis)
            cur_locus.update_gls_to_min_if_zero()
            (post_ref, post_het, post_hom) = cur_locus.get_norm_gls()

            vqual = 0
            if post_ref != 0:
                vqual = to_phred(post_ref)
            else:
                vqual = 5000
            #       print vqual
            if vqual > 50:
                res = (
                    cur_pos + 1, 0 if post_hom > post_het else 1, cur_locus.ref, cur_locus.alt, vqual, post_ref, post_het,
                post_hom)
                my_snp.append(cur_locus)
                probs = cur_locus.history
                if print_results:
                    print(res)
                    # print probs
                new_rec = base_rec.copy()
                new_rec.info.clear()  # Erase old values
                new_rec.ref = cur_locus.ref
                new_rec.alts = [cur_locus.alt]
                new_rec.qual = vqual
                new_rec.pos = cur_pos + 1
                new_rec.contig = my_chr
                the_sample = new_rec.samples[0]
                if post_hom > post_het:
                    the_sample['GT'] = (1, 1)
                else:
                    the_sample['GT'] = (0, 1)
                the_sample['DP'] = cur_locus.dp
                the_sample['AD'] = (cur_locus.ro, cur_locus.ao)
                #            the_sample['RO'] = cur_locus.ro
                #            the_sample['QR'] = cur_locus.qr
                #            the_sample['AO'] = cur_locus.ao
                #            the_sample['QA'] = cur_locus.qa
                #            the_sample['GL'] = cur_locus.get_log_gls()
                the_sample['PL'] = cur_locus.get_pl()
                the_sample['GQ'] = cur_locus.get_gq()
                out_vcf.write(new_rec)
                hit_counter+=1
        else:
            pass
            print("Locus {} not found".format(cur_pos))
    print("{} variants saved".format(hit_counter))

# 63025520
# 13025520

def call_test():
    in_vcf = pysam.VariantFile("/Users/siakhnin/data/mnist_na12878_chrom20.gatk4.vcf", "r")
    base_rec = None
    for rec in in_vcf.fetch():
        base_rec = rec.copy()
        break
    out_vcf = pysam.VariantFile("/Users/siakhnin/data//mnist_na12878_chrom20_50000000-63025520_test.rheos.kafka.vcf", "w",
                                header=in_vcf.header)
    call_variants(my_start=50000000,
                  my_stop=63025520,
                  my_chr="20",
                  input_file="/Users/siakhnin/data/giab/mnist_na12878_chrom20.bam",
                  reference_file="/Users/siakhnin/data/reference/genome.fa",
                  output_file=None,
                  out_vcf=out_vcf,
                  base_rec=base_rec)


print(timeit.timeit(call_test, number=1))