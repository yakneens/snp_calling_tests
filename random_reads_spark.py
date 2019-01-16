import pysam
from math import log10
from random import shuffle
import timeit
import sys
import dill as pickle

MIN_GL = sys.float_info.min


class Locus:
    def __init__(self, pos, ref, alt, gl_ref, gl_het, gl_hom, history):
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.gl_ref = gl_ref
        self.gl_het = gl_het
        self.gl_hom = gl_hom
        self.history = history
        # Read Depth
        self.dp = 0
        # Reference observation count
        self.ro = 0
        # Sum of reference observation qualities
        self.qr = 0
        self.ao = 0
        # Sum of alternate observation qualities
        self.qa = 0

    def add_to_history(self, new_val):
        self.history.append(new_val)

    def get_gls(self):
        return (self.gl_ref, self.gl_het, self.gl_hom)

    def get_log_gls(self):
        return (log10(self.gl_ref), log10(self.gl_het), log10(self.gl_hom))

    def get_pl(self):
        vals = (int(to_phred(self.gl_ref)), int(to_phred(self.gl_het)), int(to_phred(self.gl_hom)))
        return map(lambda x: x - min(vals), vals)

    def get_gq(self):
        pl = self.get_pl()
        if (pl[1] == 0):
            return pl[2]
        elif (pl[2] == 0):
            return pl[1]
        else:
            print "Something went wrong in PL calculation"
            exit(1)

    def get_norm_gls(self):
        norm_const = self.gl_ref + self.gl_het + self.gl_hom
        return (self.gl_ref / norm_const, self.gl_het / norm_const, self.gl_hom / norm_const)

    def update_gls_to_min_if_zero(self):
        self.gl_hom = min_gl_if_zero(self.gl_hom)
        self.gl_het = min_gl_if_zero(self.gl_het)
        self.gl_ref = min_gl_if_zero(self.gl_ref)


def to_phred(x):
    return -10 * (log10(x))


def from_phred(x):
    return pow(10, -x / 10)


def get_priors():
    heterozygosity = 0.8 * 0.001
    prior_1 = heterozygosity
    prior_0 = pow(heterozygosity, 2)
    prior_2 = 1 - prior_1 - prior_0
    return (prior_0, prior_1, prior_2)


(prior_0, prior_1, prior_2) = get_priors()


def min_gl_if_zero(val):
    if val == 0:
        return MIN_GL
    else:
        return val


def handle_match(read, read_idx, match_len, ref_offset, all_loci, ref, my_start, my_stop):
    seq = read.query_alignment_sequence[read_idx:read_idx + match_len]
    qual = read.query_alignment_qualities[read_idx:read_idx + match_len]
    # read_ref = read.get_reference_sequence()[ref_offset:ref_offset+match_len]
    ref_start = read.reference_start + ref_offset
    ref_stop = read.reference_start + ref_offset + match_len
    if ref_stop < my_start:
        return
    else:
        read_ref = ref[ref_start:ref_stop]
        for i in range(len(seq)):
            cur = read.pos + read_idx + i
            cur_ref = ref_start + i
            if cur_ref < my_start:
                continue
            elif cur_ref > my_stop:
                return
            else:
                my_locus = all_loci.get(cur_ref)
                if not my_locus:
                    my_locus = Locus(pos=cur_ref + 1, ref=ref[cur_ref - my_start], alt=None, gl_ref=prior_2,
                                     gl_het=prior_1, gl_hom=prior_0, history=[])
                p_error = from_phred(qual[i])
                my_locus.dp += 1
                my_locus.gl_het *= 0.5
                if (seq[i] != ref[cur_ref - my_start]):
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
                # my_locus.add_to_history(my_locus.gl_ref)
                my_locus.update_gls_to_min_if_zero()
                all_loci[cur_ref] = my_locus

def process_read(read_str, alignment_header, all_loci, ref, my_start, my_stop):
    read = pysam.AlignedSegment.fromstring(read_str, alignment_header)
    if not read.is_duplicate and not read.is_qcfail and not read.mapping_quality == 0:
        cigar_els = read.cigartuples
        if not cigar_els:
            print("Read {} has bad CIGAR {}, skipping".format(read.qname, read.cigar))
            return
        ref_offset = 0
        read_idx = 0
        for el in cigar_els:
            el_type = el[0]
            el_len = el[1]
            if el_type == 0 or el_type == 2:
                if el_type == 0:
                    handle_match(read, read_idx, el_len, ref_offset, all_loci, ref, my_start, my_stop)
                ref_offset += el_len
            if el_type != 2 and el_type != 4:
                read_idx += el_len


def process_reads(my_reads, all_loci, ref, my_start, my_stop):
    idx = 0
    for read in my_reads:
        idx += 1
        if idx % 10000 == 0:
            print("Processing read #: " + str(idx))
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
                        handle_match(read, read_idx, el_len, ref_offset, all_loci, ref, my_start, my_stop)
                    ref_offset += el_len
                if el_type != 2 and el_type != 4:
                    read_idx += el_len


def get_reads_and_ref(input_file,
              reference_file,
              my_start,
              my_stop,
              my_chr, ):
    print("About to start reading reads")
    ref = pysam.FastaFile(reference_file).fetch(region="{}:{}-{}".format(my_chr, my_start + 1, my_stop + 1))
    samfile = pysam.AlignmentFile(input_file, "rb")
    reads = samfile.fetch(region="{}:{}-{}".format(my_chr, my_start + 1, my_stop + 1))
    my_reads = [a for a in reads]
    # shuffle(my_reads, lambda: 0.4)
    print("About to start shuffling reads")
    shuffle(my_reads)
    return (my_reads, ref)

all_loci = {}
my_start=9999999
my_stop=10001000
my_chr="20"
(my_reads,ref) = get_reads_and_ref(my_start=my_start,
                     my_stop=my_stop,
                     my_chr=my_chr,
                     input_file="/Users/siakhnin/data/giab/mnist_na12878_chrom20.bam",
                     reference_file="/Users/siakhnin/data/reference/genome.fa")

def test_func(read_str, header):
    log = logging.getLogger("py4j")
    log.warn("Got here")
    my_read = pysam.AlignedSegment(header)
    #my_read = pysam.AlignedSegment.fromstring(read_str, header)
    return my_read.tostring()

my_alignment_header = my_reads[0].header

my_reads = [read.tostring() for read in my_reads]

read_rdd = sc.parallelize(my_reads)
new_rdd = read_rdd.map(lambda x: process_read(x, my_alignment_header, all_loci,ref,my_start,my_stop,))

read_rdd.map(lambda x: test_func(x, my_alignment_header)).collect()


def call_variants(input_file,
                  reference_file,
                  output_file,
                  my_start,
                  my_stop,
                  my_chr,
                  out_vcf,
                  base_rec,
                  print_results=True, ):
    print("About to start reading reads")
    ref = pysam.FastaFile(reference_file).fetch(region="{}:{}-{}".format(my_chr, my_start + 1, my_stop + 1))
    samfile = pysam.AlignmentFile(input_file, "rb")
    reads = samfile.fetch(region="{}:{}-{}".format(my_chr, my_start + 1, my_stop + 1))
    my_reads = [a for a in reads]
    # shuffle(my_reads, lambda: 0.4)
    print("About to start shuffling reads")
    shuffle(my_reads)
    all_loci = {}
    my_snp = []
    print("About to start processing reads")
    process_reads(my_reads, all_loci, ref, my_start, my_stop)
    print("Finished processing")
    for cur_pos in sorted(all_loci):
        cur_locus = all_loci[cur_pos]
        (post_ref, post_het, post_hom) = cur_locus.get_norm_gls()

        vqual = 0
        if post_ref != 0:
            vqual = to_phred(post_ref)
        else:
            vqual = 5000
        #       print vqual
        if vqual > 20:
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


# 63025520
# 13025520

def call_test():
    in_vcf = pysam.VariantFile("/Users/siakhnin/data/mnist_na12878_chrom20.gatk4.vcf", "r")
    base_rec = None
    for rec in in_vcf.fetch():
        base_rec = rec.copy()
        break
    out_vcf = pysam.VariantFile("/Users/siakhnin/data//mnist_na12878_chrom20_10mb.rheos.vcf", "w",
                                header=in_vcf.header)
    call_variants(my_start=9999999,
                  my_stop=20000000,
                  my_chr="20",
                  input_file="/Users/siakhnin/data/giab/mnist_na12878_chrom20.bam",
                  reference_file="/Users/siakhnin/data/reference/genome.fa",
                  output_file=None,
                  out_vcf=out_vcf,
                  base_rec=base_rec)


print(timeit.timeit(call_test, number=1))
