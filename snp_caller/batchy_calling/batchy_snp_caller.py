import pysam
from math import log10
from kafka import KafkaConsumer, KafkaProducer

from model.locus import Locus, q_chars_to_p_errors, q_chars_to_qual

import cProfile

import time
import datetime

import pickle
import zlib
import ujson


SKIP_ALL = False
KAFKA_SERVER="localhost:9092"
all_loci = {}
#rev_comp_lookup = string.maketrans(u'ACTG', u'TGAC')
rev_comp_lookup = {ord(c): ord(t) for c,t in zip(u'ACTG', u'TGAC')}

def comp(seq):
    return seq.translate(rev_comp_lookup)

def rev_comp(seq):
    return comp(seq)[::-1]

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

def min_gl_if_zero(val):
    if val == 0:
        return Locus.MIN_GL
    else:
        return val

def handle_match(read, read_idx, match_len, ref_offset, ref, my_start, my_stop):

    strand = read['strand']

    seq = read['seq'][read_idx:read_idx + match_len]
    qual = read['qual'][read_idx:read_idx + match_len]

    updated_loci = {}

    if strand == 1:
        ref_start = read['pos'] + ref_offset
    else:
        seq = rev_comp(seq)
        qual = read['qual'][read_idx:read_idx+match_len] #TODO Should this be reversed?
        ref_start = read['end'] + ref_offset - match_len

    ref_stop = ref_start + match_len

    if ref_stop < my_start:
        return
    else:
        read_ref = ref[ref_start:ref_stop]

        for i in range(len(seq)):
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

                p_error = q_chars_to_p_errors[qual[i]]
                my_locus.dp += 1

                my_locus.gl_het *= 0.5

                if (seq[i] != ref[cur_ref - my_start]):
                    my_locus.alt = seq[i]
                    my_locus.gl_ref *= p_error

                    my_locus.gl_hom *= (1 - p_error)

                    my_locus.ao += 1
                    my_locus.qa += q_chars_to_qual[qual[i]]

                else:
                    my_locus.gl_ref *= (1 - p_error)
                    my_locus.gl_hom *= p_error
                    my_locus.ro += 1
                    my_locus.qr += q_chars_to_qual[qual[i]]

                # my_locus.add_to_history(my_locus.gl_ref)

                #my_locus.update_gls_to_min_if_zero()
                #TODO Maybe ok to take this out?


                all_loci[cur_ref] = my_locus
                updated_loci[cur_ref] = my_locus

    return updated_loci

def deserialize_read_msg(m):
    return ujson.loads(m.decode('utf-8'))

def get_reads_source():
    consumer = KafkaConsumer('mapped_reads_batchy',
                                group_id='snp_caller',
                                bootstrap_servers=[KAFKA_SERVER],
                                value_deserializer=deserialize_read_msg)
    return consumer

def serialize_locus_msg(m):
    return zlib.compress(pickle.dumps(m, protocol=2))

def get_locus_saver():
    producer = KafkaProducer(bootstrap_servers=[KAFKA_SERVER],
                             value_serializer=serialize_locus_msg,
                             max_request_size=100000000)
    return producer

def process_reads(read_batch, ref, my_start, my_stop):
    updated_loci = {}

    for read in read_batch:
        if SKIP_ALL:
            continue

        strand = read['strand']
        cigar_els = read['cigartuples'] if strand == 1 else read['cigartuples'][::-1]

        if not cigar_els:
            print(f"Read {read.qname} has bad CIGAR {read.cigar}, skipping")
            return

        ref_offset = 0

        read_idx = read['qstart']

        for el in cigar_els:
            el_type = el[1]
            el_len = el[0]

            #0 is a match, 2 is a deletion, these consume reference sequence
            if el_type == 0 or el_type == 2:
                if el_type == 0:
                    updated_loci.update(handle_match(read, read_idx, el_len, ref_offset, ref, my_start, my_stop))
                ref_offset += el_len * strand
            if el_type != 2 and el_type != 4:
                read_idx += el_len
    return updated_loci

def reverse_neg_strand_read(read):
    return True

def process_messages():
    print("Starting Message Processing")

    read_source = get_reads_source()
    locus_saver = get_locus_saver()



    reference_file = "/Users/siakhnin/data/reference/genome.fa"
    ref = pysam.FastaFile(reference_file).fetch(region="20")
    counter = 1

    for message in read_source:
        my_read_batch = message.value

        new_locus_updates = process_reads(my_read_batch, ref, 0, len(ref))

        ts = time.time()
        locus_saver.send('snp_loci_batchy', new_locus_updates, partition=0,)
        te = time.time()
        print(f"{datetime.datetime.now()} - Took {te - ts} secs to send updated loci for {len(my_read_batch)} reads.")
        print(f"Processed {counter} messages. Updating shared dictionary")
        counter += 1


def main():
    cProfile.runctx("process_messages()", globals(), locals(),
                    'profile-batchy-snp-caller.out')

if __name__ == "__main__":
    main()
