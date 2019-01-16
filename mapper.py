import mappy as mp
from kafka import KafkaProducer
from json import dumps

my_start = 9999999
my_stop = 20000000
my_chr = "20"
reference_file = "/Users/siakhnin/data/reference/genome.mmi"
sample_file = "/Users/siakhnin/data/giab/RMNISTHS_30xdownsample_9999999_10100000.fastq"
print("Loading reference")
ref = mp.Aligner(reference_file)
print("Reference Loaded")

def get_partition(rec):
    if rec['pos'] > 30000000:
        return 1
    else:
        return 0

producer = KafkaProducer(bootstrap_servers=['localhost:9092'],
                         value_serializer=lambda x: dumps(x).encode('utf-8'))


for name, seq, qual in mp.fastx_read(sample_file):
    hit_count = 1
    for hit in ref.map(seq):
        if hit.ctg == "20":
            hit_dict = {}
            hit_dict['contig'] = hit.ctg
            hit_dict['cigar'] = hit.cigar_str
            hit_dict['seq'] = seq
            hit_dict['qname'] = name
            hit_dict['qual'] = qual
            hit_dict['strand'] = hit.strand
            hit_dict['cigartuples'] = hit.cigar
            hit_dict['pos'] = hit.r_st
            hit_dict['end'] = hit.r_en
            hit_dict['qstart'] = hit.q_st
            hit_dict['qend'] = hit.q_en
            hit_dict['mapq'] = hit.mapq

            print("Hit number " + str(hit_count))
            print("{}\t{}\t{}\t{}".format(hit.ctg, hit.r_st, hit.r_en, hit.cigar_str))
            hit_count+=1

            producer.send('mapped_reads', hit_dict, partition=get_partition(hit_dict), key="{}:{}".format(hit_dict['contig'], hit_dict['pos']))
            break
        else:
            print "Read maps to {}. Discarding".format(hit.ctg)
