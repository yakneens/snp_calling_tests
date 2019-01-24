import mappy as mp
from kafka import KafkaProducer
from json import dumps
import msgpack

my_start = 9999999
my_stop = 20000000
my_chr = "20"
reference_file = "/Users/siakhnin/data/reference/genome.mmi"
sample_file = "/Users/siakhnin/data/giab/RMNISTHS_30xdownsample_9999999_10100000.name_sorted.fastq"
print("Loading reference")
ref = mp.Aligner(reference_file, preset="sr")
print("Reference Loaded")

def get_partition(rec):
    if rec['pos'] > 30000000:
        return 1
    else:
        return 0

producer = KafkaProducer(bootstrap_servers=['localhost:9092'],
                         value_serializer=lambda x: dumps(x).encode('utf-8'))


accum = []

hits = []

counter = 0
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

            #print("Hit number " + str(hit_count))
            #print("{}\t{}\t{}\t{}".format(hit.ctg, hit.r_st, hit.r_en, hit.cigar_str))
            hit_count+=1
            hits.append(hit_dict)
            counter += 1

            if(len(hits) % 1500 == 0):
                producer.send('mapped_reads_batchy', hits, partition=get_partition(hit_dict))
                hits = []
                print("{} reads processed".format(counter))

            break
        else:
            pass
            #print "Read maps to {}. Discarding".format(hit.ctg)

if(len(hits) > 0):
    producer.send('mapped_reads_batchy', hits, partition=get_partition(hit_dict))
    hits = []
    print("{} reads processed".format(counter))

print(counter)
print(accum)
# out_file = open("/Users/siakhnin/data/RMNISTHS_30xdownsample_9999999_11000000.mapped.sr.msgpack","w")
# msgpack.dump(hits, out_file, encoding='utf-8')
# out_file.close()