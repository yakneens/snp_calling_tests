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
z = [u'HISEQ1:11:H8GV6ADXX:1:1110:9483:49701',
     u'HWI-D00360:5:H814YADXX:1:2208:1129:93747',
     u'HWI-D00360:5:H814YADXX:1:1213:19452:64030',
     u'HWI-D00360:5:H814YADXX:1:1209:12208:98620',
     u'HISEQ1:12:H8GVUADXX:1:2206:11424:59322',
     u'HWI-D00360:8:H88U0ADXX:2:1208:10699:74467',
     u'HISEQ1:9:H8962ADXX:2:1108:6208:16687',
     u'HISEQ1:9:H8962ADXX:2:2108:9624:93936',
     u'HWI-D00360:7:H88WKADXX:1:2202:6276:57737',
     u'HWI-D00360:7:H88WKADXX:2:1208:19783:4676']

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

            # if name.split("/")[0] in z:
            #     accum.append(hit_dict)
            producer.send('mapped_reads', hit_dict, partition=get_partition(hit_dict), key="{}:{}".format(hit_dict['contig'], hit_dict['pos']))

            counter +=1
            break
        else:
            print "Read maps to {}. Discarding".format(hit.ctg)
print(counter)
print(accum)
# out_file = open("/Users/siakhnin/data/RMNISTHS_30xdownsample_9999999_11000000.mapped.sr.msgpack","w")
# msgpack.dump(hits, out_file, encoding='utf-8')
# out_file.close()