"""Batchy Mapper

Usage:
    batchy_mapper.py file -i <input> [-g <chr>:<start>-<stop>] -r <ref> -k <kafka_url> [--bin_size=<bin_size>] [--num_bins=<num_bins>] [--check_freq=<check_freq>]
    batchy_mapper.py service -r <ref> -k <kafka_url> [--bin_size=<bin_size>] [--num_bins=<num_bins>] [--check_freq=<check_freq>]
    batchy_mapper.py (-h | --help)

Options:
    -h, --help  show help
    -i <input> path to input file
    -g <region> genome region using <chr>:<start>-<stop> format
    -r <ref> path to reference file
    -k <kafka_url> kafka broker URL
    --bin_size=<bin_size>  how many items fit in a bin before it is sent
    --num_bins=<num_bins>  number of bins
    --check_freq=<check_freq>  how often to check whether any bins are ready to send (measured in number of reads mapped)

"""


import mappy as mp
from kafka import KafkaProducer
import ujson
import cProfile
from math import floor
from bashplotlib.histogram import plot_hist
import util.logging as my_log
import logging
import datetime
from docopt import docopt


KAFKA_MAX_REQUEST_SIZE=100000000
KAFKA_TOPIC_NAME='mapped_reads_batchy'
KAFKA_BROKER_URL='localhost:9092'

BIN_SIZE=10000
NUM_BINS=100
CHECK_FREQ=1000

def get_partition(rec):
    return 0 #TODO: remove this
    if rec['pos'] > 30000000:
        return 1
    else:
        return 0

def serialize_read_msg(m):
    return ujson.dumps(m).encode("utf-8")


class Reservoir(object):

    def __init__(self, left, right, num_bins, max_bin_size):
        self.max_bin_size = max_bin_size
        self.left = left
        self.right = right
        self.num_bins = num_bins
        self.bin_width = floor((right-left)/num_bins)
        self.bins = [[] for _ in range(num_bins)]

    def put_in_bin(self, obj, place):
        bin_index = self.get_bin_index(place)

        self.bins[bin_index].append(obj)

    def get_bin_index(self, place):
        true_index = floor((place - self.left) / self.bin_width)

        if true_index < 0:
            return 0
        elif true_index >= self.num_bins:
            return self.num_bins - 1
        else:
            return floor((place - self.left) / self.bin_width)

    def empty_bin(self, bin_index):
        self.bins[bin_index] = []

    def empty_bins(self, bin_index_list):
        for bin_index in bin_index_list:
            self.bins[bin_index] = []

    def get_bin_by_index(self, bin_index):
        return self.bins[bin_index]

    def get_all_full_bins(self):
        bins_to_return = []
        for ind, bin in enumerate(self.bins):
            if len(bin) >= self.max_bin_size:
                bins_to_return.append((ind, bin))

        return bins_to_return

    def get_all_bins(self):
        return self.bins

    def get_bin_fullness(self):
        return [len(bin) for bin in self.bins]

def run_mapper(my_start, my_stop, my_contig, ref_path, sample_path, logger):

    logger.info("Loading reference")
    ref = mp.Aligner(ref_path, preset="sr")
    logger.info("Reference Loaded")

    bin_counter = 0
    producer = KafkaProducer(bootstrap_servers=[KAFKA_BROKER_URL],
                             value_serializer=serialize_read_msg,
                             max_request_size=KAFKA_MAX_REQUEST_SIZE)

    accum = []

    hits = []

    counter = 0

    max_bin_size = BIN_SIZE
    num_bins = NUM_BINS
    my_reservoir = Reservoir(my_start, my_stop, num_bins, max_bin_size)
    # fig, ax = plt.subplots(1, 1)
    # ax.set_xlim(0, int(max_bin_size*1.1))
    # plt.show()
    for name, seq, qual in mp.fastx_read(sample_path):
        hit_count = 1

        for hit in ref.map(seq):
            if hit.ctg == my_contig:
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

                my_reservoir.put_in_bin(hit_dict, hit.r_st)

                if(counter % CHECK_FREQ == 0):
                    plot_hist(my_reservoir.get_bin_fullness(), bincount=int(num_bins/10), xlab=True, height=10)
                    # ax.cla()
                    # #ax.set_xlim(0, int(max_bin_size * 1.1))
                    # ax.hist(my_reservoir.get_bin_fullness(), int(num_bins/10))
                    # fig.canvas.draw()
                    full_bins = my_reservoir.get_all_full_bins()
                    bin_indices_to_empty = []
                    for ind, bin in full_bins:
                        producer.send(KAFKA_TOPIC_NAME, bin, partition=get_partition(hit_dict))
                        producer.flush()
                        bin_indices_to_empty.append(ind)
                        bin_counter+=1
                        logger.info(f"Message {bin_counter}. Sending full bin {ind} to the broker")

                    my_reservoir.empty_bins(bin_indices_to_empty)

                # if(len(hits) % 10000 == 0):
                #     producer.send('mapped_reads_batchy', hits, partition=get_partition(hit_dict))
                #     hits = []
                #     print("{} reads processed".format(counter))

                break
            else:
                pass
                #print "Read maps to {}. Discarding".format(hit.ctg)

    for ind, bin in enumerate(my_reservoir.get_all_bins()):
        producer.send(KAFKA_TOPIC_NAME, bin, partition=get_partition(hit_dict))
        logger.info(f"Sending bin {ind} with {len(bin)} items to the broker")
        producer.flush()


    # if(len(hits) > 0):
    #     producer.send('mapped_reads_batchy', hits, partition=get_partition(hits[0]))
    #     producer.flush()
    #     hits = []
    #     print("{} reads processed".format(counter))

    logger.info(counter)
    logger.info(accum)

def main(start, stop, contig, ref_path, sample_path):
    logger = my_log.SetupLogger("mapper")
    logger.setLevel(logging.INFO)
    logger.info(f"now is {datetime.datetime.now()}", )

    run_mapper(start, stop, contig, ref_path, sample_path, logger)
    # cProfile.runctx("run_mapper()", globals(), locals(),
    #                 'profile-batchy-mapper.out')

if __name__ == "__main__":
    my_start = 50000000
    my_stop = 63025520
    my_chr = "20"
    reference_file = "/Users/siakhnin/data/reference/genome.mmi"
    sample_file = "/Users/siakhnin/data/giab/RMNISTHS_30xdownsample_50000000-63025520.name_sorted.fastq"
    args = docopt(__doc__, version="Batchy Mapper 0.0.1")

    if(args['file']):
        sample_file = args['-i']
        reference_file = args['-r']
        region = args['-g']
        split_region = region.split(":")
        my_chr = split_region[0]
        start_stop_region = split_region[1].split("-")
        my_start = start_stop_region[0]
        my_stop = start_stop_region[1]

        BIN_SIZE = args.get('--bin_size', BIN_SIZE)
        NUM_BINS = args.get('--num-bins', NUM_BINS)
        CHECK_FREQ = args.get('--check-freq', CHECK_FREQ)

        KAFKA_BROKER_URL = args['-k']

        main(my_start, my_stop, my_chr, reference_file, sample_file)
    else:
        print("Service mode not yet supported")
# out_file = open("/Users/siakhnin/data/RMNISTHS_30xdownsample_9999999_11000000.mapped.sr.msgpack","w")
# msgpack.dump(hits, out_file, encoding='utf-8')
# out_file.close()