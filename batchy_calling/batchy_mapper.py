import mappy as mp
from kafka import KafkaProducer
from json import dumps
import msgpack
import ujson
import cProfile
from math import floor
from bashplotlib.histogram import plot_hist
import matplotlib.pyplot as plt

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
        self.bins = [[] for i in range(num_bins)]

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

def run_mapper():
    my_start = 50000000
    my_stop = 63025520
    my_chr = "20"
    reference_file = "/Users/siakhnin/data/reference/genome.mmi"
    sample_file = "/Users/siakhnin/data/giab/RMNISTHS_30xdownsample_50000000-63025520.name_sorted.fastq"
    print("Loading reference")
    ref = mp.Aligner(reference_file, preset="sr")
    print("Reference Loaded")
    bin_counter = 0
    producer = KafkaProducer(bootstrap_servers=['localhost:9092'],
                             value_serializer=serialize_read_msg,
                             max_request_size=100000000)

    accum = []

    hits = []

    counter = 0

    max_bin_size = 10000
    num_bins = 100
    my_reservoir = Reservoir(my_start, my_stop, num_bins, max_bin_size)
    # fig, ax = plt.subplots(1, 1)
    # ax.set_xlim(0, int(max_bin_size*1.1))
    # plt.show()
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

                my_reservoir.put_in_bin(hit_dict, hit.r_st)

                if(counter % 10000 == 0):
                    plot_hist(my_reservoir.get_bin_fullness(), bincount=int(num_bins/10), xlab=True, height=10)
                    # ax.cla()
                    # #ax.set_xlim(0, int(max_bin_size * 1.1))
                    # ax.hist(my_reservoir.get_bin_fullness(), int(num_bins/10))
                    # fig.canvas.draw()
                    full_bins = my_reservoir.get_all_full_bins()
                    bin_indices_to_empty = []
                    for ind, bin in full_bins:
                        producer.send('mapped_reads_batchy', bin, partition=get_partition(hit_dict))
                        producer.flush()
                        bin_indices_to_empty.append(ind)
                        bin_counter+=1
                        print(f"Message {bin_counter}. Sending full bin {ind} to the broker")

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
        producer.send('mapped_reads_batchy', bin, partition=get_partition(hit_dict))
        print(f"Sending bin {ind} with {len(bin)} items to the broker")
        producer.flush()


    # if(len(hits) > 0):
    #     producer.send('mapped_reads_batchy', hits, partition=get_partition(hits[0]))
    #     producer.flush()
    #     hits = []
    #     print("{} reads processed".format(counter))

    print(counter)
    print(accum)

def main():
    cProfile.runctx("run_mapper()", globals(), locals(),
                    'profile-batchy-mapper.out')

if __name__ == "__main__":
    main()

# out_file = open("/Users/siakhnin/data/RMNISTHS_30xdownsample_9999999_11000000.mapped.sr.msgpack","w")
# msgpack.dump(hits, out_file, encoding='utf-8')
# out_file.close()