import timeit
from kafka import KafkaConsumer, KafkaProducer
import redis
import time
import datetime
import pickle
import zlib
import cProfile

SKIP_ALL = False

def run_locus_saver():
    redis_cli = redis.Redis(host='localhost',port=6379)
    redis_pipe = redis_cli.pipeline()

    locus_source = KafkaConsumer('snp_loci_batchy',
                                 group_id='locus_saver',
                                 bootstrap_servers=['localhost:9092'],
                                 value_deserializer=lambda m: pickle.loads(zlib.decompress(m)))
    msg_counter = 0
    for message in locus_source:
        if SKIP_ALL:
            continue
        ts = time.time()
        my_locus_batch = message.value
        for pos, my_locus in my_locus_batch.items():
            redis_pipe.set(pos, pickle.dumps(my_locus))
        redis_pipe.execute()
        te = time.time()
        msg_counter+=1
        print(f"Message {msg_counter}. {datetime.datetime.now()} - Processed {len(my_locus_batch)} loci in {te-ts} secs")

def main():
    cProfile.runctx("run_locus_saver()", globals(), locals(),
                    'profile-locus-saver.out')

if __name__ == "__main__":
    main()
