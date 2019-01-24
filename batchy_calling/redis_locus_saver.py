
import timeit
from kafka import KafkaConsumer, KafkaProducer
import redis
import time
import datetime
import pickle
import zlib


redis_cli = redis.Redis(host='localhost',port=6379)
redis_pipe = redis_cli.pipeline()

locus_source = KafkaConsumer('snp_loci_batchy',
                             group_id='locus_saver',
                             bootstrap_servers=['localhost:9092'],
                             value_deserializer=lambda m: pickle.loads(zlib.decompress(m)))

for message in locus_source:
    ts = time.time()
    my_locus_batch = message.value
    for pos, my_locus in my_locus_batch.items():
        redis_pipe.set(pos, pickle.dumps(my_locus))
    redis_pipe.execute()
    te = time.time()
    print(f"{datetime.datetime.now()} - Processed {len(my_locus_batch)} loci in {te-ts} secs")
