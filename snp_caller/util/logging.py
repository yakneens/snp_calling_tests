import os
import logging
import time

def SetupLogger(name_prefix):
    if not os.path.exists("log"):
        os.makedirs("log")

    recfmt = logging.Formatter('%(asctime)s.%(msecs)03d %(levelname)s %(message)s')

    handler = logging.FileHandler(time.strftime(f"log/{name_prefix}.%y%m%d.log"))
    handler.setFormatter(recfmt)
    handler.setLevel(logging.DEBUG)

    logger = logging.getLogger(f"{name_prefix} {__name__}")
    logger.addHandler(handler)

    return logger