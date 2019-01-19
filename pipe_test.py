from multiprocessing import Process, Pipe
from apscheduler.schedulers.background import BackgroundScheduler
import datetime

def do_save(storer_pipe):
    print("Starting saving at {}".format(datetime.datetime.now()))


    if storer_pipe.poll():

        msg = storer_pipe.recv()
        print(msg)
        if msg == "NEW":
            print("Saving happens here")
            storer_pipe.send("SAVED")
            print("Notifying processor.")
        else:
            print("No new messages reported. Not saving this time.")
    else:
        print("No messages from process. Not saving this time.")


def saver_job(storer_pipe):
    p = Process(target=do_save, args=[storer_pipe])
    p.start()
    p.join()


def the_processor(processor_pipe):
    print("Starting Message Processing")

    saver_notified = False

    while True:
        #print("Got a message here.")

        if processor_pipe.poll():
            msg = processor_pipe.recv()
            if msg == "SAVED":
                saver_notified = False

        if not saver_notified:
            processor_pipe.send("NEW")
            saver_notified = True
            print("Notifying saver")


processor_pipe, storer_pipe = Pipe()

scheduler = BackgroundScheduler()
job = scheduler.add_job(saver_job, 'interval', args=[storer_pipe],seconds=10, )
scheduler.start()

p = Process(target=the_processor, args=[processor_pipe])
p.start()
p.join()


