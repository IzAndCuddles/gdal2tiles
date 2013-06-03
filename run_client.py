import Queue, multiprocessing
from multiprocessing.managers import SyncManager

IP = "192.168.0.178"
PORTNUM = 5000
AUTHKEY = "abc"

def factorize_naive(n):
    """ A naive factorization method. Take integer 'n', return list of
        factors.
    """
    if n < 2:
        return []
    factors = []
    p = 2

    while True:
        if n == 1:
            return factors

        r = n % p
        if r == 0:
            factors.append(p)
            n = n / p
        elif p * p >= n:
            factors.append(n)
            return factors
        elif p > 2:
            # Advance in steps of 2 over odd numbers
            p += 2
        else:
            # If p == 2, get to 3
            p += 1
    assert False, "unreachable"

def factorizer_worker(job_q,result_q,entier,event):
    """ A worker function to be launched in a separate process. Takes jobs from
        job_q - each job a list of numbers to factorize. When the job is done,
        the result (dict mapping number -> list of factors) is placed into
        result_q. Runs until job_q is empty.
    """
    e = entier._getvalue()
    i=0
    while True:
        print e, entier._getvalue(),event.is_set()
        i+=1
        while job_q.empty():
            if event.is_set() or entier._getvalue() != e:
                print "Done ", i
                break
            else:
                event.wait(0.1)
        job = job_q.get()
        outdict = dict((n,factorize_naive(n)) for n in job)
        job_q.task_done()
        result_q.put(outdict)



def mp_factorizer(shared_job_q,shared_result_q,entier,event, nprocs):
    """ Split the work with jobs in shared_job_q and results in
        shared_result_q into several processes. Launch each process with
        factorizer_worker as the worker function, and wait until all are
        finished.
    """
    procs = []
    for i in range(nprocs):
        p = multiprocessing.Process(
                target=factorizer_worker,
                args=(shared_job_q, shared_result_q,entier,event))
        procs.append(p)
        p.start()

    done = False
    while not done:
        alive_count = len(procs)
        for p in procs:
            p.join(1.0)
            alive = p.is_alive()
            if not alive:
                alive_count -= 1
            print "entier :",entier,", entier(_getvalue) :", entier._getvalue(),", event :", event.is_set(), ", alive :", alive_count
        done = (alive_count == 0)

def make_client_manager(ip, port, authkey):
    """ Create a manager for a client. This manager connects to a server on the
        given address and exposes the get_job_q and get_result_q methods for
        accessing the shared queues from the server.
        Return a manager object.
    """
    class ServerQueueManager(SyncManager):
        pass

    ServerQueueManager.register('get_job_q')
    ServerQueueManager.register('get_entier')
    ServerQueueManager.register('get_result_q')
    ServerQueueManager.register('get_event')

    manager = ServerQueueManager(address=(ip, port), authkey=authkey)
    manager.connect()

    print 'Client connected to %s:%s' % (ip, port)
    return manager

def runclient():
    manager = make_client_manager(IP, PORTNUM, AUTHKEY)
    job_q = manager.get_job_q()
    result_q = manager.get_result_q()
    entier = manager.get_entier()
    event = manager.get_event()
    
    mp_factorizer(job_q,result_q,entier,event, 4)
    
if __name__=='__main__':
    runclient()
    