import Queue, multiprocessing
from multiprocessing.managers import SyncManager

IP = "192.168.0.178"
PORTNUM = 5000
AUTHKEY = "abc"

def mp_factorizer(shared_job_q, shared_result_q, batchId, nprocs):
    """ Split the work with jobs in shared_job_q and results in
        shared_result_q into several processes. Launch each process with
        factorizer_worker as the worker function, and wait until all are
        finished.
    """
    id = batchId()
    print id.value

def make_client_manager(ip, port, authkey):
    """ Create a manager for a client. This manager connects to a server on the
        given address and exposes the get_job_q and get_result_q methods for
        accessing the shared queues from the server.
        Return a manager object.
    """
    class ServerQueueManager(SyncManager):
        pass

    ServerQueueManager.register('get_job_q')
    ServerQueueManager.register('get_result_q')
    ServerQueueManager.register('get_batchId')

    manager = ServerQueueManager(address=(ip, port), authkey=authkey)
    manager.connect()

    print 'Client connected to %s:%s' % (ip, port)
    return manager

def runclient():
    manager = make_client_manager(IP, PORTNUM, AUTHKEY)
    job_q = manager.get_job_q()
    result_q = manager.get_result_q()
    mp_factorizer(job_q, result_q,manager.get_batchId, 4)
    
if __name__=='__main__':
    runclient()
    