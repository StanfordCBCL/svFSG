from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
size = comm.Get_size() # new: gives number of ranks in comm
rank = comm.Get_rank()

if rank == 0:
    np.savetxt('numrank', [size], fmt='%i')