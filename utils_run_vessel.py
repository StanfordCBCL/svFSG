from mpi4py import MPI
import numpy as np
import pickle
import os
from ctypes import cdll
from ctypes import c_char_p
from ctypes import c_int
from ctypes import c_double

lib = cdll.LoadLibrary('./libgnr.so')
run = lib.run
run.argtypes = [c_char_p, c_char_p, c_int, c_int, c_int, c_int, \
                c_double,c_double,c_double, c_double, c_double, \
                c_double,c_double,c_double,c_double,c_double, \
                c_double,c_double,c_double,c_double]
#run(dir_prefix.encode('ascii'), filename.encode('ascii'),restart_arg,iter_arg,gnr_arg,num_days,step_size,sigma_arg,tauw_arg,anysm_arg,tevg_arg,F0, F1, F2, F3, F4,F5, F6, F7, F8)
run.restype = c_int

comm = MPI.COMM_WORLD
size = comm.Get_size() # new: gives number of ranks in comm
rank = comm.Get_rank()

input_file = open("input_array.dat","rb")
input_array = pickle.load(input_file)
input_file.close()
numData = len(input_array)

if rank == 0:


    array_list = np.arange(numData)
    output_array = np.zeros([numData,59]) #Create output array of same size
    split = np.array_split(array_list,size) #Split input array by the number of available cores

    split_sizes = []

    for i in range(0,len(split),1):
        split_sizes = np.append(split_sizes, len(split[i]))

    split_sizes_input = split_sizes
    displacements_input = np.insert(np.cumsum(split_sizes_input),0,0)[0:-1]

    split_sizes_output = split_sizes*59
    displacements_output = np.insert(np.cumsum(split_sizes_output),0,0)[0:-1]


else:
#Create variables on other cores
    split_sizes_input = None
    displacements_input = None
    split_sizes_output = None
    displacements_output = None
    split = None
    array_list = None
    output_array = None


split = comm.bcast(split, root=0) #Broadcast split array to other cores
split_sizes = comm.bcast(split_sizes_input, root = 0)
displacements = comm.bcast(displacements_input, root = 0)
split_sizes_output = comm.bcast(split_sizes_output, root = 0)
displacements_output = comm.bcast(displacements_output, root = 0)

output_chunk = np.array([0] * len(split[rank])) #Create array to receive subset of data on each core, where rank specifies the core
comm.Scatterv([array_list,split_sizes_input, displacements_input,MPI.DOUBLE],output_chunk,root=0)

output = np.zeros([len(output_chunk),59]) #Create output array on each core

for i in range(0,np.shape(output_chunk)[0],1):
    inputData = input_array[output_chunk[i]]
    out_array_type = c_double * 59  # equiv. to C double[3] type
    out_array = out_array_type(*([0]*59))        # equiv. to double arr[3] = {...} instance
    prefix = inputData[2].encode('ascii')
    suffix = ('python_'+str(inputData[0])+'_'+str(inputData[1])).encode('ascii')
    run(prefix, suffix, inputData[3], inputData[4], inputData[5],inputData[6], inputData[7],\
        inputData[8], inputData[9], inputData[10], inputData[11], inputData[12], inputData[13],\
        inputData[14],inputData[15],inputData[16],inputData[17],inputData[18], inputData[19],\
        inputData[20],out_array)
    output[i] = np.array(out_array)

comm.Barrier()

comm.Gatherv(output,[output_array,split_sizes_output,displacements_output,MPI.DOUBLE], root=0) #Gather output data together

if rank == 0:
    parse_file = open("parse_array.dat","wb")
    pickle.dump(output_array, parse_file)
    parse_file.close()

