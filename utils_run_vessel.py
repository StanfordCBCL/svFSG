from mpi4py import MPI
import numpy as np
import pickle
from cvessel import cvessel
import os
import ctypes
import zlib


lib = ctypes.cdll.LoadLibrary('./libgnr.so')
run = lib.run
run.restype = ctypes.c_char_p #ctypes.POINTER(ctypes.c_char)
run.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, \
                ctypes.c_double,ctypes.c_double,ctypes.c_double, ctypes.c_double, ctypes.c_double, \
                ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double, \
                ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double]


comm = MPI.COMM_WORLD
size = comm.Get_size() # new: gives number of ranks in comm
rank = comm.Get_rank()


if rank == 0:
    cvessel_file = open("cvessel_array.dat","rb")
    input_array = np.array(pickle.load(cvessel_file))
    cvessel_file.close()
    numData = len(input_array)

    split = np.array_split(input_array,size) #Split input array by the number of available cores

else:
#Create variables on other cores
    split = None

data = comm.scatter(split,root=0)

for i in range(0,np.shape(data)[0],1):

    out_array_type = ctypes.c_double * 59  # equiv. to C double[3] type
    out_array = out_array_type(*([0]*59))  # equiv. to double arr[3] = {...} instance

    savearray = pickle.loads(zlib.decompress(data[i].savearray))

    loadstring = ""
    array_length = np.shape(savearray)[0]
    for j in range(array_length):
        number_length = np.shape(savearray[j])[0]
        for k in range(number_length):
            number = savearray[j][k]
            if number.is_integer():
                loadstring+=str(int(number))
            else:
                loadstring+=str(number)
            if k != number_length - 1:
                loadstring+=" "
        if j != array_length:
            loadstring+='\n'

    loadstring = loadstring.encode('utf-8')
    prefix = data[i].prefix.encode('utf-8')
    suffix = data[i].name.encode('utf-8')
    
    buffersize = 4 * 128 * 1024
    allocation = ctypes.create_string_buffer(buffersize)

    savestring = run(loadstring, prefix, suffix, data[i].restart, data[i].iteration, data[i].simulate, \
                    data[i].num_days, data[i].step_size, data[i].sigma_inv, data[i].tauw_wss, \
                    data[i].aneurysm, data[i].tevg, data[i].F[0], data[i].F[1], data[i].F[2], \
                    data[i].F[3], data[i].F[4], data[i].F[5], data[i].F[6], data[i].F[7], \
                    data[i].F[8], out_array, allocation, buffersize)

    if savestring is None:
        raise RuntimeError('No return from vessel executable. Check return buffer size in utils_run_vessel.py')

    savestring_decode = savestring.decode("utf-8")
    savearray = [[float(digit) for digit in line.split()] for line in savestring_decode.splitlines()]

    data[i].out_array = np.array(out_array)
    data[i].savearray = zlib.compress(pickle.dumps(savearray))

comm.Barrier()

output_array = comm.gather(data,root=0)

if rank == 0:
    cvessel_array = []

    for i in output_array:
        for j in i:
            cvessel_array.append(j)

    cvessel_file = open("cvessel_array.dat","wb")
    pickle.dump(cvessel_array, cvessel_file)
    cvessel_file.close()
