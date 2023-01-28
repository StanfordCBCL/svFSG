import ctypes

lib = ctypes.cdll.LoadLibrary('FolderExecutables/libgnr.so')

class cvessel():
    def __init__(self):
        # Declare input and output types for each method you intend to use
        lib.init.argtypes = [ctypes.c_void_p]
        lib.init.restype = ctypes.c_void_p

        lib.initializeVesselHandshake.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double]
        lib.initializeVesselHandshake.restype = ctypes.c_void_p

        lib.updateVesselHandshake.argtypes = [ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_double,ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double)]
        lib.updateVesselHandshake.restype = ctypes.c_int

        self.id = None

        self.obj = lib.init(1)

    def initializeVesselHandshake(self, prefix, id_num, num_days, step_size, anysm_arg, tevg_arg):
        prefix = prefix.encode('ascii')
        suffix = ('python_'+str(id_num)).encode('ascii')
        self.id = id_num
        lib.initializeVesselHandshake(self.obj, prefix, suffix, num_days, step_size, anysm_arg, tevg_arg)
        return
    
    def updateVesselHandshake(self, restart_arg, iter_arg, sigma_arg, tauw_arg, F):
        F_array_type = ctypes.c_double * 9
        out_array_type = ctypes.c_double * 59  # equiv. to C double[3] type
        out_array = out_array_type(*([0]*59))        # equiv. to double arr[3] = {...} instance
        F_array = F_array_type(*F)
        lib.updateVesselHandshake(self.obj, restart_arg, iter_arg, sigma_arg, tauw_arg, F_array, out_array)  # pointer to array passed to function and modified
        return list(out_array)
