import zlib
import pickle

class cvessel():
    def __init__(self):
        self.prefix = ""
        self.name = ""
        self.restart = 0
        self.iteration = 0
        self.simulate = 0
        self.num_days = 0
        self.step_size = 0
        self.sigma_inv = 0
        self.tauw_wss = 0
        self.aneurysm = 0
        self.tevg = 0
        self.F = [0]*9
        self.out_array = [0]*59
        self.savearray = zlib.compress(pickle.dumps([]))