from lempel_ziv_complexity import lempel_ziv_complexity as lz77
from itertools import permutations, product
from HilbertCurve import HilbertCurve
import multiprocessing as mp
from statistics import mean
from random import shuffle
import numpy as np

class InterlacedTime():

    def __init__(self, lattice):

        def is_hypercubic(lattice):
            return all([i == lattice.shape[0] for i in lattice.shape[1:]])
        
        def is_pot(lattice):
            return (len(lattice) & (len(lattice)-1) == 0) and len(lattice) != 0

        if type(lattice) is not np.ndarray:
            raise TypeError('Lattice must be a numpy array, got {}.'.format(type(lattice)))
        elif lattice.ndim < 2:
            raise ValueError('Dimensions of numpy array must be at least 2, got %d' %(lattice.ndim))
        elif not is_hypercubic(lattice):
            raise ValueError('Lattice must be hypercubic, got {}.'.format(lattice.shape))
        elif not is_pot(lattice):
            raise ValueError('Length of each lattice dimension must be a power of two!')

        self.lattice = lattice
        self.order = int(np.frexp(len(lattice))[1] - 1)
        self.dim = lattice.ndim
        self.hilbert = HilbertCurve(self.order, self.dim)
        self.Q = self.isotropic_q_order()


    def isotropic_q_order(self):
        pool = mp.Pool(mp.cpu_count())
        qs = pool.map_async(self.q_order, self.full_octahedral_group())
        pool.close()
        pool.join()
        return mean(qs.get())
        

    def full_octahedral_group(self):
        for axes in permutations(range(self.dim)):
            for f0, f1, f2 in product([1, -1], repeat=self.dim):
                yield np.transpose(self.lattice[::f0, ::f1, ::f2], axes)


    def q_order(self, lattice):

        hilbert_scan = [lattice[point] for point in self.hilbert]

        def list2string(mylist):
            return ''.join(map(str, mylist))

        def cid(hscan):
            C = lz77(list2string(hscan))
            L = 1 << (self.dim * self.order)
            return (C*np.log2(C) + 2*C*np.log2(L/C)) / L
        
        def cid_shuffles(hscan):
            rand_shuffles = []
            for _ in range(1 << self.order):
                shuffle(hscan)
                rand_shuffles.append(cid(hscan))
            return mean(rand_shuffles)

        return 1 - cid(hilbert_scan) / cid_shuffles(hilbert_scan)