from lempel_ziv_complexity import lempel_ziv_complexity as lz77
from HilbertCurve import HilbertCurve
from numpy import ndarray, log2

class ComputableInformationDensity():

    def __init__(self, p, n):
        self.hilbert = HilbertCurve(p, n)

    def cid(self, lattice):

        if type(lattice) is not ndarray:
            raise TypeError('Lattice must be a numpy array, got {}.'.format(type(lattice)))
        elif len(lattice.shape) != self.hilbert.dim:
            raise ValueError('Dimension of lattice must be %d, got %d.' %(self.hilbert.dim, len(lattice.shape)))
        elif not all([i == lattice.shape[0] for i in lattice.shape[1:]]):
            raise ValueError('Lattice must be hypercubic, got {}.'.format(lattice.shape))
        elif len(lattice) != (1 << self.hilbert.order):
            raise ValueError('Length of each lattice dimension must be %d, got %d' %(1 << self.hilbert.order, len(lattice)))

        hilbert_scan = [lattice[tuple(point)] for point in self.hilbert]
        C = lz77(self.list2string(hilbert_scan))
        L = 1 << self.hilbert.bits

        return (C*log2(C) + 2*C*log2(L/C)) / L
    
    @staticmethod
    def list2string(mylist):
        return str(mylist).replace('[', '').replace(', ', '').replace(']', '')