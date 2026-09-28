from .hilbert_curve import load_hilbert_curves
from abc import ABC, abstractmethod
from multiprocessing import Pool
from .cid import cid
import numpy as np

class CID(ABC):
    
    def __init__(self, order, nshuff):
        self. order = order
        self.nshuff = nshuff
    
    @abstractmethod
    def cid_worker(self, hscan):
        pass
    
    def cid_parallel(self, hscans):
        with Pool(8) as p:
            res = p.map_async(self.cid_worker, hscans)
            p.close()
            p.join()
        
        return np.mean(res.get(), axis=0)
    

class InterlacedTime(CID):
    
    def __init__(self, order, nshuff):
        super().__init__(order, nshuff)
        
    def cid_worker(self, hscan):
        return cid(hscan, self.nshuff)
    
    def __call__(self, data):
        if type(data) is not np.ndarray:
            raise TypeError('Data must be a numpy array, got {}.'.format(type(data)))
        elif data.ndim != 3:
            raise ValueError('Dimensions of numpy-data-array must be 3, got %d' %(data.ndim))
        elif data.shape[1:] != (1 << self.order, ) * 2:
            raise ValueError('Spatial size must be {}, got {}.'.format((1 << self.order, ) * 2, data.shape[1:]))
        elif data.shape[0] % data.shape[1] != 0:
            raise ValueError('Temporal size must be a multiplum of the Spatial size; data.shape = {}'.format(data.shape))
        
        def split_data():
            temporal_sections = data.shape[0] // (1 << self.order)
            return np.split(data, temporal_sections, axis=0)
        
        def hilbert_scans():
            for hcurve in load_hilbert_curves(self.order, 3):
                yield [ cube[point] for cube in split_data() for point in hcurve ]
        
        return self.cid_parallel( hilbert_scans() ).tolist()


class SequentialTime(CID):
    
    def __init__(self, order, nshuff):
        super().__init__(order, nshuff)
    
    def cid_worker(self, hscans):
        return [cid(hscan, self.nshuff) for hscan in hscans]
    
    def __call__(self, data):
        if type(data) is not np.ndarray:
            raise TypeError('Data must be a numpy array, got {}.'.format(type(data)))
        elif data.ndim != 3:
            raise ValueError('Dimensions of numpy-data-array must be 3, got %d' %(data.ndim))
        elif data.shape[1:] != (1 << self.order, ) * 2:
            raise ValueError('Spatial size must be {}, got {}.'.format((1 << self.order, ) * 2, data.shape[1:]))
        
        def hilbert_scans():
            for hilbert_curve in load_hilbert_curves(self.order, 2):
                yield [ [lattice[point] for point in hilbert_curve] for lattice in data ]
        
        res = super().cid_parallel( hilbert_scans() )
        
        return list(zip( np.mean(res, axis=0), np.std(res, axis=0) ))
