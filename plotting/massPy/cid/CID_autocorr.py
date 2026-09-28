import sys
sys.path.insert(0,'/astro/bhan/Mass/')
from massPy.nematic.nematicPy import get_defects
from massPy.archive import loadarchive
from massPy.cid.cidlib.cid import cid
from pathlib import Path
import argparse, json
import numpy as np

def get_order(ar):
    if ar.LX != ar.LY:
        raise ValueError('Spatial size must be equal, got (LX, LY) = (%d, %d).' %(ar.LX, ar.LY))
    elif ar.LX % 2 != 0:
        raise ValueError('Spatial size must be a power of 2, got LX = %d.' %(ar.LX))
    return int( np.log2(ar.LX) - 1 )

def window(lattice):
    LX, LY = lattice.shape
    x_slice = slice(LX // 4, 3*LX // 4)
    y_slice = slice(LY // 4, 3*LY // 4)
    return lattice[x_slice, y_slice]

def defect_lattice(frame):
    lattice = np.zeros((frame.parameters['LX'], frame.parameters['LY']), dtype=int)
    for defect in get_defects(frame.QQxx, frame.QQyx, frame.LX, frame.LY):
        i, j = (defect.get('pos') - 1.5).astype(int)
        # lattice[i, j] = 1 if defect.get('charge') > 0 else 2
        lattice[i, j] = 1   # INDISTINGUISHABLE
    return window(lattice)

def defects(ar):
    return np.column_stack([ defect_lattice(frame).flatten() for frame in ar.read_frames() ])

if __name__ == '__main__':
    
    # ** Parse. Get the input for the simulation. **
    parser = argparse.ArgumentParser(description='Parameters for CID analysis:')
    parser.add_argument("-o", dest="path2archive", help="Path to output directory.", type=Path, required=True)
    parser.add_argument("-nshuff", dest="nshuff", help="Number of shuffles. Default value = 4", type=int, default=4)
    
    args = parser.parse_args()
    
    ## ** Get Archive **
    ar = loadarchive(args.path2archive)
    order = get_order(ar)
    ds = defects(ar)
    
    res = np.column_stack([ cid(d, args.nshuff) for d in ds ])
    avgs = np.mean(res, axis=1)
    stds = np.std(res, axis=1)
    
    print(res.size)
    
    output = {'zeta': ar.zeta, 
              'CID_avg': avgs[0].tolist(),
              'CID_std': stds[0].tolist(),
              'Q_avg': avgs[1].tolist(),
              'Q_std': stds[1].tolist()
              }
    
    with open(Path.joinpath(args.path2archive.parent, 'CID_autocorr.json'), 'w') as outfile:
        json.dump(output, outfile)
