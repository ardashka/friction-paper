from cidlib.InterlacedTime import InterlacedTime
from plot.archive.archive import loadarchive
from plot.plot.lyotropic import get_defects
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
    lattice = np.zeros((frame.parameters['LX'], frame.parameters['LY']), dtype=object)
    for defect in get_defects(frame):
        i, j = (defect.get('pos') - 1.5).astype(int)
        lattice[i, j] = 'p' if defect.get('charge') > 0 else 'n'
    return window(lattice)

def defects(ar):
    return np.stack([ defect_lattice(frame) for frame in ar.read_frames() ])


if __name__ == '__main__':
    
    # ** Parse. Get the input for the simulation. **
    parser = argparse.ArgumentParser(description='Parameters for CID correlation:')
    parser.add_argument("-o", dest="path2archive", help="Path to output directory.", type=Path, required=True)
    parser.add_argument("-nshuff", dest="nshuff", help="Number of shuffles. Default value = 4", type=int, default=4)
    
    args = parser.parse_args()
    
    ## ** Get Archive **
    ar = loadarchive(args.path2archive)
    order = get_order(ar)
    ds = defects(ar)
    
    # ** CID analysis **
    it = InterlacedTime(order, args.nshuff, True)
    
    output = {'zeta': ar.zeta, 
              'Q_ITC_distinguishable': it( ds ),
              'Q_ITC_indistinguishable': it( np.where(ds=='n', 'p', ds) )}
    
    with open(Path.joinpath(args.path2archive.parent, 'CID_correlation.json'), 'w') as outfile:
        json.dump(output, outfile)
