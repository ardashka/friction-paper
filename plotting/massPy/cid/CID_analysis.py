import sys
sys.path.insert(0,'/astro/bhan/Mass/')
from massPy.cid.cidlib.interlacings import InterlacedTime, SequentialTime
from massPy.nematic.nematicPy import get_defects
from massPy.archive import loadarchive
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
    lattice = np.zeros((frame.LX, frame.LY), dtype=int)
    for defect in get_defects(frame.QQxx, frame.QQyx, frame.LX, frame.LY):
        i, j = defect['pos'].astype(int) % (frame.LX, frame.LY)
        # lattice[i, j] = 1 if defect.get('charge') > 0 else 2
        lattice[i, j] = 1   # INDISTINGUISHABLE
    return window(lattice)

def defects(ar):
    return np.stack([ defect_lattice(frame) for frame in ar.read_frames() ])

if __name__ == '__main__':
    
    # ** Parse. Get the input for the simulation. **
    parser = argparse.ArgumentParser(description='Parameters for CID analysis:')
    parser.add_argument("-o", dest="path2archive", help="Path to output directory.", type=Path, required=True)
    parser.add_argument("-nshuff", dest="nshuff", help="Number of shuffles. Default value = 4", type=int, default=4)
    # parser.add_argument("-nsample", dest="nsample", help="Number of frames to skip in CID_SEQ. Default value = 16", type=int, default=16)
    
    args = parser.parse_args()
    
    print('PATH2ARCHIVE : ', args.path2archive)
    
    ## ** Get Archive **
    ar = loadarchive(args.path2archive)
    order = get_order(ar)
    ds = defects(ar)
    
    print('--- STARTING CID ANALYSIS ---')
    
    # ** CID analysis **
    it = InterlacedTime(order, args.nshuff)
    sq = SequentialTime(order, args.nshuff)
    
    cid_itc, q_itc = it(ds)
    cid_seq, q_seq = sq(ds)
    
    output = {'zeta': ar.zeta, 
              'CID_ITC': cid_itc,
              'CID_SEQ': cid_seq,
              'Lambda_ITC': q_itc,
              'Lambda_SEQ': q_seq,
              'NUM_DEFECTS': np.sum(ds, axis=(1,2)).tolist()
              }
    
    print('OUTPUT : ', output)
    
    with open(Path.joinpath(args.path2archive.parent, 'CID_analysis.json'), 'w') as outfile:
        json.dump(output, outfile)
        
    print('--- COMPLETE ---')
