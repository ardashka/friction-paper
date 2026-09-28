from itertools import chain
import numpy as np

def hilbert_curve(p, d):
    """ Principle Hilbert curve. \n
    This is a Pythonized version of the **principle** Hilbert curve implementation
    by John Skilling as described in:
    
    Skilling, J. (2004, April). Programming the Hilbert curve. In AIP Conference
    Proceedings (Vol. 707, No. 1, pp. 381-387). American Institute of Physics.

    Algorithm:
    A single global Gray code [i.e. Reflected binary Code (RBC)] is being applied
    to the pd-bit binary representation of the Hilbert distance/index H. This over-
    transforms the distance H and the **excess work** is undone by a single traverse
    through the pd-bit Gray code representation of H.

    Args:
        p: The number of bits for each dimension (int).
        d: The dimensionality of the hypercube (int).

    Returns:
        A list of tuples: The i'th tuple contains the coordinates (int) to the
        i'th point along the principle Hilbert curve.
    """
    if type(d) is not int:
        raise TypeError('dim must be an integer, got {}.'.format(type(p)))
    elif d <= 1:
        raise ValueError('dim must be >= 2, got %d.' %(d))
    
    if type(p) is not int:
        raise TypeError('order must be an integer, got {}.'.format(type(p)))
    elif p <= 0:
        raise ValueError('order must be >= 1, got %d.' %(p))
    
    def points_from_index():
        for gray in [idx ^ (idx >> 1) for idx in range(1 << p*d)]:
            point = [extract_bits(gray, n) for n in range(d)]
            yield undo_excess_work(point)
    
    def extract_bits(num, n):
        binary = format(num, '0'+str(p*d)+'b')
        extracted_bits = binary[n::d]
        return int(extracted_bits, 2)

    def undo_excess_work(point):
        for pot in (2 << p for p in range(p-1)):
            for n in reversed(range(d)):
                if (point[n] & pot) != 0:
                    point[0] ^= pot - 1
                else:
                    exchange_low_bits = (point[0] ^ point[n]) & pot - 1
                    point[0] ^= exchange_low_bits
                    point[n] ^= exchange_low_bits
        return tuple(point)
    
    return [point for point in points_from_index()]


def load_hilbert_curves(p, d):
    """ Loads all 8 possible hilbert curves. \n
    Args:
        p: The number of bits for each dimension (int).
        d: The dimensionality of the hypercube (int).
    Returns:
        A list of of 8 lists of tuples: The i'th tuple in each list, contains the
        coordinates (int) to the i'th point along the respective Hilbert curve.
    """
    
    if d not in [2, 3]:
        raise ValueError('dim must be either 2 or 3, got %d.' %(d))
    
    def rotate_hilbert_curve(H, k):
        """ Rotate Hilbert curve by 90 degrees. \n
        Args:
            H: Hilbert curve (must be either 2D or 3D).
            k: Number of times the Hilbert curve is rotated by 90 degrees.
        Returns:
            A rotated copy of H and its swapped counterpart. 
        """
        def swap(hcurve):
            """ swap last two columns """
            hcurve[:, [-2, -1]] = hcurve[:, [-1, -2]]
            return hcurve
        
        def flip(hcurve, axis):
            """ flip along axis """
            hcurve[:, axis] = (1 << p) - 1 - hcurve[:, axis]
            return hcurve
        
        def rot(hcurve, k):
            """ rotate k times by 90 degrees """
            k %= 4
            if k == 0:
                return hcurve
            elif k == 1:
                return swap(flip(hcurve, -2))
            elif k == 2:
                return flip(flip(hcurve, -2), -1)
            elif k == 3:
                return flip(swap(hcurve), -2)
        
        H_rot = rot(np.copy(H), k)
        H_swap = swap(np.copy(H_rot))
        
        return list(map(tuple, H_rot)), list(map(tuple, H_swap))
    
    HPC = hilbert_curve(p, d)   # principle hilbert curve
    
    return list(chain.from_iterable( rotate_hilbert_curve(HPC, k) for k in range(4) ))
