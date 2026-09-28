class HilbertCurve():
    """
    This is a Pythonized version of the **principle** Hilbert curve implementation
    by John Skilling as described in:
    
    Skilling, J. (2004, April). Programming the Hilbert curve. In AIP Conference
        Proceedings (Vol. 707, No. 1, pp. 381-387). American Institute of Physics.

    Algorithm:
    ----------
    A single global Gray code [i.e. Reflected binary Code (RBC)] is being applied
    to the np-bit binary representation of the Hilbert distance/index H. This over-
    transforms the distance H and the **excess work** is undone by a single traverse
    through the np-bit Gray code representation of H.

    Parameters:
    -----------
        p   -   The number of bits for each dimension (integer).
        n   -   The dimensionality of the hypercube (integer).

    Return:
    -------
        points  -   A list of lists, where the i'th sublist contains the coordinates
                    (integers) to the i'th point along the Hilbert curve.
    """

    def __init__(self, p, n):

        if type(p) is not int:
            raise TypeError('p must be an integer: Got type {}.'.format(type(p)))
        elif p <= 0:
            raise ValueError('p must be >= 1. Got p = %d.' %(p))

        if type(n) is not int:
            raise TypeError('n must be an integer: Got type {}.'.format(type(n)))
        elif n <= 1:
            raise ValueError('n must be >= 2. Got n = %d.' %(n))

        self.order = p
        self.dim = n
        self.bits = n*p
        self.points = self.points_from_distances()

    def points_from_distances(self):
        # points = []
        for gray in [dist ^ (dist >> 1) for dist in range(pow(2, self.bits))]:
            point = [self.extract_bits(gray, n) for n in range(self.dim)]
            yield self.undo_excess_work(point)
            # points.append(self.undo_excess_work(point))
        # return points
    
    def extract_bits(self, num, n):
        binary = format(num, '0'+str(self.bits)+'b')
        extracted_bits = binary[n::self.dim]
        return int(extracted_bits, 2)

    def undo_excess_work(self, point):
        for p in range(self.order - 1):
            pot = 2 << p
            for n in reversed(range(self.dim)):
                if (point[n] & pot) != 0:
                    point[0] ^= pot - 1
                else:
                    exchange_low_bits = (point[0] ^ point[n]) & pot - 1
                    point[0] ^= exchange_low_bits
                    point[n] ^= exchange_low_bits
        return point