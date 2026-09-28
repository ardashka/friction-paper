from numpy import linalg

def _dist1(x1, x2, L):
    """ one-dimenaional distance. \n
    Args:
        x1: Integer in [0,L].
        x2: Integer in [0,L].
        L: Size of the interval [0,L].
    Returns:
        1D distance on the periodic interval [0,L].
    """
    return min((x1 - x2)%L, (x2 - x1)%L)


def euclidean(p1, p2, shape):
    """ Euclidean metric on a flat square torus. \n
    Args:
        p1: Tuple of the form (x1, x2, ..., xn).
        p2: Tuple of the form (y1, y2, ..., yn).
        shape: Size of the periodic domain (list).
    Returns:
        The Euclidean distance between p1 and p2 on the flat torus.
    """
    return linalg.norm([ _dist1(x1, x2, L) for x1, x2, L in zip(p1, p2, shape) ])


def taxicab(p1, p2, shape):
    """ Taxicab metric on a flat square torus. \n
    Args:
        p1: Tuple of the form (x1, x2, ..., xn).
        p2: Tuple of the form (y1, y2, ..., yn).
        shape: Size of the periodic domain (list).
    Returns:
        The taxicab distance between p1 and p2 on the flat torus.
    """
    return sum([ _dist1(x1, x2, L) for x1, x2, L in zip(p1, p2, shape) ])
