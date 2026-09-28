import numpy as np

def get_clusters(w, charge):
    
    LX, LY = w.shape
    clusters = []
    
    def neighbours(i, j):
        return [
            ((i-1)%LX, j), (i, (j+1)%LY), 
            ((i+1)%LX, j), (i, (j-1)%LY)
        ]
    
    def is_boundary(idx):
        pass
    
    def get_cluster(idx):
        cluster, queue = [idx], [idx]
        while queue:
            i, j = queue.pop(0)
            # ** grow cluster **
            for nbr in neighbours(i, j):
                if ( w[nbr]==charge and nbr not in cluster ):
                    cluster.append(nbr)
                    queue.append(nbr)
        return cluster
    
    for idx in np.ndindex(w.shape):
        if w[idx] == charge:
            clusters.append(get_cluster(idx))
    
    return clusters

