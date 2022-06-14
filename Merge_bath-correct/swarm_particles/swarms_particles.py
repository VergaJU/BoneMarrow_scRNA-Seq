# Import modules
from zlib import Z_BEST_COMPRESSION
import numpy as np
import random
import scanpy as sc

# Import PySwarms
import pyswarms as ps
from pyswarms.utils.functions import single_obj as fx


class swarm():
    def __init__(self, A) -> None:
        """
        :param A: The initial matrix with the cells and gene expressions

        """
        self.A = A
        rows = A.shape[0]
        cols = A.shape[1]
        self.w = random.random()
        self.B = np.random.rand(2,cols)
         
        pass
    def evaluate(self, B):
        """
            :param B: Random generated matrix B with size m x 2
            :return The evaluation of this matrix given the zScore
        """
        C = np.matmul(B,self.A)
        return self.zScore(C)

    def zScore(self, C):
        """
            :Params "Corrected" matrix
            :return: The zScore value
        """
        from numpy import linalg as LA
        x = LA.eigvals(C)
        return x[0]
    def run_optimizer(self):
        options = {'B':self.B, 'w':self.w}
        options = {'c1': 0.5, 'c2': 0.3, 'w':0.9, 'k': 2, 'p': 2}
        optimizer = ps.single.LocalBestPSO(n_particles=10, dimensions=2, options=options)
        cost, pos = optimizer.optimize(self.evaluate, iters=1000)


def pre_processing(file = "./SMM.h5ad"):
    adata = sc.read(file)
    # filter genes expressed in less than 3 cells
    sc.pp.filter_genes(adata, min_cells = 3)
    # normalize data
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.subsample(adata, n_obs=10)
    # new matrix used to be modified by birds
    adata_norm = adata.X
    return adata_norm

if __name__ == "__main__":
    A = pre_processing()
    birds = swarm(A)
    birds.run_optimizer()

