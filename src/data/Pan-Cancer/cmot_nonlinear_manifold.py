# import maninetcluster here
import numpy as np
import pandas as pd
from sklearn import preprocessing

from ManinetCluster_python.alignment import manifold_nonlinear
from ManinetCluster_python.correspondence import Correspondence
from ManinetCluster_python.distance import SquaredL2
from ManinetCluster_python.neighborhood import neighbor_graph

class nonlinear_Manifold:
    def __init__(self,file1,file2,d,m,W=None,neighbors = 5):
        self.file1 = file1
        self.file2 = file2
        self.d = d
        self.mu = m
        self.W = W
        self.neighbors = neighbors
        
    def run(self,as_file=True):
        
        if as_file:
            X = pd.read_csv(self.file1)
    #         X.index = X["Unnamed: 0_x"]
    #         X = X.drop(["Unnamed: 0_x"], axis=1)
            #X = X.T


            Y = pd.read_csv(self.file2)
            #Y.index = Y["Unnamed: 0"]
            #Y = Y.drop(["Unnamed: 0"], axis=1)
            #Y=Y.T
        
            X = X.to_numpy()
            Y = Y.to_numpy()
        else:
            X = self.file1.to_numpy()
            Y = self.file2.to_numpy()

        nx = len(X)
        ny = len(Y)
        d = self.d
        
        
        
#         actual = generate_correlation_map(X, Y)
#         actual = np.zeros([nx, ny], dtype = int)
        
#         if self.percent_align:
#             for i in self.percent_align:
#                 actual[i][i] = 1
#             corr = Correspondence(matrix=actual)
#         else:
#             corr = Correspondence(matrix=actual)
        corr = Correspondence(matrix=self.W)
        X_normalized = preprocessing.normalize(X, norm='l2')
        Y_normalized = preprocessing.normalize(Y, norm='l2')
        

        
        Wx = neighbor_graph(X_normalized, k=self.neighbors)
        Wy = neighbor_graph(Y_normalized, k=self.neighbors)
        
        Xnew, Ynew = manifold_nonlinear(X, Y, corr, d, Wx, Wy,mu=self.mu,eps=1e-8)
        Ynew_df = pd.DataFrame(Ynew)
        Xnew_df = pd.DataFrame(Xnew)
        
#         Ynew_df.to_csv("nonlinear_manifold_projected_points_mod2_GE_corr_knn_"+str(neighbors)+"_"+str(Y.shape[1])+"_d"+str(d)+".csv", index=False)
#         pd.DataFrame(Xnew).to_csv("nonlinear_manifold_projected_points_mod1_SNP_corr_knn_"+str(neighbors)+"_"+str(X.shape[1])+"_d"+str(d)+".csv", index=False) 
        
        #Ynew_df.to_csv("nonlinear_manifold_projected_points_mod2_"+str(Y.shape[1])+"_d"+str(d)+".csv", index=False)
        #pd.DataFrame(Xnew).to_csv("nonlinear_manifold_projected_points_mod1_"+str(X.shape[1])+"_d"+str(d)+".csv", index=False)
        return Xnew_df,Ynew_df
        
def generate_correlation_map(x, y):
    """Correlate each n with each m.

    Parameters
    ----------
    x : np.array
      Shape N X T.

    y : np.array
      Shape M X T.

    Returns
    -------
    np.array
      N X M array in which each element is a correlation coefficient.

    """
    mu_x = x.mean(1)
    mu_y = y.mean(1)
    n = x.shape[1]
    if n != y.shape[1]:
        raise ValueError('x and y must ' +
                         'have the same number of timepoints.')
    s_x = x.std(1, ddof=n - 1)
    s_y = y.std(1, ddof=n - 1)
    cov = np.dot(x,
                 y.T) - n * np.dot(mu_x[:, np.newaxis],
                                  mu_y[np.newaxis, :])
    deno = np.dot(s_x[:, np.newaxis], s_y[np.newaxis, :])
    final = cov / np.dot(s_x[:, np.newaxis], s_y[np.newaxis, :])
    final[np.isnan(final)] = 0    

    return final
#     return cov / np.dot(s_x[:, np.newaxis], s_y[np.newaxis, :])
        
        
