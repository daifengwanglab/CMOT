import os, sys, argparse, time, random, math

import ot
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold, StratifiedKFold, train_test_split
from sklearn.cluster import KMeans, AgglomerativeClustering 
import scipy.cluster.hierarchy as sch
from scipy import stats
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import Normalizer
from sklearn.metrics import roc_curve, auc,mean_squared_error,explained_variance_score,accuracy_score,balanced_accuracy_score, f1_score

from cmot_nonlinear_manifold import *
from cmot_functions import *

runtimes = []

def main(args):
       
    start = time.time()
    sourceX = args.sourceX
    sourceY = args.sourceY
    targetY = args.targetY
    labels =args.labels
    K = args.K
    d = args.d
    W = args.W
    hc = args.hc
    reg_e =args.reg_e
    reg_cl = args.reg_cl
    topFeat = args.topFeat
    k = args.K
    outdir = args.outdir
    
    X = pd.read_csv(sourceX)
    Y = pd.read_csv(sourceY)
    Y_hat = pd.read_csv(targetY)
    if labels:
        labels = pd.read_csv(labels)
    W = np.load(W)
    
    os.system('mkdir -p '+ outdir)

    # STEP A)
    X_aligned,Y_aligned = align_samples(X,Y,W,K=K,d=d,mu=0.5)
    best_match = get_best_match(X_aligned,Y_aligned)
    X_NMABestMatch, X_BestMatch = find_best_match_df(best_match,X_aligned, X)
    
    # STEP B)
    hc = AgglomerativeClustering(n_clusters = hc, affinity = 'euclidean', linkage ='ward')
    hclabels=hc.fit_predict(Y)
    l = pd.DataFrame(hclabels, columns = ["HC"])

    # STEP C)
    Xs = Y.to_numpy()
    ys = l["HC"].to_numpy()
    Xt = Y_hat.to_numpy() 
    method = "lpl1_reg"

    # print ("Starting Optimal transport using:", method)
    transp_Xs = optimal_transport(Xs,Xt,reg_e, reg_cl,ys, method)

    # print ("Predicting phenotype ...")
    X_hat_pred = predict_phenotype(transp_Xs,Y_hat,X_BestMatch,ys,k,topFeat)
    
    print ("Completed prediction!")
    print ("Total runtime:",time.time()-start)
    runtimes.append(time.time()-start)
    
    X_hat_pred.to_csv(outdir+"/Norm_ModalityXhat.csv", index=False)

    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Optional app description');
    parser.add_argument('--sourceX', type=str, help='source modality X"', required=True);
    parser.add_argument('--sourceY', type=str, help='source modality Y"', required=True);
    parser.add_argument('--labels', type=str, help='prior knowledge of cells labels in source modality Y"', default = None)
    parser.add_argument('--targetY', type=str, help='target modality Y_hat"', required=True);
    parser.add_argument('--K', type=int, help='specify nearest neighbors for NMA (Step A)"', default = 5);
    parser.add_argument('--d', type=int, help='specify latent dibinary matrix specifying correspondence between cells of X and Y (Step A)"', default = 10);
    parser.add_argument('--W', type=str, help='binary matrix specifying correspondence between cells of modalities X and Y"', default = 1e00);
    parser.add_argument('--hc', type=int, help='specify number of clusters as cell label information (Step B)"', default = 2);
    parser.add_argument('--reg_e', type=float, help='entropy regularization of OT (Step B)"', default = 1e00);
    parser.add_argument('--reg_cl', type=float, help='label regularizatoin of OT (Step B)"', default = 1e00);
    parser.add_argument('--topFeat', type=int, help='top variable features of Y to use for K-Nearest Neighbors (Step C)"', default = 100);
    parser.add_argument('--k', type=int, help='k-nearest neighbors for cross-modality inference (Step C)"', default = 1e00);
    parser.add_argument('--outdir', type=str, help='output directory to save results', default='./results')
    
    args = parser.parse_args();
    main(args);