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

def main(args):
    
    targetXhat = args.targetX_hat
    predXhat = args.predX_hat
    outdir = args.outdir
    
    X_hat = pd.read_csv(targetXhat)
    X_hat_pred = pd.read_csv(predXhat)
    
    prediction_Corr=  check_correlation(X_hat_pred,X_hat)
    
    with open(outdir+"/cellWisePearsonCorr.txt", "w") as file:
        file.write(str(prediction_Corr))
        
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Optional app description');
    parser.add_argument('--targetX_hat', type=str, help='test modality X_hat"', required=True);
    parser.add_argument('--predX_hat', type=str, help='inferred modality X_hat"', required=True);
    parser.add_argument('--outdir', type=str, help='output directory to save results', default='./results')
    args = parser.parse_args();
    main(args);