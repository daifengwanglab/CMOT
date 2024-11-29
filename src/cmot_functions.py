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

def align_samples(X,Y,W,K,d=5,mu=0.5, method="NMA"):
    random.seed(30)
    
    if method == "NMA":
        obj = nonlinear_Manifold(X, Y,d,mu,W,K)
        X_aligned,Y_aligned = obj.run(as_file=False)
    elif method == "Unioncom":
        d=5
        obj = Unioncom(X, Y,d,epoch_pd=700, epoch_DNN=900)
        X_aligned,Y_aligned = obj.run(as_file=False)
    return X_aligned,Y_aligned

def kmeans(X,Y,k=6):
    X_Y = pd.concat([X,Y]).to_numpy()

    kmeanModel = KMeans(n_clusters=6).fit(X)
    kmeanModel.fit(X_Y)
    clusterLabels = kmeanModel.labels_  
    return clusterLabels

def check_mapping_quality(ot_lpl1, transp_Xs_lpl1):
    
    coupling = ot_lpl1.coupling_
    cost=ot_lpl1.cost_
    prob_cost=coupling#np.multiply(coupling,cost)
    prob_cost = np.transpose(prob_cost)

    pca = PCA()
    components = pca.fit_transform(prob_cost)

    var_ratio = np.cumsum(pca.explained_variance_ratio_)
    idx = np.where(var_ratio<=0.95)[0][-1]+1

    clf = IsolationForest(random_state = 1, contamination= 'auto')
    preds = clf.fit_predict(components[:,:idx])
    
    poor_mapped_cell_idx=np.argwhere(preds==-1).reshape(-1)
    values=poor_mapped_cell_idx
    with open("./poorly_mapped_cells.txt", "w") as output:
        for row in values:
            output.write(str(row) + '\n')
  
    poor_mapped_cell_percent = (len(np.argwhere(preds==-1).reshape(-1))/prob_cost.shape[0])*100
    print ("Warning: Found {} poorly mapped cells from the target modaity. Users can remove these cells to \
    avoid incorrect imputation. Their cell ids are stored in poorly_mapped_cells.txt".format(poor_mapped_cell_percent))

def optimal_transport(Xs,Xt,reg_e, reg_cl,ys=None,method="lpl1_reg"):
   
    if method == "emd":
        ot_emd = ot.da.EMDTransport()
        ot_emd.fit(Xs=Xs, Xt=Xt)
        transp_Xs_emd = ot_emd.transform(Xs=Xs)        
        return transp_Xs_emd

    elif method == "sinkhorn":
        ot_sinkhorn = ot.da.SinkhornTransport(reg_e=1e-1)
        ot_sinkhorn.fit(Xs=Xs, Xt=Xt)
        transp_Xs_sinkhorn = ot_sinkhorn.transform(Xs=Xs)
        return transp_Xs_sinkhorn

    elif method == "lpl1_reg":
        ot_lpl1 = ot.da.SinkhornLpl1Transport(reg_e=reg_e, reg_cl=reg_cl)
        ot_lpl1.fit(Xs=Xs, ys=ys, Xt=Xt)
        transp_Xs_lpl1 = ot_lpl1.transform(Xs=Xs)
        check_mapping_quality(ot_lpl1, transp_Xs_lpl1)
        return transp_Xs_lpl1
       
    elif method == "emd_laplace":
        ot_emd_laplace = ot.da.EMDLaplaceTransport(reg_lap=10, reg_src=10)
        ot_emd_laplace.fit(Xs=Xs, Xt=Xt)
        transp_Xs_emd_laplace = ot_emd_laplace.transform(Xs=Xs)
        return transp_Xs_emd_laplace
        
    elif method == "l1l2_reg":
        ot_l1l2 = ot.da.SinkhornL1l2Transport(reg_e=1e-1, reg_cl=2e0, max_iter=20,
                                      verbose=True)
        ot_l1l2.fit(Xs=Xs, ys=ys, Xt=Xt)
        transp_Xs_l1l2 = ot_l1l2.transform(Xs=Xs)
        return transp_Xs_l1l2
    
def get_nn(X,ys,k,components,targetcomponent):
    knn = NearestNeighbors(n_neighbors=k)
    knn.fit(components)
    neighDist, neigh = knn.kneighbors(targetcomponent, return_distance=True)
    nn = X.iloc[neigh[0]]
    nnL = ys[neigh[0]]
    
    return nn, nnL, neighDist,neigh

def get_corr_with_n(nn,targetGE):
    for idx,row in nn.iterrows():
        b = np.array(row).reshape(1,-1)
        corr = stats.pearsonr(targetGE[0][:150], b[0][:150])
        nnCorr.append(corr[0])
        
    return nnCorr

def get_phenotype(X,ys,transp_Xs,targetSNP,k,topFeat,pca=False):
    
#     targetSNP = Xt[targetIdx][:topFeat].reshape(1,-1)
    Xp = transp_Xs[:,:topFeat]

    if pca:
        pca = PCA(n_components=4)
        components = pca.fit_transform(Xp)
        targetcomponent = pca.transform(targetSNP)
    else:
        components = Xp
        targetcomponent = targetSNP
    
    nn,nnL,neighDist,neighIdx = get_nn(X,ys,k,components,targetcomponent)

    neighDist = neighDist[0]
    for i in range(len(neighDist)):
        neighDist[i] = math.exp(-neighDist[i])
    
    nn = nn.multiply(neighDist, axis=0)
    
    averageGE = nn.mean()
    return averageGE

def predict_phenotype(transp_Xs,Y_hat,X,ys,k,topFeat):

    predGEdf = pd.DataFrame(columns=X.columns)
    for i in range(len(Y_hat)):
        Y_hat_i = Y_hat.to_numpy()[i][:topFeat].reshape(1,-1)
        predGE = get_phenotype(X,ys,transp_Xs,Y_hat_i,k,topFeat)
        predGEdf = pd.concat([predGEdf,pd.DataFrame(predGE).T])
        
    predGEdf=predGEdf.T
    for col in predGEdf.columns:
        predGEdf[col] = NormalizeData(predGEdf[col])
    predGEdf = predGEdf.T
    return predGEdf

def check_correlation(predicted_p,X_hat):
    predCorr = []
    for i in range(len(predicted_p)):
        a = predicted_p.iloc[i][:]
        b = X_hat.iloc[i][:]
        corr = stats.pearsonr(a,b)
        predCorr.append(corr[0])
    return predCorr

def NormalizeData(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))

def get_best_match(geNMA,snpNMA):
    best_match =[]
    for i in range(len(snpNMA)):
        X = geNMA.to_numpy()
        y_i = snpNMA.iloc[i].to_numpy().reshape(1,-1)
        knn = NearestNeighbors(n_neighbors=1).fit(X)
        dist,indices = knn.kneighbors(y_i)
        best_match.append(indices[0][0])
    return best_match


def find_best_match_df(best_match,geNMA, geXTrain):
    geNMABestMatch = geNMA.copy()
    geNMABestMatch = geNMABestMatch.reset_index()
    geNMABestMatch = geNMABestMatch.reindex(best_match)
    geNMABestMatch = geNMABestMatch.drop(["index"], axis=1)
    geNMABestMatch = geNMABestMatch.reset_index(drop=True)
    
    geXTrainBestMatch = geXTrain.copy()
    geXTrainBestMatch = geXTrainBestMatch.reset_index()
    geXTrainBestMatch = geXTrainBestMatch.reindex(best_match)
    geXTrainBestMatch = geXTrainBestMatch.drop(["index"], axis=1)
    geXTrainBestMatch = geXTrainBestMatch.reset_index(drop=True)
    
    return geNMABestMatch, geXTrainBestMatch

def cmot_parameter_search(X,Y,X_aligned,Y_aligned,labels): 
    
    skf = StratifiedKFold(n_splits=2)
    X=X.reset_index(drop=True)
    Y=Y.reset_index(drop=True)
    X_aligned=X_aligned.reset_index(drop=True)
    Y_aligned=Y_aligned.reset_index(drop=True)
    labels=labels.reset_index(drop=True)
    
    for train_index, test_index in skf.split(X, labels):

        X_ps, X_hat_ps = X.loc[train_index], X.loc[test_index]
        Y_ps, Y_hat_ps = Y.loc[train_index], Y.loc[test_index]

        X_aligned_ps, X_hat_aligned_ps = X_aligned.loc[train_index], X_aligned.loc[test_index]
        Y_aligned_ps, Y_hat_aligned_ps = Y_aligned.loc[train_index], Y_aligned.loc[test_index]
        
        best_match = get_best_match(X_aligned_ps,Y_aligned_ps)
        X_NMABestMatch, X_BestMatch = find_best_match_df(best_match,X_aligned_ps, X_ps)
        
        if labels is not None:
            labelYTrain_ps, labelYTest_ps = labels.loc[train_index], labels.loc[test_index]
            
        min_err=np.inf
        best_reg_e = None
        for reg_e in [10.0,100.0,200.0,500.0,1000.0,5000.0]:
            ot_lpl1 = ot.da.SinkhornLpl1Transport(reg_e=reg_e, reg_cl=1e00,log=True)
            ot_lpl1.fit(Xs=Y_ps.to_numpy(), ys=labelYTrain_ps["HC"].to_numpy(), Xt=Y_hat_ps.to_numpy())
            transp_Xs_lpl1 = ot_lpl1.transform(Xs=Y_ps.to_numpy())  
            if ot_lpl1.log_['err'][-1]<min_err:
                min_err=ot_lpl1.log_['err'][-1]
                best_reg_e = reg_e
                
        ot_lpl1 = ot.da.SinkhornLpl1Transport(reg_e=best_reg_e, reg_cl=1e00)
        ot_lpl1.fit(Xs=Y_ps.to_numpy(), ys=labelYTrain_ps["HC"].to_numpy(), Xt=Y_hat_ps.to_numpy())
        transp_Xs_lpl1 = ot_lpl1.transform(Xs=Y_ps.to_numpy())  
        
        k_vs_topfeat=pd.DataFrame()
        for k in list(range(20,min(Y_ps.shape[0], Y_ps.shape[1]),20)):
            for_each_k=[]
            for topFeat in list(range(20,Y_ps.shape[1],20)):
                
                X_hat_ps_pred = predict_phenotype(transp_Xs_lpl1,Y_hat_ps,X_BestMatch,labelYTrain_ps.to_numpy(),k,topFeat)

                predCorr = []
                for i in range(len(X_hat_ps_pred)):
                    a = X_hat_ps_pred.iloc[i][:]
                    b = X_hat_ps.iloc[i][:]
                    corr = stats.pearsonr(a,b)
                    corr_sp = stats.spearmanr(a,b)
                    predCorr.append(corr[0])

                for_each_k.append(np.mean(predCorr))

            k_vs_topfeat[str(k)] = for_each_k
        
        k_vs_topfeat.index=list(range(20,Y_ps.shape[1],20))

        a, b = k_vs_topfeat.stack().idxmin()
        best_top_feat = min(k_vs_topfeat.loc[[a], [b]].columns.to_list())
        best_k = min(k_vs_topfeat.loc[[a], [b]].index.to_list())
        
        return best_reg_e, best_top_feat, best_k
