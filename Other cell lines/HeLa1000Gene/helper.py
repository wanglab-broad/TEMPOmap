import numpy as np
import scanpy as sc


def separate_adata_by_sample(adata, sample_list, adata_layers=None):
    '''for separating AnnData that has all samples (with layers)
        where var is kept unchanged while obs & layers are separated
        - one or multiple samples
    '''
    X = adata.X[np.isin(adata.obs['sample'], sample_list)]
    obs = adata.obs.loc[np.isin(adata.obs['sample'], sample_list)]
    adata_new = sc.AnnData(X=X, obs=obs, var=adata.var)
    if adata_layers:
        for layer in adata_layers:
            adata_new.layers[layer] = adata.layers[layer][np.isin(adata.obs['sample'], sample_list)]
    return adata_new