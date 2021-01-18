from scipy import sparse
from scipy.spatial.distance import cdist
import pandas as pd
import numpy as np
import anndata
from scHCLpy import main
import tqdm


def process_adata(adata, reference_df):

    adata_genes = set(adata.raw.var_names)
    refer_genes = set(reference_df.columns)
    shared_genes = sorted(list(adata_genes & refer_genes))
    adata_missing_genes = sorted(list(refer_genes - adata_genes))

    print(f'{len(shared_genes)} reference genes found, {len(adata_missing_genes)} not found!')

    var_shared = pd.Series(shared_genes)
    var_missing = pd.DataFrame(adata_missing_genes)
    newVar = pd.concat([var_shared, var_missing]).set_index(0)

    X_shared = adata.raw[:, shared_genes].X
    X_missing = sparse.csr_matrix((X_shared.shape[0], len(adata_missing_genes)))
    newX = sparse.hstack([X_shared, X_missing])

    assert newX.shape[1] == reference_df.shape[1]

    # transform the X

    # normalizing row is a bit cumbersome with sparse matrices!
    Z = np.array(newX.sum(1)).flatten()  # sum over columns
    scalingMat = sparse.diags(1e5 / (Z), offsets=0)
    X_scaled = scalingMat @ newX

    # log(1+x):
    # note that original zeros in the matrix stay zero
    X_scaled.data = np.log(1+X_scaled.data)

    new_adata = anndata.AnnData(X_scaled, var=newVar, obs=adata.obs.copy())
    return new_adata


def calc_correlation_in_batches(adata, reference_df, BATCHSIZE=1000):
    """
    unfortunately, sklearn and scipy dont have correlation distance for sparse matrices.
    but we can just batch the query adata into several sets of cells, calculate
    their correlation with all references and yield the result.
    Hence we never have to load the entire adata.X into memory
    :params BATCHSIZE: number of cells processed simultaniously, anything around 1000 should fit into mem easily
    """
    BATCHSIZE = 1000  # cells being done at once

    adata_sorted = adata[:, reference_df.columns]
    for i in tqdm.tqdm(range(0, len(adata_sorted), BATCHSIZE)):
        _tmp_adata = adata_sorted[i:i+BATCHSIZE]

        X_query = _tmp_adata.X.A
        C = 1 - cdist(X_query, reference_df, 'correlation')
        C = pd.DataFrame(C, index=_tmp_adata.obs.index, columns=reference_df.index)
        yield C


"""
from scHCLpy import main
reference_df = main.load_reference()
pd.read_csv(pathlib.Path(main.__file__).parent / 'data' / 'scHCL_ref.expr.csv.gz')
"""
def scHCL_adata(adata, verbose=False):
    ref_df = main.load_reference()
    transformed_adata = process_adata(adata, ref_df)

    scHCL_df = []
    for C_batch in calc_correlation_in_batches(transformed_adata, ref_df):
        scHCL_batch = main.call_celltypes(C_batch)
        scHCL_df.append(scHCL_batch)

    scHCL_df = pd.concat(scHCL_df)
    return scHCL_df
