from scipy import sparse
import pandas as pd
import numpy as np
import anndata
from scHCLpy import main
import tqdm
from sklearn.metrics import pairwise_distances
from scHCLpy import reference_hcl
import obonet

def construct_reference(adata, celltype_field):
    # TODO
    """
    given an annotated dataset of celltypes, construct the reference:
    - for each celltype, subsample 100 cells
    - form the average expression, normalized to 100k umis
    - round
    - between all celltypes: calculate top20 DE genes
    """
    n_cells = 100
    n_repeats = 3
    topN_genes = 20

    adata_ref_x = []
    adata_ref_obs = []
    for ct in adata.obs[celltype_field].unique():
        A = adata[adata.obs[celltype_field]==ct]
        for i in range(n_repeats):
            ix = np.random.choice(len(A), size=np.minimum(n_cells, len(A)), replace=False)
            X = A[ix].X.toarray()
            X = X.sum(0)
            X = 1e5 * X/X.sum()
            X = np.log(1+X)
            adata_ref_x.append(X)
            adata_ref_obs.append({'celltype': ct, 'repeat': i})


class HCLReference_adata():
    def __init__(self, reference_df):
        self.reference_df = reference_df

    def query(self, adata, n_cores=1):
        """
        query all cells within the AnnData again this reference
        :params adata:
        :params n_cores:
        :returns: two pd.DataFrames, the first with some basic cell type annotation (best match)
                  the second dataFrame has the top 10 matches
        """
        transformed_adata = process_adata(adata, self.reference_df)
        scHCL_df, scHCL_df_extended_Celltypes = call_celltypes(transformed_adata, self.reference_df, n_cores)
        return scHCL_df, scHCL_df_extended_Celltypes

    def get_celltypes(self):
        """
        Returns a list of all cell types in the reference
        """
        return self.reference_df.index.values

    def get_genes(self):
        """
        Returns a list of all genes in the reference
        """
        return self.reference_df.columns.values


def process_adata(adata, reference_df):
    """
    matches the set of genes between query adata and reference to a common
    subset. Also transforms the query data into normalized log counts
    (required to compare against the reference, which had the same transformation).

    :returns: AnnData object with a subset of genes, .X is normalized log counts
    """

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


def calc_correlation_in_batches(adata, reference_df, BATCHSIZE=1000, n_cores=1):
    """
    unfortunately, sklearn and scipy dont have correlation distance for sparse matrices.
    but we can just batch the query adata into several sets of cells, calculate
    their correlation with all references and yield the result.
    Hence we never have to load the entire adata.X into memory

    :params adata: query AnnData. Each cell is matched against the reference
    :params reference_df: DataFrame of reference expression values
    :params BATCHSIZE: number of cells processed simultaniously, anything around 1000 should fit into mem easily
    :params n_cores: CPUs used to do the correlation calculation
    """
    BATCHSIZE = 1000  # cells being done at once

    # to get rid of these
    """
    anndata.py:1094: FutureWarning: is_categorical is deprecated and will be removed in a future version.  Use is_categorical_dtype instead
    """
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter(action='ignore', category=FutureWarning)

        adata_sorted = adata[:, reference_df.columns]
        for i in tqdm.tqdm(range(0, len(adata_sorted), BATCHSIZE)):
            _tmp_adata = adata_sorted[i:i+BATCHSIZE]

            X_query = _tmp_adata.X.toarray()
            # C = 1 - cdist(X_query, reference_df, 'correlation')
            C = 1 - pairwise_distances(X_query, reference_df, metric='correlation', n_jobs=n_cores)
            C = pd.DataFrame(C, index=_tmp_adata.obs.index, columns=reference_df.index)
            yield C


def call_celltypes(transformed_adata, ref_df, n_cores):
    """
    given transformed query data, for each cell in there, grind out the most
    correlated reference celltypes
    """
    scHCL_df = []
    scHCL_df_extended_Celltypes = []
    for C_batch in calc_correlation_in_batches(transformed_adata, ref_df, n_cores=n_cores):
        scHCL_batch = main.call_celltypes(C_batch)
        scHCL_batch_extended = main.call_celltypes2(C_batch, n_best=10)
        scHCL_df.append(scHCL_batch)
        scHCL_df_extended_Celltypes.append(scHCL_batch_extended)

    scHCL_df = pd.concat(scHCL_df)
    scHCL_df_extended_Celltypes = pd.concat(scHCL_df_extended_Celltypes)

    return scHCL_df, scHCL_df_extended_Celltypes


def create_ontology_id_to_name_mapping():
    graph = obonet.read_obo('http://purl.obolibrary.org/obo/cl/cl-basic.obo')
    id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data=True)}
    return id_to_name

def scHCL_adata(adata, verbose=False, n_cores=1, n_min=10):
    """
    previous main function
    """

    ref_df = reference_hcl.load_HCL_reference()
    transformed_adata = process_adata(adata, ref_df)

    scHCL_df, scHCL_df_extended_Celltypes = call_celltypes(transformed_adata, ref_df, n_cores)
    scHCL_df = annotate_refined(scHCL_df, n_min)

    # cell ontology
    id_to_name = create_ontology_id_to_name_mapping()
    scHCL_df['CLid'] = scHCL_df['hcl_refined'].apply(lambda x: reference_hcl.refined_celltypes_to_cell_ontology[x] if x in reference_hcl.refined_celltypes_to_cell_ontology else 'unknown')
    scHCL_df['CL_name'] = scHCL_df['CLid'].apply(lambda x: id_to_name[x] if x in id_to_name else 'unknown')

    return scHCL_df, scHCL_df_extended_Celltypes


def annotate_refined(scHCL_df, n_min: int):
    """
    The celltypes in the original scHCL are a bit nasty, with nonstandard names
    and the same celltype being annotated with 3 diferent strings/names
    We simplify the nomenclature here, and also get rid of "rare" celltypes
    that are probably mistakes in annotation.
    """
    CORRELATION_CUTOFF = 0.25 # any cell with correlation less than this will be classified as other

    scHCL_df['hcl_refined'] = scHCL_df['hcl_celltype'].apply(reference_hcl.celltype_rename)

    # filter away the infrequent annotations
    common_cts = scHCL_df['hcl_refined'].value_counts()
    common_cts = common_cts[common_cts > n_min].index
    scHCL_df['hcl_refined'] = scHCL_df['hcl_refined'].apply(lambda x: x if x in common_cts else 'low frequency')

    # also mark anything thats got shady correlations as "other" since it doesnt
    # really match any reference celltypes
    scHCL_df.loc[scHCL_df['hcl_score'] < CORRELATION_CUTOFF, 'hcl_refined'] = 'low correlation'

    return scHCL_df
