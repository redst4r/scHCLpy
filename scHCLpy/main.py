import pandas as pd
import numpy as np
# import seaborn as sns
from sklearn.metrics import pairwise_distances
from scHCLpy import reference_hcl

class HCLReference():
    def __init__(self, reference_df):
        self.reference_df = reference_df

    def query(self, query_df, n_cores=1):
        transformed_query = process_query_df(query_df, self.reference_df)
        C = calculate_correlation_matrix(transformed_query, self.reference_df, n_cores)
        scHCL_df = call_celltypes(C)
        return scHCL_df, C


# def test_case():
#     """
#     the test case they also provide on their website with the lung cells
#     """
#     df_lung = pd.read_csv(DATADIR / 'scHCL_lung.expr.csv.gz', sep=',', index_col=0).T
#     # need to convert to upper case gene names (as in the reference!)
#     df_lung.columns = df_lung.columns.map(lambda x: x.upper())

"""
df_bio_synonym = biomart_mapping.biomart_query_all(verbose=False, extra_fields=['external_synonym'], force_download=True)
df_bio_synonym.query('external_synonym=="IL8"')
synonym_2symbol = {row['external_synonym']: row['hgnc_symbol'] for _,row in df_bio_synonym[['hgnc_symbol', 'external_synonym']].drop_duplicates().iterrows()}
renames = {_:synonym_2symbol[_] for _ in not_found_genes if _ in synonym_2symbol}
"""


def adata_to_df(adata, use_raw):
    query_df = pd.DataFrame(adata.raw.X.toarray() if use_raw else adata.X.toarray(),
                            index=adata.obs.index,
                            columns=adata.raw.var.index if use_raw else adata.var.index)
    return query_df


def process_query_df(query_df, reference_df):
    """
    apply the same preprocessing steps to the query
    as has been done to the reference
    """

    test_df = []

    # find the genes from the reference in the query
    not_found_genes = []
    found_genes = []
    for gene in reference_df.columns:
        if gene in query_df.columns:
            found_genes.append(gene)
            test_df.append(query_df[gene])
        else:
            # if the gene is not measured in the query, just assume it's zero
            not_found_genes.append(gene)
            _t = pd.Series(0, index=query_df.index, name=gene)
            test_df.append(_t)

    test_df = pd.DataFrame(test_df).T
    print(f'{len(found_genes)} reference genes found, {len(not_found_genes)} not found!')

    # transforming: library size normalization and log+1 transform
    test_df = 1e5 * test_df / test_df.values.sum(1, keepdims=True)
    test_df = np.log(test_df+1)

    return test_df


def call_celltypes2(C, n_best):
    # nastyu, since argsort is ascending! firsr reverse, then take the topN
    ix_sorted = np.argsort(C.values, axis=1)[:,::-1][:, :n_best]
    name_closest_ref_cell = C.columns.values[ix_sorted]

    # the indexing gets a little funky, better do it one by one:
    # 1. top Hit
    # 2. 2nd top Hit
    # ....
    score = []
    for i in range(n_best):
        s = C.values[np.arange(len(ix_sorted)), ix_sorted[:,i]]
        score.append(s)
    score = np.vstack(score).T

    hcl_celltype_call = pd.DataFrame(name_closest_ref_cell, index=C.index, columns=[f'hcl_celltype_top{1+_}' for _ in range(n_best)])
    hcl_celltype_score = pd.DataFrame(score, index=C.index, columns=[f'hcl_score_top{1+_}' for _ in range(n_best)])
    scHCL_df = pd.concat([hcl_celltype_call, hcl_celltype_score], axis=1)
    return scHCL_df


def call_celltypes(C):
    ix_closest_ref_cell = np.argmax(C.values, axis=1)
    name_closest_ref_cell = C.columns.values[ix_closest_ref_cell]
    score = C.values[np.arange(len(ix_closest_ref_cell)), ix_closest_ref_cell]
    hcl_celltype_call = pd.Series(name_closest_ref_cell, index=C.index, name='hcl_celltype')
    hcl_celltype_score = pd.Series(score, index=C.index, name='hcl_score')
    scHCL_df = pd.DataFrame({'hcl_score': hcl_celltype_score,
                             'hcl_celltype': hcl_celltype_call})
    return scHCL_df


def calculate_correlation_matrix(transformed_query_df, ref_df, n_cores):
    assert np.all(transformed_query_df.columns == ref_df.columns)

    # 1- since cdist calculates the correlation distance, i.e. if perfectly correlated cdist=0
    # C = 1 - cdist(transformed_query_df, ref_df, 'correlation')

    C = 1 - pairwise_distances(transformed_query_df, ref_df, metric='correlation', n_jobs=n_cores)
    C = pd.DataFrame(C, index=transformed_query_df.index, columns=ref_df.index)

    return C


def scHCL(query_df, verbose=False, n_cores=1):

    ref_df = reference_hcl.load_HCL_reference()
    transformed_query_df = process_query_df(query_df, ref_df)

    if verbose:
        print(f'Reference:\t {ref_df.shape[0]} cells, {ref_df.shape[1]} genes')
        print(f'Query:\t {query_df.shape[0]} cells, {query_df.shape[1]} genes')

    # calculate the Pearson correlation between each query and reference cell
    if verbose:
        print('Calculating Pearson correlation matrix')

    C = calculate_correlation_matrix(transformed_query_df, ref_df, n_cores)

    # find clostest matching reference type
    if verbose:
        print('Finding closest ref celltypes')

    scHCL_df = call_celltypes(C)

    relevant_cts = get_closest_cells(C, n_closest=3)
    if verbose:
        print('clustermap')
    # sns.clustermap(C.loc[:, relevant_cts])

    return scHCL_df, C


def get_closest_cells(C, n_closest):
    """
    across all query cells get the closest referenc cell types (n_closest per query)
    and compile this list of present ref celltypes
    """
    ref_cells = set()
    for cellid, row in C.iterrows():
        ix_best = np.argsort(row.values)[::-1][:n_closest]
        names = C.columns.values[ix_best]
        ref_cells = ref_cells | set(names)  # add them to the reference cells
    return list(ref_cells)
