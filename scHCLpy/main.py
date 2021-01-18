import pandas as pd
import numpy as np
import scanpy as sc
# import seaborn as sns
from scipy.spatial.distance import cdist
import pathlib


# class HCL():
#     def __init__(self):
#         self.reference_df = load_reference()

DATADIR = pathlib.Path(__file__).parent / 'data'

def test_case():
    """
    the test case they also provide on their website with the lung cells
    """
    df_lung = pd.read_csv(DATADIR / 'scHCL_lung.expr.csv.gz', sep=',', index_col=0).T
    # need to convert to upper case gene names (as in the reference!)
    df_lung.columns = df_lung.columns.map(lambda x: x.upper())

"""
df_bio_synonym = biomart_mapping.biomart_query_all(verbose=False, extra_fields=['external_synonym'], force_download=True)
df_bio_synonym.query('external_synonym=="IL8"')
synonym_2symbol = {row['external_synonym']: row['hgnc_symbol'] for _,row in df_bio_synonym[['hgnc_symbol', 'external_synonym']].drop_duplicates().iterrows()}
renames = {_:synonym_2symbol[_] for _ in not_found_genes if _ in synonym_2symbol}
"""

# there's a few genes that need to be renamed from the reference data
manual_renames = {
    'WISP1': 'CCN4',
    'RNR1': 'MT-RNR1',
    'COX1': 'MT-CO1',
    'COX2': 'MT-CO2',
    'NOV': 'CCN3',
    'AGPAT9': 'GPAT3',
}
renames = {
    'WISP2': 'CCN5',
    'NOV': 'PLXNA1',
    'PRKCDBP': 'CAVIN3',
    'AC079767.4': 'LINC01857',
    'ACRC': 'GCNA',
    'IKBKAP': 'ELP1',
    'CXorf36': 'DIPK2B',
    'FAM150B': 'ALKAL2',
    'CTGF': 'CCN2',
    'HCP5': 'CYCSP5',
    'C1orf186': 'RHEX',
    'AIM1': 'CRYBG1',
    'SDPR': 'CAVIN2',
    'FAM26F': 'CALHM6',
    'LINC01272': 'PELATON',
    'AMICA1': 'JAML',
    'KIAA0125': 'FAM30A',
    'GNB2L1': 'RACK1',
    'HMHA1': 'ARHGAP45',
    'C1orf228': 'ARMH1',
    'ADRBK1': 'GRK2',
    'C10orf10': 'DEPP1',
    'NRD1': 'NRDC',
    'C10orf11': 'LRMDA',
    'RNF219-AS1': 'OBI1-AS1',
    'HMP19': 'NSG2',
    'FAM101B': 'RFLNB',
    'C1orf64': 'SRARP',
    'LINC00599': 'MIR124-1HG',
    'LINC00936': 'ATP2B1-AS1',
    'MUM1L1': 'PWWP3B',
    'KIRREL': 'KIRREL1',
    'LINC00152': 'CYTOR',
    'LINC00675': 'TMEM238L',
    'LINC01330': 'LRRC77P',
    'IL10RB-AS1': 'IL10RB-DT',
    'MKL2': 'MRTFB',
    'TMEM27': 'CLTRN',
    'C10orf128': 'TMEM273',
    'FAM134B': 'RETREG1',
    'KIAA1462': 'JCAD',
    'C8orf4': 'TCIM',
    'APCDD1L-AS1': 'APCDD1L-DT',
    'C2orf82': 'SNORC',
    'FAM212B': 'INKA2',
    'C6orf25': 'MPIG6B',
    'AC079776.2': 'LINC02572',
    'PPP1R2P9': 'PPP1R2C',
    'MLTK': 'MAP3K20',
    'C17orf76-AS1': 'SNHG29',
    'IL8': 'CXCL8',
    'PRMT10': 'PRMT9',
    'C19orf10': 'MYDGF',
    'TMEM66': 'SARAF',
    'PTPLAD2': 'HACD4',
    'NOTCH2NL': 'NOTCH2NLA',
    'C19orf59': 'MCEMP1',
    'GPR116': 'ADGRF5',
    'ELTD1': 'ADGRL4',
    'PPAP2A': 'PLPP1',
    'PPAP2B': 'PLPP3',
    'FAM132B': 'ERFE',
    'CD97': 'ADGRE5',
    'FIGF': 'VEGFD',
    'C1orf63': 'RSRP1',
    'WBP5': 'TCEAL9',
    'MRP63': 'MRPL57',
    'SGOL1': 'SGO1',
    'LSMD1': 'NAA38',
    'SMEK2': 'PPP4R3B',
    'IGJ': 'JCHAIN',
    'C9orf89': 'CARD19',
    'DARC': 'ACKR1',
    'TENC1': 'TNS2',
    'LINC00657': 'NORAD',
    'CYR61': 'CCN1',
    'FAM65B': 'RIPOR2',
    'FAM159B': 'SHISAL2B',
    'NKX2.2': 'NKX2-2',
    'KIAA1244': 'ARFGEF3',
    'GPR56': 'ADGRG1',
    'ERO1LB': 'ERO1B',
    'HNRPDL': 'HNRNPDL',
    'KIAA1377': 'CEP126',
    'CTSL1': 'CTSL',
    'MIR143HG': 'CARMN',
    'FAIM3': 'FCMR',
    'PPAP2C': 'PLPP2',
    'LEPRE1': 'P3H1',
    'CCDC109B': 'MCUB',
    'RP11-396F22.1': 'CPNE8-AS1',
    'DPCR1': 'MUCL3',
    'GIF': 'CBLIF',
    'FAM101A': 'RFLNA',
    'C11orf70': 'CFAP300',
    'COX3': 'MT-CO3',
    'ND1': 'IVNS1ABP',
    'ATP6': 'MT-ATP6',
    'ND2': 'MT-ND2',
    'ND3': 'MT-ND3',
    'ND4': 'MT-ND4',
    'RNR1': 'NR4A2',
    'COX1': 'PTGS1',
    'COX2': 'PTGS2',
    'CYTB': 'MT-CYB',
    'ND6': 'MT-ND6',
    'ATP8': 'MT-ATP8',
    'ND5': 'MT-ND5',
    'PVRL2': 'NECTIN2',
    'GRAMD2': 'GRAMD2A',
    'KIAA1161': 'MYORG',
    'C10orf25': 'ZNF22-AS1',
    'RNU12': 'RNU12-2P',
    'AC091814.2': 'LINC02470',
    'ZNF788': 'ZNF788P',
    'CYBASC3': 'CYB561A3',
    'RPS17L': 'RPS17',
    'SELS': 'SELENOS',
    'HLA-DQB': 'HLA-DQB1',
    'MTND5': 'MT-ND5',
    'TCRBV3S1': 'TRBV28',
    'GUCY1B3': 'GUCY1B1',
    'FAM69B': 'DIPK1B',
    'AC092159.2': 'TMEM18-DT',
    'LINC01158': 'PANTR1',
    'RP11-480I12.5': 'TUBA5P',
    'CTC-534A2.2': 'SHLD3',
    'FAM64A': 'PIMREG',
    'NGFRAP1': 'BEX3',
    'LPPR1': 'PLPPR1',
    'LECT1': 'CNMD',
    'C3orf79': 'LINC02877',
    'APITD1-CORT': 'CENPS-CORT',
    'MRVI1-AS1': 'IRAG1-AS1',
    'C3orf58': 'DIPK2A',
    'KIAA0020': 'PUM3',
    'NHP2L1': 'SNU13',
    'FAM178A': 'SLF2',
    'ZMYM6NB': 'TMEM35B',
    'AGPAT9': 'LPCAT1',
    'C1orf168': 'FYB2',
    'C10orf32': 'BORCS7',
    'FLJ22763': 'C3orf85',
    'SSSCA1-AS1': 'ZNRD2-AS1',
    'C14orf105': 'CCDC198',
    'KIAA1456': 'TRMT9B',
    'FAM60A': 'SINHCAF',
    'GUCY1A3': 'GUCY1A1',
    'KIAA1024': 'MINAR1',
    'RP11-532F12.5': 'SPINT1-AS1',
    'FAM65C': 'RIPOR3',
    'KIAA0101': 'PCLAF',
    'LINC01314': 'CTXND1',
    'FAM46B': 'TENT5B',
    'MIR4435-1HG': 'MIR4435-2HG',
    'GPR98': 'ADGRV1',
    'TMEM8C': 'MYMK',
    'KIAA1524': 'CIP2A',
    'DOPEY2': 'DOP1B',
    'LINC00704': 'MANCR',
    'DFNB59': 'PJVK',
    'SNHG23': 'MEG8',
    'NCK1-AS1': 'NCK1-DT',
    'BZRAP1': 'TSPOAP1',
    'ADCK3': 'COQ8A',
    'RQCD1': 'CNOT9',
    'CBWD7': 'CBWD6',
    'HRSP12': 'RIDA',
    'CCDC23': 'SVBP',
    'ERO1L': 'ERO1A',
    'C1orf86': 'FAAP20',
    'PTPLAD1': 'HACD3',
    'UTP11L': 'UTP11',
    'ACN9': 'SDHAF3',
    'PTPLA': 'HACD1',
    'ND4L': 'MT-ND4L',
    'C12orf39': 'SPX',
    'CGB': 'CGB3',
    'FAM195A': 'MCRIP2',
    'TUBB2C': 'TUBB4B',
    'C6orf221': 'KHDC3L',
    'TUBB4Q': 'TUBB7P',
    'PLAC1L': 'OOSP2',
    'C7orf70': 'FAM220A',
    'C11orf73': 'HIKESHI',
    'MYEOV2': 'COPS9',
    'C5orf48': 'TEX43',
    'LINC00684': 'FAM236A',
    'TEX40': 'CATSPERZ',
    'SPTY2D1-AS1': 'SPTY2D1OS'
}
HLA_renames = {
    'HLA.A': 'HLA-A',
    'HLA.B': 'HLA-B',
    'HLA.C': 'HLA-C',
    'HLA.DMB': 'HLA-DMB',
    'HLA.DPA1': 'HLA-DPA1',
    'HLA.DPB1': 'HLA-DPB1',
    'HLA.DQA1': 'HLA-DQA1',
    'HLA.DQB1': 'HLA-DQB1',
    'HLA.DRA': 'HLA-DRA',
    'HLA.DRB1': 'HLA-DRB1',
    'HLA.DRB5': 'HLA-DRB5',
    'HLA.E': 'HLA-E',
    'HLA.F': 'HLA-F',
}

# for some symbols, its just not clear what the map to because ambiguities
unkown_renames = ['RNR1', 'COX1', 'COX2']


def _load_reference_raw():
    REF_FILE = DATADIR / 'scHCL_ref.expr.csv.gz'
    # print(REF_FILE)
    assert REF_FILE.exists()
    df_ref = pd.read_csv(REF_FILE, sep=',', index_col=0).T
    return df_ref


def load_reference():
    """
    loading the refernce data from the HCL paper. This is a dataFrame where each
    row corresponds to a particular cell/celltype and each column is a particular gene.
    Note that the genes are not exhaustive, as only genes that are somehow differential
    in at least on celltype will appear.

    Note that a log transform is applied (consistent with the query processing)
    """
    df_ref = _load_reference_raw()

    # the renaming is a little iffy!
    # drop the ones we cant map
    df_ref.drop(unkown_renames, axis=1, inplace=True)
    _tmp_renames = renames.copy()
    _tmp_renames.update(HLA_renames)
    _tmp_renames.update(manual_renames)  # this should overwrite some of the wrong automatic renames!

    df_ref = df_ref.rename(_tmp_renames, axis=1)

    # now some of genes might be present in mutiple columns (one col from the correct name, one from the rename)
    # so lets just add up
    cols = []
    all_genes = df_ref.columns.unique()
    for gene in all_genes:
        series = df_ref[[gene]].sum(1)  # [[]] needed: if just single occurance [genename] would be a series that we cant sum across
        series.name = gene
        cols.append(series)
    df_ref = pd.DataFrame(cols).T

    # log transform the ref
    df_ref = np.log(df_ref+1)

    return df_ref


def adata_to_df(adata, use_raw):
    query_df = pd.DataFrame(adata.raw.X.A if use_raw else adata.X.A,
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


def scHCL(query_df, verbose=False):

    ref_df = load_reference()
    transformed_query_df = process_query_df(query_df, ref_df)

    if verbose:
        print(f'Reference:\t {ref_df.shape[0]} cells, {ref_df.shape[1]} genes')
        print(f'Query:\t {query_df.shape[0]} cells, {query_df.shape[1]} genes')

    # calculate the Pearson correlation between each query and reference cell
    if verbose:
        print('Calculating Pearson correlation matrix')
    C = 1 - cdist(transformed_query_df, ref_df, 'correlation')  # 1- since cdist calculates the correlation distance, i.e. if perfectly correlated cdist=0
    C = pd.DataFrame(C, index=transformed_query_df.index, columns=ref_df.index)

    # find clostest matching reference type
    if verbose:
        print('Finding closest ref celltypes')

    ix_closest_ref_cell = np.argmax(C.values, axis=1)
    name_closest_ref_cell = C.columns.values[ix_closest_ref_cell]
    score = C.values[np.arange(len(ix_closest_ref_cell)), ix_closest_ref_cell]
    hcl_celltype_call = pd.Series(name_closest_ref_cell, index=C.index, name='hcl_celltype')
    hcl_celltype_score = pd.Series(score, index=C.index, name='hcl_score')

    scHCL_df = pd.DataFrame({'hcl_score': hcl_celltype_score,
                             'hcl_celltype': hcl_celltype_call})
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
