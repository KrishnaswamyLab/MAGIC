#!/usr/bin/env python


from __future__ import print_function, division, absolute_import
import magic
import pandas as pd
import numpy as np
try:
    import anndata
except (ImportError, SyntaxError):
    # anndata not installed
    pass


def test_scdata():
    scdata = pd.read_csv("../data/test_data.csv")
    scdata_norm = magic.preprocessing.library_size_normalize(scdata)
    assert scdata.shape == scdata_norm.shape

    fast_magic_operator = magic.MAGIC(t='auto', a=20, k=10)

    str_gene_magic = fast_magic_operator.fit_transform(
        scdata_norm, genes=['VIM', 'ZEB1'])
    int_gene_magic = fast_magic_operator.fit_transform(
        scdata_norm, genes=[-2, -1])
    assert str_gene_magic.shape[0] == scdata_norm.shape[0]
    assert np.all(str_gene_magic == int_gene_magic)

    pca_magic = fast_magic_operator.fit_transform(
        scdata_norm, genes="pca_only")
    assert pca_magic.shape[0] == scdata_norm.shape[0]
    assert pca_magic.shape[1] == fast_magic_operator.n_pca

    fast_magic = fast_magic_operator.fit_transform(scdata_norm,
                                                   genes="all_genes")
    assert scdata_norm.shape == fast_magic.shape


def test_anndata():
    try:
        anndata
    except (ImportError, SyntaxError):
        # anndata not installed
        return
    scdata = anndata.read_csv("../data/test_data.csv")
    fast_magic_operator = magic.MAGIC(t='auto', a=None, k=10)
    sc_magic = fast_magic_operator.fit_transform(
        scdata, genes="all_genes")
    assert np.all(sc_magic.var_names == scdata.var_names)
    assert np.all(sc_magic.obs_names == scdata.obs_names)
    sc_magic = fast_magic_operator.fit_transform(
        scdata, genes=['VIM', 'ZEB1'])
    assert np.all(sc_magic.var_names.values == np.array(['VIM', 'ZEB1']))
    assert np.all(sc_magic.obs_names == scdata.obs_names)
