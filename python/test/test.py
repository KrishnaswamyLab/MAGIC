#!/usr/bin/env python


from __future__ import print_function, division, absolute_import
import matplotlib as mpl
mpl.use("agg")
import magic
import numpy as np
import scprep
try:
    import anndata
except (ImportError, SyntaxError):
    # anndata not installed
    pass


def test_scdata():
    scdata = scprep.io.load_csv("../data/test_data.csv", cell_names=False)
    scdata = scprep.filter.filter_empty_cells(scdata)
    scdata = scprep.filter.filter_empty_genes(scdata)
    scdata_norm = scprep.normalize.library_size_normalize(scdata)
    scdata_norm = scprep.transform.sqrt(scdata_norm)
    assert scdata.shape == scdata_norm.shape
    np.random.seed(42)
    magic_op = magic.MAGIC(t='auto', a=20, k=10)
    str_gene_magic = magic_op.fit_transform(
        scdata_norm, genes=['VIM', 'ZEB1'])
    int_gene_magic = magic_op.fit_transform(
        scdata_norm, genes=[-2, -1])
    assert str_gene_magic.shape[0] == scdata_norm.shape[0]
    assert np.all(str_gene_magic == int_gene_magic)
    pca_magic = magic_op.fit_transform(
        scdata_norm, genes="pca_only")
    assert pca_magic.shape[0] == scdata_norm.shape[0]
    assert pca_magic.shape[1] == magic_op.n_pca
    magic_all_genes = magic_op.fit_transform(scdata_norm,
                                             genes="all_genes")
    assert scdata_norm.shape == magic_all_genes.shape
    dremi = magic_op.knnDREMI("VIM", "ZEB1", plot=True)
    np.testing.assert_allclose(dremi, 1.573619, atol=0.0000005)


def test_anndata():
    try:
        anndata
    except NameError:
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
