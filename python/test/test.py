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

import os

data_path = os.path.join("..", "data", "test_data.csv")
if not os.path.isfile(data_path):
    data_path = os.path.join("..", data_path)
scdata = scprep.io.load_csv(data_path, cell_names=False)
scdata = scprep.filter.filter_empty_cells(scdata)
scdata = scprep.filter.filter_empty_genes(scdata)
scdata = scprep.filter.filter_duplicates(scdata)
scdata_norm = scprep.normalize.library_size_normalize(scdata)
scdata_norm = scprep.transform.sqrt(scdata_norm)


def test_genes_str_int():
    magic_op = magic.MAGIC(t="auto", decay=20, knn=10, verbose=False)
    str_gene_magic = magic_op.fit_transform(scdata_norm, genes=["VIM", "ZEB1"])
    int_gene_magic = magic_op.fit_transform(
        scdata_norm, graph=magic_op.graph, genes=[-2, -1]
    )
    assert str_gene_magic.shape[0] == scdata_norm.shape[0]
    np.testing.assert_array_equal(str_gene_magic, int_gene_magic)


def test_pca_only():
    magic_op = magic.MAGIC(t="auto", decay=20, knn=10, verbose=False)
    pca_magic = magic_op.fit_transform(scdata_norm, genes="pca_only")
    assert pca_magic.shape[0] == scdata_norm.shape[0]
    assert pca_magic.shape[1] == magic_op.n_pca


def test_all_genes():
    magic_op = magic.MAGIC(t="auto", decay=20, knn=10, verbose=False, random_state=42)
    int_gene_magic = magic_op.fit_transform(scdata_norm, genes=[-2, -1])
    magic_all_genes = magic_op.fit_transform(scdata_norm, genes="all_genes")
    assert scdata_norm.shape == magic_all_genes.shape
    int_gene_magic2 = magic_op.transform(scdata_norm, genes=[-2, -1])
    np.testing.assert_allclose(int_gene_magic, int_gene_magic2, rtol=0.015)


def test_all_genes_approx():
    magic_op = magic.MAGIC(
        t="auto", decay=20, knn=10, verbose=False, solver="approximate", random_state=42
    )
    int_gene_magic = magic_op.fit_transform(scdata_norm, genes=[-2, -1])
    magic_all_genes = magic_op.fit_transform(scdata_norm, genes="all_genes")
    assert scdata_norm.shape == magic_all_genes.shape
    int_gene_magic2 = magic_op.transform(scdata_norm, genes=[-2, -1])
    np.testing.assert_allclose(int_gene_magic, int_gene_magic2, atol=0.003, rtol=0.008)


def test_dremi():
    magic_op = magic.MAGIC(t="auto", decay=20, knn=10, verbose=False)
    # test DREMI: need numerical precision here
    magic_op.set_params(random_state=42)
    magic_op.fit(scdata_norm)
    dremi = magic_op.knnDREMI("VIM", "ZEB1", plot=True)
    np.testing.assert_allclose(dremi, 1.466004, atol=0.0000005)


def test_solver():
    # Testing exact vs approximate solver
    magic_op = magic.MAGIC(
        t="auto", decay=20, knn=10, solver="exact", verbose=False, random_state=42
    )
    data_imputed_exact = magic_op.fit_transform(scdata_norm)
    # should have exactly as many genes stored
    assert magic_op.X_magic.shape[1] == scdata_norm.shape[1]
    # should be nonzero
    assert np.all(data_imputed_exact >= 0)

    magic_op = magic.MAGIC(
        t="auto",
        decay=20,
        knn=10,
        n_pca=150,
        solver="approximate",
        verbose=False,
        random_state=42,
    )
    # magic_op.set_params(solver='approximate')
    data_imputed_apprx = magic_op.fit_transform(scdata_norm)
    # should have n_pca genes stored
    assert magic_op.X_magic.shape[1] == 150
    # make sure they're close-ish
    np.testing.assert_allclose(data_imputed_apprx, data_imputed_exact, atol=0.15)
    # make sure they're not identical
    assert np.any(data_imputed_apprx != data_imputed_exact)


def test_anndata():
    try:
        anndata
    except NameError:
        # anndata not installed
        return
    scdata = anndata.read_csv(data_path)
    fast_magic_operator = magic.MAGIC(
        t="auto", solver="approximate", decay=None, knn=10, verbose=False
    )
    sc_magic = fast_magic_operator.fit_transform(scdata, genes="all_genes")
    assert np.all(sc_magic.var_names == scdata.var_names)
    assert np.all(sc_magic.obs_names == scdata.obs_names)
    sc_magic = fast_magic_operator.fit_transform(scdata, genes=["VIM", "ZEB1"])
    assert np.all(sc_magic.var_names.values == np.array(["VIM", "ZEB1"]))
    assert np.all(sc_magic.obs_names == scdata.obs_names)
