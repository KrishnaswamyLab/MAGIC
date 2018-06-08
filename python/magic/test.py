#!/usr/bin/env python

# Generating random fractal tree via DLA
from __future__ import print_function, division
import doctest
import pandas as pd
import phate
import magic
import numpy as np


clusters = pd.read_csv("../data/MAP.csv", header=None)
clusters.columns = pd.Index(['wells', 'clusters'])
bmmsc = pd.read_csv("../data/BMMC_myeloid.csv.gz", index_col=0)

C = clusters['clusters']  # using cluster labels from original publication

# library_size_normalize performs L1 normalization on each cell
bmmsc_norm = phate.preprocessing.library_size_normalize(bmmsc)
bmmsc_norm = np.sqrt(bmmsc_norm)
magic_operator = magic.MAGIC(
    t='auto', a=20, k=10)
phate_operator = phate.PHATE(
    n_components=2, t='auto', a=200,
    k=10, mds='metric', mds_dist='euclidean',
    n_landmark=1000)

y_magic = magic_operator.fit_transform(bmmsc_norm)
y_phate = phate_operator.fit_transform(bmmsc_norm)


def test_magic():
    doctest.testmod()


def test_bmmsc():

    return 0


if __name__ == "main":
    test_magic()
