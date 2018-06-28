#!/usr/bin/env python


from __future__ import print_function, division
import doctest
import magic
import pandas as pd
import numpy as np


scdata = pd.read_csv("../data/HMLE_TGFb_day_8_10.csv")
scdata_norm = magic.preprocessing.library_size_normalize(scdata)
assert scdata.shape == scdata_norm

fast_magic_operator = magic.MAGIC(t='auto', a=20, k=10)
classic_magic_operator = magic.MAGIC(t=10, a=20, k=10, n_pca=None)

fast_magic = fast_magic_operator.fit_transform(scdata_norm)
classic_magic = classic_magic_operator.fit_transform(scdata_norm)



def test_magic():
    doctest.testmod()


def test_scdata():

    return 0


if __name__ == "main":
    test_magic()
