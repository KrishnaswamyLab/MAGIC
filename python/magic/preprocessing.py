# author: Daniel Burkhardt <daniel.burkhardt@yale.edu>
# (C) 2017 Krishnaswamy Lab GPLv2

from __future__ import print_function, division
import warnings
import scprep


def library_size_normalize(data, verbose=False):
    """Performs L1 normalization on input data
    Performs L1 normalization on input data such that the sum of expression
    values for each cell sums to 1
    then returns normalized matrix to the metric space using median UMI count
    per cell effectively scaling all cells as if they were sampled evenly.

    Parameters
    ----------
    data : ndarray [n,p]
        2 dimensional input data array with n cells and p dimensions

    Returns
    -------
    data_norm : ndarray [n, p]
        2 dimensional array with normalized gene expression values
    """
    warnings.warn("magic.preprocessing is deprecated. "
                  "Please use scprep.normalize instead. "
                  "Read more at http://scprep.readthedocs.io",
                  FutureWarning)
    return scprep.normalize.library_size_normalize(data)
