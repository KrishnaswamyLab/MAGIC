import numbers
import numpy as np
import pandas as pd
from scipy import sparse


def check_positive(**params):
    """Check that parameters are positive as expected

    Raises
    ------
    ValueError : unacceptable choice of parameters
    """
    for p in params:
        if params[p] <= 0:
            raise ValueError(
                "Expected {} > 0, got {}".format(p, params[p]))


def check_int(**params):
    """Check that parameters are integers as expected

    Raises
    ------
    ValueError : unacceptable choice of parameters
    """
    for p in params:
        if not isinstance(params[p], numbers.Integral):
            raise ValueError(
                "Expected {} integer, got {}".format(p, params[p]))


def check_if_not(x, *checks, **params):
    """Run checks only if parameters are not equal to a specified value

    Parameters
    ----------

    x : excepted value
        Checks not run if parameters equal x

    checks : function
        Unnamed arguments, check functions to be run

    params : object
        Named arguments, parameters to be checked

    Raises
    ------
    ValueError : unacceptable choice of parameters
    """
    for p in params:
        if params[p] is not x:
            [check(p=params[p]) for check in checks]


def check_in(choices, **params):
    """Checks parameters are in a list of allowed parameters

    Parameters
    ----------

    choices : array-like, accepted values

    params : object
        Named arguments, parameters to be checked

    Raises
    ------
    ValueError : unacceptable choice of parameters
    """
    for p in params:
        if params[p] not in choices:
            raise ValueError(
                "{} value {} not recognized. Choose from {}".format(
                    p, params[p], choices))


def check_between(v_min, v_max, **params):
    """Checks parameters are in a specified range

    Parameters
    ----------

    v_min : float, minimum allowed value (inclusive)

    v_max : float, maximum allowed value (inclusive)

    params : object
        Named arguments, parameters to be checked

    Raises
    ------
    ValueError : unacceptable choice of parameters
    """
    for p in params:
        if params[p] < v_min or params[p] > v_max:
            raise ValueError("Expected {} between {} and {}, "
                             "got {}".format(p, v_min, v_max, params[p]))


def matrix_is_equivalent(X, Y):
    """
    Checks matrix equivalence with numpy, scipy and pandas
    """
    return isinstance(X, Y.__class__) and X.shape == Y.shape and \
        np.sum((X != Y).sum()) == 0


def select_cols(data, idx):
    if isinstance(data, pd.DataFrame):
        try:
            data = data.loc[:, idx]
        except KeyError:
            if isinstance(idx, numbers.Integral) or \
                    issubclass(np.array(idx).dtype.type, numbers.Integral):
                data = data.loc[:, np.array(data.columns)[idx]]
            else:
                raise
    else:
        if isinstance(data, (sparse.coo_matrix,
                             sparse.bsr_matrix,
                             sparse.lil_matrix,
                             sparse.dia_matrix)):
            data = data.tocsr()
        data = data[:, idx]
    return data


def convert_to_same_format(data, target_data, columns=None):
    try:
        if isinstance(target_data, pd.SparseDataFrame):
            data = pd.SparseDataFrame(data)
        elif isinstance(target_data, pd.DataFrame):
            data = pd.DataFrame(data)
        else:
            # nothing to do
            return data
        target_columns = target_data.columns
        if columns is not None:
            target_columns = target_columns[columns]
        data.columns = target_columns
        data.index = target_data.index
    except NameError:
        # pandas not installed
        pass
    return data
