import numbers
import numpy as np
import pandas as pd
try:
    import anndata
except (ImportError, SyntaxError):
    # anndata not installed
    pass


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
        if params[p] is not x and params[p] != x:
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
    return X is Y or (isinstance(X, Y.__class__) and X.shape == Y.shape and
                      np.sum((X != Y).sum()) == 0)


def convert_to_same_format(data, target_data, columns=None):
    # create new data object
    if is_pandas(target_data, sparse=True):
        data = pd.SparseDataFrame(data)
        pandas = True
    elif is_pandas(target_data):
        data = pd.DataFrame(data)
        pandas = True
    elif is_anndata(target_data):
        data = anndata.AnnData(data)
        pandas = False
    else:
        # nothing to do
        return data
    # retrieve column names
    target_columns = target_data.columns if pandas else target_data.var
    # subset column names
    try:
        if columns is not None:
            if pandas:
                target_columns = target_columns[columns]
            else:
                target_columns = target_columns.iloc[columns]
    except (KeyError, IndexError):
        # keep the original column names
        target_columns = columns
    # set column names on new data object
    if pandas:
        data.columns = target_columns
        data.index = target_data.index
    else:
        data.var = target_columns
        data.obs = target_data.obs
    return data


def in_ipynb():
    """Check if we are running in a Jupyter Notebook

    Credit to https://stackoverflow.com/a/24937408/3996580
    """
    __VALID_NOTEBOOKS = ["<class 'google.colab._shell.Shell'>",
                         "<class 'ipykernel.zmqshell.ZMQInteractiveShell'>"]
    try:
        return str(type(get_ipython())) in __VALID_NOTEBOOKS
    except NameError:
        return False


def is_pandas(data, sparse=None):
    try:
        if sparse is True:
            return isinstance(data, pd.SparseDataFrame)
        elif sparse is False:
            return isinstance(data, pd.DataFrame) and not \
                isinstance(data, pd.SparseDataFrame)
        elif sparse is None:
            return isinstance(data, pd.DataFrame)
        else:
            raise ValueError("Expected sparse in [True, False, None]. "
                             "Got {}".format(sparse))
    except NameError:
        # pandas not installed
        return False


def is_anndata(data):
    try:
        return isinstance(data, anndata.AnnData)
    except NameError:
        # anndata not installed
        return False


def has_empty_columns(data):
    try:
        return np.any(np.array(data.sum(0)) == 0)
    except AttributeError:
        if is_anndata(data):
            return np.any(np.array(data.X.sum(0)))
        else:
            raise
