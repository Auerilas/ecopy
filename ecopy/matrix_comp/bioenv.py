import numpy as np
from pandas import DataFrame
from scipy.spatial.distance import pdist, squareform
from scipy.stats import spearmanr
from itertools import combinations

def bioenv(dist, vars_df, columns=None):
    '''
    Docstring for function ecopy.simper
    ====================
    Find best subset of environmental variables maximally correlated with 
        community distances
    
    Use
    ----
    biodiv(dist, vars_df, columns=False)
    
    Returns a pandas.DataFrame containing the "best" subset of variables
    at each subset size, as well as the correlation coefficient of each

    Parameters
    ----------
    dist: A numpy.ndarray containing distances between objects
    vars_df: A pandas.DataFrame containing columns of variables 
        associated with the objects in `distance_matrix`
    columns: An iterable of strs, optional column names in `vars_df` to include 
        as variables in the calculations. If not provided, defaults to all 
        columns in `vars_df`

    Examples
    --------
    import numpy as np
    import pandas as pd
    import ecopy as ep

    dm = np.array([[0.0, 0.5, 0.25, 0.75],
                  [0.5, 0.0, 0.1, 0.42],
                  [0.25, 0.1, 0.0, 0.33],
                  [0.75, 0.42, 0.33, 0.0]])

    df = pd.DataFrame([[7.0, 400],
                       [8.0, 530],
                       [7.5, 450],
                       [8.5, 810]],
                      index=['A','B','C','D'],
                      columns=['pH', 'Elevation'])

    be = ep.bioenv(dm, df)
    print(be)
    '''
    if not isinstance(dist, (np.ndarray, DataFrame)):
        msg = 'Must provide a numpy.ndarray or pandas.DataFrame as input'
        raise TypeError(msg)
    if dist.shape[0] != dist.shape[1]:
        msg = 'Matrix dist must be a square, symmetric distance matrix'
        raise ValueError(msg)
    if not np.allclose(dist.T, dist):
        msg = 'Matrix dist must be a square, symmetric distance matrix'
        raise ValueError(msg)
    if np.any(dist < 0):
        msg = 'Distance matrix cannot have negative values'
        raise ValueError(msg)

    if not isinstance(vars_df, DataFrame):
        msg = 'Must provide a pandas.DataFrame as input'
        raise TypeError(msg)

    if columns is None:
        columns = vars_df.columns.values.tolist()

    if len(set(columns)) != len(columns):
        msg = 'Duplicate column names are not supported'
        raise ValueError(msg)

    if len(columns) < 1:
        msg = 'Must provide at least one column'
        raise ValueError(msg)

    for column in columns:
        if column not in vars_df:
            raise ValueError("Column '%s' not in data frame." % column)

    try:
        vars_df = vars_df.astype(float)
    except ValueError:
        raise TypeError("All specified columns in the data frame must be "
                        "numeric.")

    # Scale the vars and extract the underlying numpy array from the data
    # frame. We mainly do this for performance as we'll be taking subsets of
    # columns within a tight loop and using a numpy array ends up being ~2x
    # faster.
    vars_array = _scale(vars_df).values
    dm_flat = squareform(dist, force='tovector', checks=False)

    num_vars = len(columns)
    var_idxs = np.arange(num_vars)

    # For each subset size, store the best combination of variables:
    #     (string identifying best vars, subset size, rho)
    max_rhos = np.empty(num_vars, dtype=[('vars', object),
                                         ('size', int),
                                         ('correlation', float)])
    for subset_size in range(1, num_vars + 1):
        max_rho = None
        for subset_idxs in combinations(var_idxs, subset_size):
            # Compute Euclidean distances using the current subset of
            # variables. pdist returns the distances in condensed form.
            vars_dm_flat = pdist(vars_array[:, subset_idxs],
                                 metric='euclidean')
            rho = spearmanr(dm_flat, vars_dm_flat)[0]

            # If there are ties for the best rho at a given subset size, choose
            # the first one in order to match vegan::bioenv's behavior.
            if max_rho is None or rho > max_rho[0]:
                max_rho = (rho, subset_idxs)

        vars_label = ', '.join([columns[i] for i in max_rho[1]])
        max_rhos[subset_size - 1] = (vars_label, subset_size, max_rho[0])

    return DataFrame.from_records(max_rhos, index='vars')

def _scale(df):
    df = df.copy()
    df -= df.mean()
    df /= df.std()

    if df.isnull().any().any():
        raise ValueError('Column(s) in the data frame could not be scaled, '
                         'likely because the column(s) had no variance.')
    return df