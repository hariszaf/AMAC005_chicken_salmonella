from variables import *

import numpy as np
import pandas as pd

import scipy.stats as stats
import statsmodels.api as sm

from scipy.stats import rankdata, spearmanr


def wls(df):
    """Weighted Least Squares (WLS) Regression"""

    # Fit Weighted Least Squares (WLS) model
    weights = (
        1 / df[NEIGHBORS_NUMBER].value_counts().reindex(df[NEIGHBORS_NUMBER]).values
    )
    X = sm.add_constant(df[NEIGHBORS_NUMBER])
    y = df[RATIO]
    model = sm.WLS(y, X, weights=weights).fit()

    # Create regression line for plot
    x_range           = np.linspace(df[NEIGHBORS_NUMBER].min(), df[NEIGHBORS_NUMBER].max(), 100)
    X_pred            = sm.add_constant(x_range)
    y_pred            = model.predict(X_pred)
    wls_regression_df = pd.DataFrame({NEIGHBORS_NUMBER: x_range, RATIO: y_pred})

    # Compute R² and p-value
    r_squared = model.rsquared
    p_value   = model.pvalues[NEIGHBORS_NUMBER]

    return wls_regression_df, r_squared, p_value


def ols(df):
    """Ordinary Least Squares (OLS) Regression"""

    # Fit Ordinary Least Squares (OLS) model
    X = sm.add_constant(df[NEIGHBORS_NUMBER])
    y = df[RATIO]
    model = sm.OLS(y, X).fit()

    # Create regression line for plot
    x_range           = np.linspace(df[NEIGHBORS_NUMBER].min(), df[NEIGHBORS_NUMBER].max(), 100)
    X_pred            = sm.add_constant(x_range)
    y_pred            = model.predict(X_pred)
    ols_df = pd.DataFrame({NEIGHBORS_NUMBER: x_range, RATIO: y_pred})

    # Compute R² and p-value
    r_squared = model.rsquared
    p_value   = model.pvalues[NEIGHBORS_NUMBER]

    return ols_df, r_squared, p_value


def mean_genome_abd(abd_df, meta_df, bin_id, type, case):
    """Function to get mean and std for a genome's abundance on a certain type-case combination"""

    samples = meta_df.columns[meta_df.loc[type] == case]
    sub_df = abd_df[["genome"] + samples.tolist() + ["classification"]]

    findings  = sub_df[sub_df["genome"] == bin_id].select_dtypes(include=['number'])
    mean      = findings.mean(axis=1)
    std       = findings.std(axis=1)
    formatted = mean.round(4).astype(str) + " ± " + std.round(4).astype(str)

    return formatted


def fit_polynomial_regression(df, degree=2):
    # Create polynomial features
    X = np.vander(df[NEIGHBORS_NUMBER], N=degree + 1, increasing=True)

    # Fit the polynomial regression model
    y = df[RATIO]
    model = sm.OLS(y, X).fit()

    # Create regression curve
    x_range = np.linspace(df[NEIGHBORS_NUMBER].min(), df[NEIGHBORS_NUMBER].max(), 100)
    X_pred  = np.vander(x_range, N=degree + 1, increasing=True)
    y_pred  = model.predict(X_pred)

    poly_regression_df = pd.DataFrame({NEIGHBORS_NUMBER: x_range, RATIO: y_pred})

    # Compute R² and p-values
    r_squared = model.rsquared
    p_values  = model.pvalues  # p-values for each coefficient

    return poly_regression_df, r_squared, p_values


def weighted_pearson(df, col_x, col_y, weight_col):
    """Compute weighted Pearson correlation for pandas DataFrame columns."""

    x = df[col_x]
    y = df[col_y]
    w = df[weight_col]

    x_mean = np.average(x, weights=w)
    y_mean = np.average(y, weights=w)
    cov_xy = np.sum(w * (x - x_mean) * (y - y_mean))
    std_x  = np.sqrt(np.sum(w * (x - x_mean) ** 2))
    std_y  = np.sqrt(np.sum(w * (y - y_mean) ** 2))

    return cov_xy / (std_x * std_y)


def weighted_spearman(df, col_x, col_y, weights_col):
    """Compute weighted Spearman correlation for two columns in a DataFrame."""

    x = df[col_x].values  # Get values from the first column
    y = df[col_y].values  # Get values from the second column
    weights = df[weights_col].values  # Get weights from the specified weights column

    x_rank = rankdata(x)  # Rank the x values
    y_rank = rankdata(y)  # Rank the y values

    return weighted_pearson(x_rank, y_rank, weights)


def weighted_pearson_with_pvalue(df, col_x, col_y, weight_col):
    """Compute weighted Pearson correlation and its p-value for pandas DataFrame columns."""

    x = df[col_x]
    y = df[col_y]
    w = df[weight_col]

    # Calculate weighted mean
    x_mean = np.average(x, weights=w)
    y_mean = np.average(y, weights=w)

    # Calculate weighted covariance
    cov_xy = np.sum(w * (x - x_mean) * (y - y_mean))

    # Calculate weighted standard deviations
    std_x = np.sqrt(np.sum(w * (x - x_mean) ** 2))
    std_y = np.sqrt(np.sum(w * (y - y_mean) ** 2))

    # Weighted Pearson correlation
    r = cov_xy / (std_x * std_y)

    # Calculate the number of non-zero weights
    n = np.sum(w > 0)

    # Calculate t-statistic
    t_stat = r * np.sqrt((n - 2) / (1 - r**2))

    # Calculate p-value from t-distribution
    p_value = 2 * (1 - stats.t.cdf(np.abs(t_stat), df=n - 2))

    return r, p_value


def spearman_corr(df, col_x, col_y):
    """Compute Spearman correlation for two columns in a DataFrame."""

    x = df[col_x].values  # Get values from the first column
    y = df[col_y].values  # Get values from the second column

    x_rank = rankdata(x)  # Rank the x values
    y_rank = rankdata(y)  # Rank the y values

    # Compute the Spearman correlation using the ranked values
    return spearmanr(x_rank, y_rank).correlation
