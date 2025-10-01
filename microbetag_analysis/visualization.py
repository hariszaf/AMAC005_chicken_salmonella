import seaborn as sns
from venn import venn
from plotnine import (
    ggplot,
    aes,
    geom_point,
    theme_minimal,
    labs,
    ggsave,
    geom_line,
    annotate,
)
import scipy.stats as stats

from variables import *
from stats import *


def plot_neighbors_per_seed_compl(df, condition, model="WLS"):
    """
    Plots the total number of seed complements a species is invovlved in, against the number of 
    its neighboring taxa.
    """

    if model == "WLS":
        regression_df, r_squared, p_value = wls(df)

    elif model == "OLS":
        regression_df, r_squared, p_value = ols(df)

    else:
        print("Model not recognized. Please use 'WLS' or 'OLS'.")
        # regression_df, r_squared, p_value =

    # Generate the plot
    plot = (
        ggplot(df, aes(
            x     = NEIGHBORS_NUMBER,
            y     = COMPLEMENTS_NUMBER,
            color = "order")
        )
        + geom_point(size=3)
        + geom_line(regression_df, aes(x=NEIGHBORS_NUMBER, y=RATIO), color="red")
        + theme_minimal()
        + labs(
            title = f"{condition} at the family level",
            x     = "Neighbors Number",
            y     = "Seed Complements Number",
        )
        + annotate(
            "text",
            x     = df[NEIGHBORS_NUMBER].max() * 0.8,
            y     = df[RATIO].max() * 0.9,
            label = f"R² = {r_squared:.3f}\nP-value = {p_value:.3g}",
            size  = 10,
            color = "black",
        )
    )

    return plot, df


def plot_regression(ax, df, x_col, y_col, title):
    """Plots a regression plot on the given Axes object and includes p-value and R-squared."""

    # Perform linear regression to get p-value and R-squared
    slope, intercept, r_value, p_value, std_err = stats.linregress(df[x_col], df[y_col])
    r_squared = r_value**2

    # Create the regression plot
    sns.regplot(x=df[x_col], y=df[y_col], scatter_kws={"alpha": 0.7}, line_kws={"color": "red"}, ax=ax)

    # Add p-value and R-squared to the title
    ax.set_title(f"{title}\n$R^2$ = {r_squared:.2f}, p-value = {p_value:.4f}")


def plot_dual_regression(ax, df, x_col, y1_col, y2_col, title, color1="blue", color2="red"):
    """Plots two regressions with separate y-axes: one for y1_col (left) and y2_col (right)."""

    # Left Y-axis (cooperation vs cooccurrence)
    ax1 = ax
    (
        slope1,
        intercept1,
        r_value1,
        p_value1,
        std_err1
    ) = stats.linregress(df[x_col], df[y1_col])

    r_squared1 = r_value1**2

    sns.regplot(
        x           = df[x_col],
        y           = df[y1_col],
        scatter_kws = {"alpha": 0.7},
        line_kws    = {"color": color1},
        ax          = ax1
    )
    ax1.set_ylabel(y1_col, color=color1)
    ax1.tick_params(axis="y", colors=color1)

    # Right Y-axis (competition vs cooccurrence)
    ax2 = ax1.twinx()
    (
        slope2,
        intercept2,
        r_value2,
        p_value2,
        std_err2
    ) = stats.linregress(df[x_col], df[y2_col])

    r_squared2 = r_value2**2

    sns.regplot(
        x           = df[x_col],
        y           = df[y2_col],
        scatter_kws = {"alpha": 0.3, "color": "red"},
        line_kws    = {"color": color2},
        ax          = ax2
    )
    ax2.set_ylabel(y2_col, color=color2)
    ax2.tick_params(axis="y", colors=color2)

    # Titles
    ax1.set_xlabel(x_col)
    ax1.set_title(
        f"{title}\n{y1_col}: R²={r_squared1:.2f}, p={p_value1:.4f} | {y2_col}: R²={r_squared2:.2f}, p={p_value2:.4f}"
    )


# def plot_regression(ax, df, x1_col, x2_col, y_col, title, color1="red", color2="blue"):
#     """Plots two regression plots on the same Axes object and includes p-value & R-squared."""

#     # First regression (x1_col vs y_col)
#     slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(df[x1_col], df[y_col])
#     r_squared1 = r_value1**2
#     sns.regplot(x=df[x1_col], y=df[y_col], scatter_kws={"alpha": 0.7}, line_kws={"color": color1}, ax=ax, label=f"{x1_col} (R²={r_squared1:.2f}, p={p_value1:.4f})")

#     # Second regression (x2_col vs y_col)
#     slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(df[x2_col], df[y_col])
#     r_squared2 = r_value2**2
#     sns.regplot(x=df[x2_col], y=df[y_col], scatter_kws={"alpha": 0.7}, line_kws={"color": color2}, ax=ax, label=f"{x2_col} (R²={r_squared2:.2f}, p={p_value2:.4f})")

#     # Titles and legend
#     ax.set_title(title)
#     ax.legend()
