import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Circle
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from cycler import cycler
import seaborn as sns
from anndata import AnnData
from typing import Tuple
from scipy.stats import median_abs_deviation


def plot_rna_qc(adata: AnnData, zoom=None, zoom_bins=100, outliers=None):
    fig = plt.figure(layout='constrained', figsize=(15, 4))
    l_fig, r_fig = fig.subfigures(1, 2, wspace=0.07, width_ratios=[2, 1])
    axs_left = l_fig.subplots(ncols=2)
    obs_df = adata.obs
    hist_ax = plot_hist(x_var=obs_df["total_counts"], ax=axs_left[0],
                        zoom=zoom, zoom_bins=zoom_bins, outliers=outliers)
    hist_ax.set_xlabel("Count Depth")
    _ = plot_genes_vs_total(obs_df=obs_df, ax=axs_left[1])
    plot_hist_viol(y_var=obs_df.pct_counts_mt, fig=r_fig)

    if outliers is not None:
        _add_threshold_line(adata, outliers, fig)
    return fig


def _add_threshold_line(adata: AnnData, outliers, fig: plt.Figure) -> None:
    axs = fig.get_axes()
    
    # n_gene outliers
    axs[1].hlines(min(adata.obs["n_genes_by_counts"][~outliers]), *axs[1].get_xlim(), color="red")
    axs[1].hlines(max(adata.obs["n_genes_by_counts"][~outliers]), *axs[1].get_xlim(), color="red")
    
    # total cout outliers
    axs[1].vlines(min(adata.obs["total_counts"][~outliers]), *axs[1].get_ylim(), color="red")
    axs[1].vlines(max(adata.obs["total_counts"][~outliers]), *axs[1].get_ylim(), color="red")
    
    
    axs[3].hlines(max(adata.obs["pct_counts_mt"][~outliers]), -0.5, 0.5, color="red", zorder=10)
    axs[4].hlines(max(adata.obs["pct_counts_mt"][~outliers]), *axs[4].get_xlim(), color="red")


def plot_10x_fig(adata: AnnData,):
    count_data = adata.obs["total_counts"].sort_values(ascending=False)
    order = range(1, len(count_data) + 1)
    plt.semilogy(order, count_data)
    plt.xlabel("Barcode Rank")
    plt.ylabel("Total Count Depth")


def plot_hist(x_var, n_bins=100, zoom: Tuple[int, int] = None, zoom_bins=100, ax=None, outliers=None):
    if ax is None:
        ax = plt.axes()
    hist_ax = sns.histplot(data=x_var, bins=n_bins, ax=ax)
    if zoom is not None:
        zoom_ax = hist_ax.inset_axes([0.5, 0.5, 0.47, 0.47])
        zoom_ax.hist(x_var, bins=zoom_bins)

        zoom_ax.set_xlim(zoom)
        zoom_ax.spines["right"].set_visible(False)
        zoom_ax.spines["top"].set_visible(False)
        n_ticks = np.int64(zoom_ax.get_xlim())
        steps = max(n_ticks) // 4
        zoom_ax.set_xticks([i for i in range(*n_ticks, steps)])

        if outliers is not None:
            zoom_ax.vlines(min(x_var[~outliers]), *zoom_ax.get_ylim(), color="red")

    return hist_ax


def plot_qc_violin(y_var, ax=None, outliers=None):

    len_y_var = len(y_var)
    alpha = np.array([0.6] * len_y_var)
    colors = np.repeat("black", len_y_var)

    if outliers is not None:
        alpha[outliers] = 0.1
        colors[outliers] = "grey"

    ax = sns.stripplot(y_var, color=colors, alpha=alpha, jitter=0.4, size=1.4, ax=ax)
    ax = sns.violinplot(y_var, ax=ax, inner=None)
    ax.set_xticks([])
    return ax


def plot_genes_vs_total(obs_df: pd.DataFrame, ax=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 3))

    sctr = ax.scatter(
        x=obs_df.total_counts,
        y=obs_df.n_genes_by_counts,
        c=obs_df.pct_counts_mt,
        s=1,
    )
    ax.set_ylabel("Number of Genes")
    ax.set_xlabel("Count Depth")

    plt.colorbar(sctr, ax=ax, shrink=0.6)
    return ax


def plot_hist_viol(y_var, fig=None, bins=100):
    if fig is None:
        fig = plt.figure(figsize=(6, 4))

    gs = fig.add_gridspec(nrows=1, ncols=2, width_ratios=(4, 1),
                          left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0.01, hspace=0.05)
    # For scatter Plot
    ax = fig.add_subplot(gs[0, 0])
    ax.set_title('Distribution of  Mitochondrial RNA')
    ax.set_ylabel("% Mitochondrial Genes")

    # Histogram
    ax_histy = fig.add_subplot(gs[0, 1], sharey=ax)
    ax_histy.tick_params(axis="both", labelleft=False)
    ax_histy.set_xticks([])

    # Remove borders
    for spine in ax_histy.spines:
        if spine != "left":
            ax_histy.spines[spine].set_visible(False)

    ax_histy.hist(y_var, bins=bins, orientation='horizontal', edgecolor='black', linewidth=0.25)
    plot_qc_violin(y_var, ax=ax)


def calculate_descriptive(adata: AnnData, mt_threshold: float, sample_name, is_rhapsody=False) -> pd.DataFrame:
    pd.options.display.float_format = '{:,.2f}'.format
    obs_df = adata.obs
    pct_above_threshold = np.sum(obs_df.pct_counts_mt < mt_threshold) / len(obs_df) * 100
    stats = {
            "Cell Count": len(obs_df.index),
            "Number Reads": sum(obs_df.total_counts),
            "Median Genes Per Cell": obs_df.n_genes_by_counts.median(),
            "% MT Per Cell": obs_df.pct_counts_mt.median(),
            f"Cells Pass {mt_threshold}% Threshold": f"{pct_above_threshold:.2f}%",
        }

    if is_rhapsody:
        results = obs_df["Sample_Tag"].value_counts()
        rhap_stats = {
            "Number of Multi Tagged Cells": results["Multiplet"],
            "% Multi Tagged Cells":  f"{(results['Multiplet'] / np.sum(results)) * 100:.2f}%",
            "Number of Undetermined Tags": results["Undetermined"],
            "% Undetermined Cells Tagged": f"{(results['Undetermined'] / np.sum(results)) * 100:.2f}%",
        }
        stats.update(rhap_stats)

    df = pd.DataFrame(stats, index=[sample_name])
    return df.round(2).T


def get_per_sample_stats(adata: AnnData, mt_threshold: float, samples_col: str) -> pd.DataFrame:
    sample_stats = adata.obs.groupby(samples_col)[[
        "pct_counts_mt",
        "total_counts",
        "n_genes_by_counts",
    ]].agg(["mean"])

    if sample_stats.index.isin(["Multiplet", "Undetermined"]).any():
        sample_stats.drop(["Multiplet", "Undetermined"], inplace=True)

    sample_stats.columns = ["Mean % MT Per Cell", "Mean # Reads Per Cell", "Mean # Gene Per Cell"]
    
    sample_stats[f"Cells Pass {mt_threshold}% Threshold"] = adata.obs.groupby("Sample_Tag").apply(
        lambda x: f"{(np.sum(x['pct_counts_mt'] < mt_threshold) / len(x) * 100):.2f}%"
    )

    sample_stats["Cell Count"] = adata.obs.groupby(samples_col).size()
    sample_stats["Read Depth"] = adata.obs.groupby(samples_col)["total_counts"].sum()
    return sample_stats

def df_to_fig(df: pd.DataFrame, ax:plt.Axes=None, colors: list=None):
    if ax is None:
        ax = plt.axes()
    nrows, _ = df.shape

    # Compute column widths based on the longest string in each column
    col_widths = [len(col) for col in df.columns]
    row_height = 1.5  # Increase this value for more space between rows
    for idx, col in enumerate(df.columns):
        col_widths[idx] = max(col_widths[idx], df[col].astype(str).apply(len).max())

    if colors is not None:
        col_widths.insert(0, 2)  # Add width for the color column with extra space

    # Adjust the x-axis limits
    total_width = sum(col_widths) + len(col_widths) - 1  # Including gaps for clarity
    ax.set_xlim(0, total_width)
    row_height = 1.5  # As previously set
    header_spacing = 0.75  # This adds space above and below the headers
    ax.set_ylim((-1.5 - header_spacing) * row_height, nrows * row_height)  # Adjusted y-limits

    df.reset_index(inplace=True, drop=True)
    
    for idx, row in df.iterrows():
        y_position = (idx + header_spacing) * row_height  # Adjust the y-coordinate considering header spacing
        x = col_widths[0] / 2  # start at the center of the first column
        
        
        if colors is not None:
            circle = Circle((x, y_position), 0.4, color=colors[idx])
            ax.add_patch(circle)
            x += col_widths[0]  # Move to the right edge of the color column
        
        for j, col in enumerate(df.columns):
            ax.text(x + col_widths[j + 1] / 2, y_position, row[col], ha="center", va='center', fontsize=20)  # Increase fontsize
            x += col_widths[j + 1]  # Adjust for the next column

    __add_headers(df, ax, col_widths)
    ax.invert_yaxis()
    ax.axis('off')  # Hide the axes
    ax.set_aspect("equal")
    return ax

def __add_headers(df, ax, col_widths):
    columns = df.columns
    x = col_widths[0]  # start at the right edge of the color column
    
    for j, header in enumerate(columns):
        ax.text(x + col_widths[j + 1] / 2, -1, header, weight='bold', ha='center', va='center', fontsize=20)
        x += col_widths[j + 1]  # Move to the center of the next column
    ax.axhline(y=-0.5, color='black', linewidth=2)



def __get_column_widths(df: pd.DataFrame) -> list:
    # Compute column widths based on the longest string in each column
    col_widths = []
    for col in df.columns:
        max_data_width = df[col].astype(str).apply(len).max()  # get maximum string length in the column
        header_width = len(col)
        max_width = max(max_data_width, header_width)  # take the maximum between header width and data width
        col_widths.append(max_width)

    return col_widths


if __name__ == "__main__":
    from utils.rhapsody_reader import RhapsodyReader
    import scanpy as sc

    rhaps_reader = RhapsodyReader()
    seq_file = "../../data/mu_data/markersCombined_seq1-p121C_RSEC_MolsPerCell.h5mu"
    mdata = rhaps_reader.read_mudata(seq_file)
    rna = mdata.mod["rna"]

    # mitochondrial genes
    rna.var["mt"] = rna.var_names.str.startswith("MT-")

    # ribosomal genes
    rna.var["ribo"] = rna.var_names.str.startswith(("RPS", "RPL"))

    # hemoglobin genes
    rna.var["hb"] = rna.var_names.str.contains(("^HB[^(P)]"))

    sc.pp.calculate_qc_metrics(
        rna, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
    )
