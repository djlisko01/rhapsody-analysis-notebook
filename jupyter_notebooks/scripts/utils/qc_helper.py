from anndata import AnnData
import numpy as np
from scipy.stats import median_abs_deviation
from enum import Enum
from typing import Union


class MethodOptions(Enum):
    MAD = "MAD"
    MANUAL = "MANUAL"


class OutlierFinder:
    def __init__(self):
        self._mad_thresholds = {
            "log_p_count_threshold": 5,
            "top_genes_threshold": 5,
            "log_n_genes_threshold": 5,
            "pct_counts_mt_threshold": 20,
            "pct_mt_mad_threshold": 5
        }

        self._target_thresholds = {
            "n_genes_by_count": None,
            "count_depth": None,
            "pct_count_mt": None
        }

    def find_outliers(self, adata: AnnData, method: MethodOptions):
        if method not in MethodOptions:
            raise ValueError(f"Method should be {MethodOptions.MAD} or {MethodOptions.MANUAL}")

        if method == MethodOptions.MAD:
            return self.find_outliers_mad(adata)

        return self.set_outliers_manually(adata)

    @staticmethod
    def is_mod_outlier(adata, metric: str, n_mads: int):
        M = adata.obs[metric]
        outliers = (M < np.median(M) - n_mads * median_abs_deviation(M)) | (
                np.median(M) + n_mads * median_abs_deviation(M) < M
        )
        return outliers

    def set_outliers_manually(self, adata: AnnData):
        obs_df = adata.obs
        gene_outlier = obs_df.n_genes_by_counts < self._target_thresholds["n_genes_by_count"]
        total_count_outliers = obs_df.total_counts < self._target_thresholds["count_depth"]
        mt_outliers = obs_df.pct_counts_mt > self._target_thresholds["pct_count_mt"]
        return gene_outlier | total_count_outliers | mt_outliers

    def get_target_thresholds(self):
        return self._target_thresholds

    def find_outliers_mad(self, adata: AnnData):
        is_outlier = self.is_mod_outlier
        mt_outlier = adata.obs.pct_counts_mt > self._mad_thresholds["pct_counts_mt_threshold"]
        pct_mt_mad_outlier = is_outlier(adata, "pct_counts_mt", self._mad_thresholds["pct_mt_mad_threshold"])
        log_p_count_outlier = is_outlier(adata, "log1p_total_counts", self._mad_thresholds["log_p_count_threshold"])
        top_genes_outlier = is_outlier(adata, "pct_counts_in_top_20_genes", self._mad_thresholds["top_genes_threshold"])
        log_n_genes = is_outlier(adata, "log1p_n_genes_by_counts", self._mad_thresholds["log_n_genes_threshold"])

        return mt_outlier | log_p_count_outlier | top_genes_outlier | log_n_genes | pct_mt_mad_outlier

    def set_mad_threshold(self, metric_set: dict) -> None:
        s1 = set(metric_set.keys())
        s2 = set(metric_set.keys())

        if not s2.issubset(s1):
            raise ValueError(f"Only set these metric: {s1}")
        self._mad_thresholds.update(metric_set)

    def set_manual_threshold(self, metric_set: dict) -> None:
        s1 = set(metric_set.keys())
        s2 = set(metric_set.keys())

        if not s2.issubset(s1):
            raise ValueError(f"Only set these metric: {s1}")

        self._target_thresholds.update(metric_set)

    def get_mad_option(self) -> dict:
        return self._mad_thresholds

    @staticmethod
    def calculate_threshold_values(data, n_mads):
        median_val = np.median(data)
        mad_val = median_abs_deviation(data)
        lower_threshold = median_val - n_mads * mad_val
        upper_threshold = median_val + n_mads * mad_val
        return lower_threshold, upper_threshold


