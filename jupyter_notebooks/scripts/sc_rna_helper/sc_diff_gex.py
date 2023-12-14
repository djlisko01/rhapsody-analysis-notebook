from anndata import AnnData
import scanpy as sc
import numpy as np
from typing import Sequence
import random
import pandas as pd

class PseudoBulk:
    """
    A class for generating pseudo bulk data from single-cell RNA sequencing (scRNA-seq) data.

    Attributes:
        min_number_cells (int): Minimum number of cells required for a sample to be included in the pseudo bulk data.
        sample_key (str|int): The key in the AnnData object's metadata to identify samples.
        adata (AnnData): An AnnData object containing the scRNA-seq data.
        cell_id_key (str): The key in the AnnData object's metadata to identify cell types.
        reps_per_sample (int): Number of replicates to be generated per sample.
        samples_rm (set): A set of samples removed due to not meeting the minimum cell count.

    Args:
        adata (AnnData): An AnnData object containing scRNA-seq data.
        sample_key (str|int): The key in the AnnData object's metadata to identify samples.
        cell_id_key (str): The key in the AnnData object's metadata to identify cell types.
    """
    def __init__(self, adata: AnnData, sample_key: str | int, cell_id_key: str) -> None:
        self.min_number_cells = 30
        self.sample_key =sample_key
        self.adata = adata.copy()
        self.cell_id_key = cell_id_key
        self.reps_per_sample = 3
        self.samples_rm = set()

    def set_min_number_cells(self, sample_min_cells: int):
        """
        Sets the minimum number of cells required for a sample to be included in the pseudo bulk data.

        Args:
            sample_min_cells (int): The minimum number of cells required per sample.
        """
        self.min_number_cells = sample_min_cells

    def set_replicates_per_sample(self, reps_per_sample: int):
        """
        Sets the number of replicates to be generated for each sample in the pseudo bulk data.

        Args:
            reps_per_sample (int): The number of replicates per sample.
        """
        self.reps_per_sample = reps_per_sample

    def generate_pseudo_bulk(self, condition_key: str,  cell_type: str, raw_data_layer: str = None, cols_to_keep:Sequence[str]=None) -> AnnData:
        """
        Generates pseudo bulk data for a specific cell type across different samples.

        Args:
            condition_key (str): Key in the AnnData object's metadata used to define the condition or grouping.
            cell_type (str): The cell type for which to generate pseudo bulk data.
            raw_data_layer (str, optional): Specifies the layer of AnnData object to be used. If None, uses the default layer.

        Returns:
            AnnData: An AnnData object containing the aggregated pseudo bulk data.
        """
        pbs = []
        # Based on cell identity key
        cell_subset = self.adata[self.adata.obs[self.cell_id_key] == cell_type].copy()
        samples = self.adata.obs[self.sample_key].unique()

        if raw_data_layer:
            cell_subset.X = cell_subset.layers["raw"].A

        for sample in samples:
            sample_cell_subset = cell_subset[cell_subset.obs[self.sample_key] == sample]

            if sample_cell_subset.n_obs < self.min_number_cells:
                self.samples_rm.add(sample)
                continue

            X = np.array([sample_cell_subset.X.sum(axis=0)])
            var = sample_cell_subset.var[[]]

            if cols_to_keep:
                    obs = sample_cell_subset.obs[cols_to_keep].iloc[0]
                    obs_df = pd.DataFrame(obs).T
                    rep_adata = AnnData(X=X, var=var, obs=obs_df)
            else:
                rep_adata = AnnData(X=X, var=var)

            rep_adata.obs_names = [sample]
            rep_adata.obs[condition_key] = sample_cell_subset.obs[condition_key].iloc[0]
            rep_adata.obs["num_cells"] = sample_cell_subset.n_obs
            rep_adata.obs["cell_type"] = cell_type
            pbs.append(rep_adata)

        return sc.concat(pbs)
    
    def generate_pseudo_dbulk_psuedo_rep(self, condition_key: str, cell_type: str, raw_data_layer:str = None, cols_to_keep:Sequence[str]=None) -> AnnData:
        """
        Generates pseudo bulk data with pseudo replicates for a specific cell type across different samples.

        Args:
            condition_key (str): Key in the AnnData object's metadata used to define the condition or grouping.
            cell_type (str): The cell type for which to generate pseudo bulk data with pseudo replicates.
            raw_data_layer (str, optional): Specifies the layer of AnnData object to be used. If None, uses the default layer.

        Returns:
            AnnData: An AnnData object containing the aggregated pseudo bulk data with pseudo replicates.
        """
        pbs = []
        cell_subset = self.adata[self.adata.obs[self.cell_id_key] == cell_type].copy()
        samples = self.adata.obs[self.sample_key].unique()

        if raw_data_layer:
            cell_subset.X = cell_subset.layers["raw"].A

        for sample in samples:
            sample_cell_subset = cell_subset[cell_subset.obs[self.sample_key] == sample]

            if sample_cell_subset.n_obs < self.min_number_cells:
                self.samples_rm.add(sample)
                continue
    
            indices = list(sample_cell_subset.obs.index)
            random.shuffle(indices)
            indices = np.array_split(np.array(indices), self.reps_per_sample)
            
            for i, sample_idx in enumerate(indices):
                
                rep_num = i + 1
                X = np.array([sample_cell_subset[sample_idx].X.sum(axis=0)])
                var = sample_cell_subset.var[[]]

                obs = None
                if cols_to_keep:
                    obs = sample_cell_subset.obs[cols_to_keep].iloc[0]
                    obs_df = pd.DataFrame(obs).T
                    rep_adata = AnnData(X=X, var=var, obs=obs_df)
                else:
                    rep_adata = AnnData(X=X, var=var)
            

                rep_adata.obs_names = [f"{sample}_{rep_num}"]
                rep_adata.obs[condition_key] = sample_cell_subset.obs[condition_key].iloc[0]
                rep_adata.obs["sample"] = [sample]
                rep_adata.obs["replicate"] = rep_num
                rep_adata.obs["cell_type"] = cell_type
                rep_adata.obs["num_cells"] = len(sample_idx)
                pbs.append(rep_adata)


        
        return sc.concat(pbs)



