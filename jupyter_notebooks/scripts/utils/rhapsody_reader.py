import gzip
import shutil
import warnings

import mudata
import numpy as np
from anndata import AnnData
import pandas as pd
from enum import Enum
import re
from typing import List, Union, Tuple
import os
import scanpy as sc
from scipy import io

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from muon import MuData
    import muon as mu

SPLIT_BY = "|pAbO"  # BD specific


class RhapsodyReader:
    def __init__(self):
        self.mudata = None
        self.genes_per_sample = {}
        self.abt_per_sample = {}

    def read_from_csv(self, file_path: os.PathLike) -> None:
        adata = sc.read_csv(file_path,
                            first_column_names=True,
                            dtype='float32'
                            )
        abt_df, rna_df = self.split_data_by(adata=adata, split_by=SPLIT_BY)
        self.create_mudata(abt_df, rna_df)

    def split_data_by(self, adata: AnnData, split_by: str = SPLIT_BY) -> tuple:
        abt_df = self.get_abt_data(adata, split_by)
        rna_df = self.get_rna_data(adata, split_by)
        return abt_df, rna_df

    def read_multi_csv(self, root_dir: str, regex: str, comment: str = "#", index_col: Union[str, int] = 0) -> None:
        file_paths, samples = self.search_files(root_dir, regex)
        df_list = []
        set_genes = set()
        sets_genes = []

        # TODO: SPEED UP THIS PROCESS?
        for file_path, sample in zip(file_paths, samples):
            df = pd.read_csv(file_path, comment=comment)
            # Rename set ran rename index columns
            index_col = self.__get_indices_name(df, index_col=index_col)
            df.set_index(index_col, inplace=True)
            df.index = df.index.astype(str) + "_" + sample
            df_list.append(df)
            is_pbt = df.columns.str.endswith(SPLIT_BY)
            self.genes_per_sample[sample] = np.sum(~is_pbt)
            self.abt_per_sample[sample] = np.sum(is_pbt)

        df = pd.concat(df_list, verify_integrity=True)  # verify unique indices
        abt_df, rna_df = self.split_data_by(df, split_by=SPLIT_BY)  # split data
        df.fillna(0, inplace=True)
        self.create_mudata(abt_df, rna_df)

        # Just a faster way to keep track of genes and abt per sample.
        self.mudata.mod["abt"].uns["abt_per_sample"] = self.abt_per_sample
        self.mudata.mod["rna"].uns["genes_per_sample"] = self.genes_per_sample

    def create_mudata(self, abt_data, rna_data):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.mudata = MuData({"abt": abt_data, "rna": rna_data})

    def add_obs_data(self, root, regex):
        pass

    @staticmethod
    def __get_indices_name(df, index_col: Union[str, int] = 0) -> str:
        if type(index_col) is int:
            return df.columns[index_col]
        return index_col

    @staticmethod
    def search_files(directory, pattern) -> Tuple[List[str], List[str]]:
        results = []  # List to store matching file paths
        samples = []

        # Compile the regex pattern
        regex = re.compile(pattern)

        # Walk through the directory and its subdirectories
        for root, dirs, files in os.walk(directory):
            for file in files:
                file_path = os.path.join(root, file)
                if regex.search(file):
                    sample = root.split("/")[-1]  # Assumes the regex file is a child to the sample directory
                    samples.append(sample)
                    results.append(file_path)

        return results, samples

    @staticmethod
    def read_mudata(file_path):
        return mu.read_h5mu(file_path)

    def save_mudata(self, save_dir: str, file_name: str) -> None:
        file_location = f'{save_dir}/{file_name}.h5mu'
        if self.mudata:
            self.mudata.write(file_location)
            print(f"Saved to {file_location}")

    def add_cell_annotations(self, file_dir: str, comment: str = "#") -> None:
        obs_df = pd.read_csv(file_dir, comment=comment, index_col=0)
        self.mudata.obs = obs_df

    def get_mudata_obj(self) -> mudata.MuData:
        if self.mudata:
            return self.mudata.copy()

    def clean_memory(self) -> None:
        del self.mudata
        self.mudata = None

    @staticmethod
    def get_abt_data(adata: AnnData, abt_pattern: str) -> AnnData:
        bool_lst = adata.var_names.str.endswith(abt_pattern)
        return adata[:, bool_lst]

    @staticmethod
    def get_rna_data(adata: AnnData, abt_pattern: str) -> AnnData:
        bool_lst = np.invert(adata.var_names.str.endswith(abt_pattern))
        return adata[:, bool_lst]
    
    @staticmethod
    def create_10x_files(adata, dir, raw_layer):
        RhapsodyReader.create_barcodes_table(adata=adata, save_path=f"{dir}/barcodes.tsv")
        RhapsodyReader.create_features_table(adata, f"{dir}/features.tsv")
        RhapsodyReader.create_matrix(adata, f"{dir}/matrix.mtx", layer=raw_layer)
        # RhapsodyReader.gzip_files_in_directory(dir)

    
    @staticmethod
    def create_matrix(adata, save_path, layer=None):
        if layer:
            matrix = adata.layers[layer].T
        else:
            matrix = adata.X.T
            
        io.mmwrite(save_path, matrix.astype(np.int64))

    @staticmethod
    def create_barcodes_table(adata, save_path):
        with open(save_path, "w") as file:
            for item in adata.obs_names:
                file.write(item + "\n")
    
    @staticmethod
    def create_features_table(adata, save_path):
        tab_var_names = ["\t".join([x, x, "Gene Expression"]) for x in adata.var_names]

        with open(save_path, "w") as file:
            for var_name in tab_var_names:
                file.write(var_name + "\n")

    @staticmethod
    def gzip_files_in_directory(directory_path):
        for root, _, files in os.walk(directory_path):
            for file in files:
                file_path = os.path.join(root, file)
                with open(file_path, 'rb') as f_in:
                    with gzip.open(file_path + '.gz', 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                os.remove(file_path)  # remove the original file after compression

    

class FileFormats(Enum):
    CSV = "csv"
    MUDATA = "mudata"


if __name__ == "__main__":
    rp_read = RhapsodyReader()
    root_dir = "/Users/djlisko/gitrepos/rhapsody-analysis/data/00_raw_data"
    regex = r"_?_Combined_*seq\d-p\d+\w_RSEC_Mol"
    rp_read.read_multi_csv(root_dir, regex)
    # print(rp_read.mudata.mod["rna"].uns)
