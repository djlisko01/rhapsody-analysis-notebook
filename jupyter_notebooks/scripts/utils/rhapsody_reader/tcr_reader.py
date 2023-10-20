import pandas as pd
import numpy as np
import anndata
import 


class TCR_Reader:
    def __init__(self):
        self.tcr_anndata = anndata.AnnData()

    def VDJ_perCell_loader(self, file_path: str) -> None:
       df =  pd.read_csv(file_path, comment="#")

    
    def _set_receptor_subtype(self, series: pd.Series) -> pd.Series:
        pass
