import os
import scipy.io as sio 

def adata_to_mtx(adata, dir_path):
    try:
        mtx = adata.X.T
        # File paths for features, obs and mtx
        barcode_path = os.path.join(dir_path, "barcodes.tsv")
        feature_path = os.path.join(dir_path, "features.tsv")
        mtx_path =os.path.join(dir_path, "sparse_matrix.mtx")
            
        with open(barcode_path, "w") as f:
            barcodes = "\n".join(list(adata.obs_names))
            f.writelines(barcodes)
        
        with open(feature_path, "w") as f:
            var_names = "\n".join(list(adata.var_names))
            f.writelines(var_names)
        sio.mmwrite(mtx_path, mtx)
        print(f"File saved to {dir_path}")
        
    except Exception as err:
        print(err)
        print("Failed to export data")