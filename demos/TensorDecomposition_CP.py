import numpy as np

import torch

import tensorly as tl
from tensorly.decomposition import parafac

class CP:
    def __init__(self,):
        pass



    def get_normalized_parafac_decomposition(self, data_tensor, rank):
        # CP for initial guess
        weights, factors = parafac(data_tensor, rank=rank, verbose=False)
        spatial_parafac = factors[0] # n_p x n_modes
        temp_parafac    = factors[1]
        pat_parafac     = factors[2]

        assert set(weights) == {1}

        norm_spatial = np.linalg.norm(spatial_parafac, axis=0, keepdims=True)
        # norm_temp    = np.linalg.norm(temp_parafac, axis=0, keepdims=True)
        norm_pat     = np.linalg.norm(pat_parafac, axis=0, keepdims=True)

        normalized_spatial_parafac = spatial_parafac / norm_spatial
        # normalized_temp_parafac    = temp_parafac / norm_temp
        scaled_temp_parafac        = temp_parafac * norm_spatial * norm_pat
        normalized_pat_parafac     = pat_parafac / norm_pat

        X_parafac   = tl.cp_to_tensor((weights, factors))
        mse_parafac = torch.nn.functional.mse_loss(torch.tensor(X_parafac), torch.tensor(data_tensor)).item() / torch.mean(torch.tensor(data_tensor)**2)
        print(f"MSE with parafac: {mse_parafac:.4f}")

        normalized_parafac = {"normalized_spatial_parafac" : normalized_spatial_parafac,
                            "scaled_temp_parafac"        : scaled_temp_parafac,
                            "normalized_pat_parafac"     : normalized_pat_parafac}
        
        return normalized_parafac, mse_parafac