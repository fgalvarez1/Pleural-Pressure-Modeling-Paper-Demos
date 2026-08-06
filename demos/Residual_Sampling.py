import torch

from funtd import P1DesignMatrix

class ResidualSampling:
    def __init__(self, X_observed, U_s, U_t, Phi, C, variables_combo_selected, nodes_dict):
        self.X_observed = torch.tensor(X_observed, dtype=torch.float32)
        self.U_s = torch.tensor(U_s, dtype=torch.float32)
        self.U_t = torch.tensor(U_t, dtype=torch.float32)
        self.Phi = Phi
        self.C = C
        self.variables_combo_selected = variables_combo_selected
        self.nodes_dict = nodes_dict

        self.mvn, self.Residuals = self.get_MultiVarNormal_for_sampling_residuals_using_FTD_modes()



    def get_Phi_new(self, var_dict_new):
        Phi_lst = []
        for var in self.variables_combo_selected:
            if var.count("-") == 1:
                var1, var2 = var.split("-")

                var1_values = var_dict_new[var1]
                var2_values = var_dict_new[var2]
                if var2 == "Position":
                    mask_SUP2 = (torch.tensor(var2_values) == 0)
                    mask_PRO1 = (torch.tensor(var2_values) == 1)

                    var1_SUP2 = torch.tensor(var1_values)[mask_SUP2]
                    var1_PRO1 = torch.tensor(var1_values)[mask_PRO1]

                    var1_SUP2_nodes  = self.nodes_dict[f"{var1}-SUP2"] 
                    var1_PRO1_nodes  = self.nodes_dict[f"{var1}-PRO1"]

                    Phi_var1_SUP2 = P1DesignMatrix(var1_SUP2_nodes).build(var1_values) * mask_SUP2.float().unsqueeze(1)
                    Phi_var1_PRO1 = P1DesignMatrix(var1_PRO1_nodes).build(var1_values) * mask_PRO1.float().unsqueeze(1)

                    Phi_lst.append(Phi_var1_SUP2)
                    Phi_lst.append(Phi_var1_PRO1)

                else:
                    raise ValueError("Not implemented yet")
            
            else:
                raise ValueError("Not implemented yet")
            
        Phi_new = torch.cat(Phi_lst, dim=1)
        return Phi_new



    def get_MultiVarNormal_for_sampling_residuals_using_FTD_modes(self):
        # Sample from the same modes of the FTD, using residuals
        # Residuals
        V_s = self.U_s.t() @ self.U_s
        V_t = self.U_t.t() @ self.U_t
        V   = V_s * V_t
        M_p = torch.einsum('stp,sr,tr->pr', self.X_observed, self.U_s, self.U_t)
        U_p_target = M_p @ torch.linalg.inv(V)

        Residuals = (U_p_target - self.Phi @ self.C)
        
        # corr_matrix = torch.corrcoef(Pr_pat_parafac_tensor.T)
        # print(torch.round(corr_matrix, decimals=2))
        
        mu  = torch.mean(Residuals, dim=0)
        cov = torch.cov(Residuals.T)

        mvn = torch.distributions.multivariate_normal.MultivariateNormal(loc=mu, covariance_matrix=cov)

        return mvn, Residuals



    def generate_synthetic_field(self, phi_new):
        P_data_based = torch.einsum('sr,tr,r->st', self.U_s, self.U_t, (phi_new @ self.C).flatten())
        P_sample     = torch.einsum('sr,tr,r->st', self.U_s, self.U_t, self.mvn.sample())
        P_synthetic  = P_data_based + P_sample
        return P_synthetic