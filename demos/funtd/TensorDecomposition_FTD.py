import torch
import torch.nn.functional as F

from .DesignMatrix_P1 import P1DesignMatrix

class FTD:
    def __init__(self, variables_names, variables_dict, K):
        self.variables_names = variables_names
        self.variables_dict = variables_dict
        self.K = K



    def generate_adaptive_nodes(self, var_values, min_patients_per_node=3, max_nodes=10, min_nodes=2):
        """
        Generates an adaptive finite element mesh (nodes) based on data quantiles.
        
        Args:
            var_values (torch.Tensor or list): The numerical data for the specific subgroup.
            min_patients_per_node (int): The statistical evidence required per node.
            max_nodes (int): The absolute ceiling for complexity.
            min_nodes (int): The absolute floor for complexity.
            
        Returns:
            torch.Tensor: A 1D tensor of strictly unique, sorted nodes.
        """
        # Ensure input is a float tensor
        if not isinstance(var_values, torch.Tensor):
            var_values = torch.tensor(var_values)
        var_values = var_values.float()
        
        N = len(var_values)
        assert N > 0, "The input data must contain at least one value."
            
        # Calculate number of nodes
        K = max(min_nodes, min(max_nodes, N // min_patients_per_node))
        
        # Adaptive Spacing: Create evenly spaced percentiles
        q_vals = torch.linspace(0.0, 1.0, steps=K)
        
        # Map percentiles to actual clinical values
        nodes = torch.quantile(var_values, q_vals)
        
        # Remove duplicates caused by tied data
        nodes = torch.unique(nodes)
        
        assert len(nodes) > 1, "All patients have the same value, cannot create a P1 mesh with a single node."
            
        return nodes



    def build_global_design_matrix(self):
        Phi_lst    = []
        nodes_dict = {}

        for var in self.variables_names:
            if var.count("-") == 2:
                var1, var2, var3 = var.split("-") # var2 and var3 are the categorical variables
                # Design matrix (P x K) for numerical variable
                var1_values = self.variables_dict[var1]
            
                # Get the full One-Hot status matrix (P x 2) for position, (P x 3) for group, etc.
                var2_values = self.variables_dict[var2]
                var3_values = self.variables_dict[var3]
                if var2 == "Position" and var3 == "Groupe":
                    # Position-Groupe
                    num_classes = len(set(var3_values))

                    mask_SUP2_healthy = (torch.tensor(var2_values) == 0) & (torch.tensor(var3_values) == 0)
                    mask_PRO1_healthy = (torch.tensor(var2_values) == 1) & (torch.tensor(var3_values) == 0)

                    var1_SUP2_healthy = torch.tensor(var1_values)[mask_SUP2_healthy]
                    var1_PRO1_healthy = torch.tensor(var1_values)[mask_PRO1_healthy]

                    var1_SUP2_healthy_nodes  = self.generate_adaptive_nodes(var1_SUP2_healthy, max_nodes=self.K)
                    var1_PRO1_healthy_nodes  = self.generate_adaptive_nodes(var1_PRO1_healthy, max_nodes=self.K)

                    nodes_dict[f"{var1}-SUP2-healthy"] = var1_SUP2_healthy_nodes
                    nodes_dict[f"{var1}-PRO1-healthy"] = var1_PRO1_healthy_nodes

                    Phi_var1_SUP2_healthy = P1DesignMatrix(var1_SUP2_healthy_nodes).build(var1_values) * mask_SUP2_healthy.float().unsqueeze(1)
                    Phi_var1_PRO1_healthy = P1DesignMatrix(var1_PRO1_healthy_nodes).build(var1_values) * mask_PRO1_healthy.float().unsqueeze(1)

                    Phi_lst.append(Phi_var1_SUP2_healthy)
                    Phi_lst.append(Phi_var1_PRO1_healthy)

                    if num_classes > 1:
                        mask_SUP2_asthma = (torch.tensor(var2_values) == 0) & (torch.tensor(var3_values) == 1)
                        mask_PRO1_asthma = (torch.tensor(var2_values) == 1) & (torch.tensor(var3_values) == 1)

                        var1_SUP2_asthma = torch.tensor(var1_values)[mask_SUP2_asthma]
                        var1_PRO1_asthma = torch.tensor(var1_values)[mask_PRO1_asthma]

                        var1_SUP2_asthma_nodes  = self.generate_adaptive_nodes(var1_SUP2_asthma, max_nodes=self.K)
                        var1_PRO1_asthma_nodes  = self.generate_adaptive_nodes(var1_PRO1_asthma, max_nodes=self.K)

                        nodes_dict[f"{var1}-SUP2-asthma"] = var1_SUP2_asthma_nodes
                        nodes_dict[f"{var1}-PRO1-asthma"] = var1_PRO1_asthma_nodes

                        Phi_var1_SUP2_asthma = P1DesignMatrix(var1_SUP2_asthma_nodes).build(var1_values) * mask_SUP2_asthma.float().unsqueeze(1)
                        Phi_var1_PRO1_asthma = P1DesignMatrix(var1_PRO1_asthma_nodes).build(var1_values) * mask_PRO1_asthma.float().unsqueeze(1)

                        Phi_lst.append(Phi_var1_SUP2_asthma)
                        Phi_lst.append(Phi_var1_PRO1_asthma)

                    if num_classes > 2:
                        mask_SUP2_copd = (torch.tensor(var2_values) == 0) & (torch.tensor(var3_values) == 2)
                        mask_PRO1_copd = (torch.tensor(var2_values) == 1) & (torch.tensor(var3_values) == 2)

                        var1_SUP2_copd = torch.tensor(var1_values)[mask_SUP2_copd]
                        var1_PRO1_copd = torch.tensor(var1_values)[mask_PRO1_copd]

                        var1_SUP2_copd_nodes  = self.generate_adaptive_nodes(var1_SUP2_copd, max_nodes=self.K)
                        var1_PRO1_copd_nodes  = self.generate_adaptive_nodes(var1_PRO1_copd, max_nodes=self.K)

                        nodes_dict[f"{var1}-SUP2-copd"] = var1_SUP2_copd_nodes
                        nodes_dict[f"{var1}-PRO1-copd"] = var1_PRO1_copd_nodes

                        Phi_var1_SUP2_copd = P1DesignMatrix(var1_SUP2_copd_nodes).build(var1_values) * mask_SUP2_copd.float().unsqueeze(1)
                        Phi_var1_PRO1_copd = P1DesignMatrix(var1_PRO1_copd_nodes).build(var1_values) * mask_PRO1_copd.float().unsqueeze(1)

                        Phi_lst.append(Phi_var1_SUP2_copd)
                        Phi_lst.append(Phi_var1_PRO1_copd)

                else:
                    raise ValueError(f"Combination of categorical variables still not implemented: {var2}-{var3}")

            elif var.count("-") == 1:
                var1, var2 = var.split("-") # var2 is the categorical variable
                
                if var1 in ["Position", "Groupe", "Sex"]: # var1 and var2 are categorical
                    var1_values  = self.variables_dict[var1]
                    num1_classes = len(set(var1_values))

                    var2_values  = self.variables_dict[var2]
                    num2_classes = len(set(var2_values))
                    
                    interaction_ids = (torch.tensor(var1_values) * num2_classes) + torch.tensor(var2_values)

                    Phi_var = F.one_hot(interaction_ids, num_classes=num1_classes * num2_classes).float()
                    Phi_lst.append(Phi_var)

                    nodes_dict[var] = [f"{i}-{j}" for i in range(num1_classes) for j in range(num2_classes)] # list(set(var_values))
                
                else:
                    # Design matrix (P x K) for numerical variable
                    var1_values = self.variables_dict[var1]
                    var2_values = self.variables_dict[var2]
                    if var2 == "Position":
                        mask_SUP2 = (torch.tensor(var2_values) == 0)
                        mask_PRO1 = (torch.tensor(var2_values) == 1)

                        var1_SUP2 = torch.tensor(var1_values)[mask_SUP2]
                        var1_PRO1 = torch.tensor(var1_values)[mask_PRO1]

                        var1_SUP2_nodes  = self.generate_adaptive_nodes(var1_SUP2, max_nodes=self.K)
                        var1_PRO1_nodes  = self.generate_adaptive_nodes(var1_PRO1, max_nodes=self.K)

                        nodes_dict[f"{var1}-SUP2"] = var1_SUP2_nodes
                        nodes_dict[f"{var1}-PRO1"] = var1_PRO1_nodes

                        Phi_var1_SUP2 = P1DesignMatrix(var1_SUP2_nodes).build(var1_values) * mask_SUP2.float().unsqueeze(1)
                        Phi_var1_PRO1 = P1DesignMatrix(var1_PRO1_nodes).build(var1_values) * mask_PRO1.float().unsqueeze(1)

                        Phi_lst.append(Phi_var1_SUP2)
                        Phi_lst.append(Phi_var1_PRO1)

                    elif var2 == "Groupe":
                        num_classes = len(set(var2_values))

                        mask_healthy                  = (torch.tensor(var2_values) == 0)
                        var1_healthy                  = torch.tensor(var1_values)[mask_healthy]
                        var1_healthy_nodes            = self.generate_adaptive_nodes(var1_healthy, max_nodes=self.K)
                        nodes_dict[f"{var1}-healthy"] = var1_healthy_nodes
                        Phi_var1_healthy              = P1DesignMatrix(var1_healthy_nodes).build(var1_values) * mask_healthy.float().unsqueeze(1)
                        Phi_lst.append(Phi_var1_healthy)

                        if num_classes > 1:
                            mask_asthma                  = (torch.tensor(var2_values) == 1)
                            var1_asthma                  = torch.tensor(var1_values)[mask_asthma]
                            var1_asthma_nodes            = self.generate_adaptive_nodes(var1_asthma, max_nodes=self.K)
                            nodes_dict[f"{var1}-asthma"] = var1_asthma_nodes
                            Phi_var1_asthma              = P1DesignMatrix(var1_asthma_nodes).build(var1_values) * mask_asthma.float().unsqueeze(1)
                            Phi_lst.append(Phi_var1_asthma)

                        if num_classes > 2:
                            mask_copd                  = (torch.tensor(var2_values) == 2)
                            var1_copd                  = torch.tensor(var1_values)[mask_copd]
                            var1_copd_nodes            = self.generate_adaptive_nodes(var1_copd, max_nodes=self.K)
                            nodes_dict[f"{var1}-copd"] = var1_copd_nodes
                            Phi_var1_copd              = P1DesignMatrix(var1_copd_nodes).build(var1_values) * mask_copd.float().unsqueeze(1)
                            Phi_lst.append(Phi_var1_copd)

                        # NOTE For now only cases: healthy, healthy+asthma, healthy+asthma+copd. We can add more if needed.

                    elif var2 == "Sex":
                        mask_M = (torch.tensor(var2_values) == 0)
                        mask_F = (torch.tensor(var2_values) == 1)

                        var1_M = torch.tensor(var1_values)[mask_M]
                        var1_F = torch.tensor(var1_values)[mask_F]

                        var1_M_nodes  = self.generate_adaptive_nodes(var1_M, max_nodes=self.K)
                        var1_F_nodes  = self.generate_adaptive_nodes(var1_F, max_nodes=self.K)

                        nodes_dict[f"{var1}-M"] = var1_M_nodes
                        nodes_dict[f"{var1}-F"] = var1_F_nodes

                        Phi_var1_M = P1DesignMatrix(var1_M_nodes).build(var1_values) * mask_M.float().unsqueeze(1)
                        Phi_var1_F = P1DesignMatrix(var1_F_nodes).build(var1_values) * mask_F.float().unsqueeze(1)

                        Phi_lst.append(Phi_var1_M)
                        Phi_lst.append(Phi_var1_F)
                    
                    else:
                        raise ValueError(f"Categorical variable still not implemented: {var2}")

            elif var in ["Position", "Groupe", "Sex"]:
                var_values  = self.variables_dict[var]
                num_classes = len(set(var_values))
                Phi_var     = F.one_hot(torch.tensor(var_values), num_classes=num_classes).float()
                Phi_lst.append(Phi_var)

                nodes_dict[var] = list(set(var_values))
                assert self.variables_names == [var], "This case is intended for when we only want to include one categorical variable."
            
            else:
                var_values = self.variables_dict[var]
                var_nodes  = self.generate_adaptive_nodes(var_values, max_nodes=self.K)
                Phi_var    = P1DesignMatrix(var_nodes).build(var_values)
                Phi_lst.append(Phi_var)
                
                nodes_dict[var] = var_nodes
        
        Phi = torch.cat(Phi_lst, dim=1)
        return Phi, nodes_dict



    def fit(self, X_observed, Phi, U_s_init, U_t_init, R):
        # ALS algorithm
        X_observed = torch.tensor(X_observed, dtype=torch.float32)
        U_s        = torch.tensor(U_s_init, dtype=torch.float32)
        U_t        = torch.tensor(U_t_init, dtype=torch.float32)

        X_observed_power = torch.mean(X_observed**2)
        
        # U_s = U_s_init.clone()
        # U_t = U_t_init.clone()

        C = torch.randn(Phi.shape[1], R)

        epochs = 5

        print("ALS Algorithm")
        for epoch in range(epochs):
            
            # ==========================================
            # PHASE 1: EXACT UPDATE FOR C (Age Mode)
            # ==========================================
            # 1. Compute Gram matrices
            V_s = U_s.t() @ U_s
            V_t = U_t.t() @ U_t
            V = V_s * V_t  # Hadamard product
            
            # 2. MTTKRP: Project X onto Space and Time
            M_p = torch.einsum('stp,sr,tr->pr', X_observed, U_s, U_t)
            
            # 3. Unconstrained target: U_p_target = M_p @ V^{-1}
            U_p_target = M_p @ torch.linalg.inv(V)
            
            # # 4. Constrained update for C
            # C = torch.linalg.lstsq(Phi, U_p_target).solution

            # 4. A bit of Tikhonov regularization for numerical stability
            I = torch.eye(Phi.shape[1])
            lambda_reg = 1e-6

            # Compute Ridge Regression explicitly
            Phi_T_Phi = Phi.t() @ Phi
            Phi_T_U   = Phi.t() @ U_p_target

            # Solve the regularized system: (Phi^T Phi + lambda*I) * C = Phi^T U
            C = torch.linalg.solve(Phi_T_Phi + lambda_reg * I, Phi_T_U)
            
            # 5. Rebuild the current U_p for the next phases
            norms = torch.norm(C, p=2, dim=0, keepdim=True)
            C     = C / norms
            U_p   = Phi @ C
            U_t   = U_t * norms
            
            # ==========================================
            # PHASE 2: EXACT UPDATE FOR U_s (Spatial Mode)
            # ==========================================
            V_t = U_t.t() @ U_t
            V_p = U_p.t() @ U_p
            V = V_t * V_p
            
            M_s = torch.einsum('stp,tr,pr->sr', X_observed, U_t, U_p)
            U_s = M_s @ torch.linalg.inv(V)
            
            # Normalize U_s to avoid scale explosion, push scale to C
            norms = torch.norm(U_s, p=2, dim=0, keepdim=True)
            U_s = U_s / norms
            U_t = U_t * norms
            # C = C * norms  
            # U_p = Phi @ C  # Rebuild U_p with new scale
            
            # ==========================================
            # PHASE 3: EXACT UPDATE FOR U_t (Temporal Mode)
            # ==========================================
            V_s = U_s.t() @ U_s
            V_p = U_p.t() @ U_p
            V = V_s * V_p
            
            M_t = torch.einsum('stp,sr,pr->tr', X_observed, U_s, U_p)
            U_t = M_t @ torch.linalg.inv(V)
            
            # Normalize U_t, push scale to C
            # norms = torch.norm(U_t, p=2, dim=0, keepdim=True)
            # U_t = U_t / norms
            # C = C * norms  
            # U_p = Phi @ C  # Rebuild U_p with new scale

            # ==========================================
            # EVALUATE LOSS
            # ==========================================
            # Reconstruct tensor to check convergence
            X_hat = torch.einsum('sr,tr,pr->stp', U_s, U_t, U_p)
            loss  = torch.nn.functional.mse_loss(X_hat, X_observed) / X_observed_power
            
            print(f"Epoch {epoch+1:02d} | Normalized MSE: {loss.item():.6f}")
        
        # Global R2 Score
        var_observed = torch.var(X_observed, unbiased=False)
        r2_score     = 1. - (torch.nn.functional.mse_loss(X_hat, X_observed) / var_observed)
        print (f"Final R2 Score: {r2_score.item():.4f}")
        
        ftd_matrices = {"U_s" : U_s.detach().cpu().numpy(),
                             "U_t" : U_t.detach().cpu().numpy(),
                             "C"   : C.detach().cpu().numpy()}
        
        return ftd_matrices, loss.item()