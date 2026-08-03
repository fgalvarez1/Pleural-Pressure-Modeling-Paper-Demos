import torch

class P1DesignMatrix:
    def __init__(self, var_nodes):
        self.var_nodes = var_nodes



    def build(self, var_values):
        """
        Constructs the Phi matrix for P1 finite elements.
        
        var_values: 1D tensor of length P (observed values of certain variable)
        var_nodes: 1D tensor of length K (finite element node locations)
        
        Returns:
        Phi: Tensor of shape (P, K)
        """
        P = len(var_values)
        K = len(self.var_nodes)
        Phi = torch.zeros((P, K))
        
        for i, var_value in enumerate(var_values):
            # Handle edge cases (extrapolation if age is outside the nodes)
            if var_value <= self.var_nodes[0].item():
                Phi[i, 0] = 1.0
            elif var_value >= self.var_nodes[-1].item():
                Phi[i, -1] = 1.0
            else:
                # Find which interval [node_left, node_right] the age falls into
                idx = torch.searchsorted(self.var_nodes, var_value) - 1
                n_left, n_right = self.var_nodes[idx].item(), self.var_nodes[idx+1].item()
                
                # P1 (Linear) interpolation weights
                dist = n_right - n_left
                w_left = (n_right - var_value) / dist
                w_right = (var_value - n_left) / dist
                
                Phi[i, idx] = w_left
                Phi[i, idx + 1] = w_right
                
        return Phi