import numpy as np
import math
import matplotlib.pyplot as plt

class FiguresFTD:
    def __init__(self, K, var_dict, U_s, U_t, Phi, C, variables_combo_selected, nodes_dict):
        self.K = K
        self.var_dict = var_dict
        self.U_s = U_s
        self.U_t = U_t
        self.C = C
        self.Phi = Phi
        self.variables_combo_selected = variables_combo_selected
        self.nodes_dict = nodes_dict

        self.R = U_s.shape[1]



    def plot_temporal_modes(self, filename):
        font_size = 14

        # Plot temporal modes
        plt.figure()
        for ind_r in range(self.R):
            plt.plot(self.U_t[:, ind_r], label=f'Mode {ind_r+1}')

        plt.xlim(0, self.U_t.shape[0] - 1)

        plt.xlabel('Phase', fontsize=font_size)
        plt.ylabel('Temporal Modes', fontsize=font_size)
        plt.tick_params(axis='both', labelsize=font_size)
        plt.legend(fontsize=font_size)

        plt.savefig(filename, bbox_inches='tight')
        plt.show()



    def plot_variables_modes(self, filename):
        font_size = 14

        plt.rcParams.update({
            'font.size': font_size,          # General font size
            'axes.titlesize': font_size + 2,  # Subplot titles
            'axes.labelsize': font_size,      # X and Y labels
            'xtick.labelsize': font_size - 2, # X tick numbers
            'ytick.labelsize': font_size - 2, # Y tick numbers
            'legend.fontsize': font_size - 2, # Legend text
            'figure.titlesize': font_size + 4 # Main figure title
        })

        n_rows = math.ceil((self.Phi.shape[1]/self.K)/2)
        fig, axis = plt.subplots(n_rows, 2, figsize=(12, 5*n_rows),sharey=True)
        axs = axis.flatten()

        i_var  = 0
        i_node = 0
        for var in self.variables_combo_selected:
            if var.count("-") == 2:
                var1, var2, var3 = var.split("-")

                if var2 == "Position" and var3 == "Groupe":
                    nodes = self.nodes_dict[f"{var1}-SUP2-healthy"]
                    var1_replaced = var1.replace("p_transported_concatenated_norm", r"p_{\text{norm}}")
                    var1_replaced = var1_replaced.replace("DeltaV_V", r"\Delta V / V_{\rm{exhal}}")
                    var1_replaced = var1_replaced.replace("p_end_inh_mean", r"p_{\rm{i}}")
                    # var1_replaced = var1_replaced.replace('_', r'\_')
                    for r in range(self.R):
                        axs[i_var].plot(nodes, self.C[i_node:i_node+len(nodes),r], ".-",label=f'Mode {r+1}')
                    axs[i_var].set_xlabel(f"${var1_replaced}$")
                    # axs[i_var].set_ylabel(rf'$C_{{\mathrm{{{var1_replaced}}}}}$ (SUP2-healthy)')
                    axs[i_var].set_ylabel(rf'${var1_replaced}$ modes (SUP2-healthy)')
                    axs[i_var].legend()
                    i_var  += 1
                    i_node += len(nodes)

                    nodes = self.nodes_dict[f"{var1}-PRO1-healthy"]
                    for r in range(self.R):
                        axs[i_var].plot(nodes, self.C[i_node:i_node+len(nodes),r], ".-",label=f'Mode {r+1}')
                    axs[i_var].set_xlabel(f"${var1_replaced}$")
                    # axs[i_var].set_ylabel(rf'$C_{{\mathrm{{{var1_replaced}}}}}$ (PRO1-healthy)')
                    axs[i_var].set_ylabel(rf'${var1_replaced}$ modes (PRO1-healthy)')
                    axs[i_var].legend()
                    i_var  += 1
                    i_node += len(nodes)

                    var3_values = self.var_dict[var3]
                    num_classes = len(set(var3_values))

                    if num_classes > 1:
                        nodes = self.nodes_dict[f"{var1}-SUP2-asthma"]
                        for r in range(self.R):
                            axs[i_var].plot(nodes, self.C[i_node:i_node+len(nodes),r], ".-",label=f'Mode {r+1}')
                        axs[i_var].set_xlabel(f"${var1_replaced}$")
                        # axs[i_var].set_ylabel(rf'$C_{{\mathrm{{{var1_replaced}}}}}$ (SUP2-asthma)')
                        axs[i_var].set_ylabel(rf'${var1_replaced}$ modes (SUP2-asthma)')
                        axs[i_var].legend()
                        i_var += 1
                        i_node += len(nodes)

                        nodes = self.nodes_dict[f"{var1}-PRO1-asthma"]
                        for r in range(self.R):
                            axs[i_var].plot(nodes, self.C[i_node:i_node+len(nodes),r], ".-",label=f'Mode {r+1}')
                        axs[i_var].set_xlabel(f"${var1_replaced}$")
                        # axs[i_var].set_ylabel(rf'$C_{{\mathrm{{{var1_replaced}}}}}$ (PRO1-asthma)')
                        axs[i_var].set_ylabel(rf'${var1_replaced}$ modes (PRO1-asthma)')
                        axs[i_var].legend()
                        i_var  += 1
                        i_node += len(nodes)

                    if num_classes > 2:
                        nodes = self.nodes_dict[f"{var1}-SUP2-copd"]
                        for r in range(self.R):
                            axs[i_var].plot(nodes, self.C[i_node:i_node+len(nodes),r], ".-",label=f'Mode {r+1}')
                        axs[i_var].set_xlabel(f"${var1_replaced}$")
                        # axs[i_var].set_ylabel(rf'$C_{{\mathrm{{{var1_replaced}}}}}$ (SUP2-COPD)')
                        axs[i_var].set_ylabel(rf'${var1_replaced}$ modes (SUP2-COPD)')
                        axs[i_var].legend()
                        i_var  += 1
                        i_node += len(nodes)

                        nodes = self.nodes_dict[f"{var1}-PRO1-copd"]
                        for r in range(self.R):
                            axs[i_var].plot(nodes, self.C[i_node:i_node+len(nodes),r], ".-",label=f'Mode {r+1}')
                        axs[i_var].set_xlabel(f"${var1_replaced}$")
                        # axs[i_var].set_ylabel(rf'$C_{{\mathrm{{{var1_replaced}}}}}$ (PRO1-COPD)')
                        axs[i_var].set_ylabel(rf'${var1_replaced}$ modes (PRO1-COPD)')
                        axs[i_var].legend()
                        i_var  += 1
                        i_node += len(nodes)

            elif var.count("-") == 1:
                var1, var2 = var.split("-")
                var1_replaced = var1.replace("p_transported_concatenated_norm", r"p_{\text{norm}}")
                var1_replaced = var1_replaced.replace("DeltaV_V", r"\Delta V / V_{\rm{exhal}}")
                var1_replaced = var1_replaced.replace("p_end_inh_mean", r"p_{\rm{i}}")
                # var1_replaced = var1_replaced.replace('_', r'\_')

                if var2 == "Position":
                    # Supine
                    nodes = self.nodes_dict[f"{var1}-SUP2"]
                    for r in range(self.R):
                        axs[i_var].plot(nodes, self.C[i_node:i_node+len(nodes),r], ".-",label=f'Mode {r+1}')

                    mask_SUP2 = np.array(self.var_dict[var2]) == 0
                    custom_vertical_line = [(0, 0), (0, 1)]
                    axs[i_var].scatter(np.array(self.var_dict[var1])[mask_SUP2], [0] * len(np.array(self.var_dict[var1])[mask_SUP2]), 
                            marker=custom_vertical_line, 
                            s=150,                       # Adjust 's' to change the height of the lines
                            linewidths=1.5,              # Adjust thickness of the lines
                            color='tab:grey', 
                            # label='var1 distribution', 
                            transform=axs[i_var].get_xaxis_transform(),
                            clip_on=False)               # clip_on=False ensures dots don't get cut off
                    
                    if var1 == "p_end_inh_mean":
                        axs[i_var].set_xlabel(f"${var1_replaced}$ (kPa)", fontsize=font_size)
                    else:
                        axs[i_var].set_xlabel(f"${var1_replaced}$", fontsize=font_size)
                    # axs[i_var].set_ylabel(rf'$C_{{\mathrm{{{var1_replaced}}}}}$ (SUP2)')
                    axs[i_var].set_ylabel(rf'${var1_replaced}$ modes (supine)', fontsize=font_size)
                    axs[i_var].legend(fontsize=font_size-2)
                    axs[i_var].tick_params(axis="both", labelsize=font_size)
                    # axs[i_var].set_xlim(torch.min(nodes), torch.max(nodes))
                    i_var  += 1
                    i_node += len(nodes)

                    # Prone
                    nodes = self.nodes_dict[f"{var1}-PRO1"]
                    for r in range(self.R):
                        axs[i_var].plot(nodes, self.C[i_node:i_node+len(nodes),r], ".-",label=f'Mode {r+1}')

                    mask_PRO1 = np.array(self.var_dict[var2]) == 1
                    custom_vertical_line = [(0, 0), (0, 1)]
                    axs[i_var].scatter(np.array(self.var_dict[var1])[mask_PRO1], [0] * len(np.array(self.var_dict[var1])[mask_PRO1]), 
                            marker=custom_vertical_line, 
                            s=150,                       # Adjust 's' to change the height of the lines
                            linewidths=1.5,              # Adjust thickness of the lines
                            color='tab:grey', 
                            # label='var1 distribution', 
                            transform=axs[i_var].get_xaxis_transform(),
                            clip_on=False)               # clip_on=False ensures dots don't get cut off

                    if var1 == "p_end_inh_mean":
                        axs[i_var].set_xlabel(f"${var1_replaced}$ (kPa)", fontsize=font_size)
                    else:
                        axs[i_var].set_xlabel(f"${var1_replaced}$", fontsize=font_size)
                    # axs[i_var].set_ylabel(rf'$C_{{\mathrm{{{var1_replaced}}}}}$ (PRO1)')
                    axs[i_var].set_ylabel(rf'${var1_replaced}$ modes (prone)', fontsize=font_size)
                    axs[i_var].legend(fontsize=font_size-2)
                    axs[i_var].tick_params(axis="both", labelsize=font_size)
                    # axs[i_var].set_xlim(torch.min(nodes), torch.max(nodes))
                    i_var  += 1
                    i_node += len(nodes)

                elif var2 == "Groupe":
                    var2_values = self.var_dict[var2]
                    num_classes = len(set(var2_values))

                    nodes = self.nodes_dict[f"{var1}-healthy"]
                    for r in range(self.R):
                        axs[i_var].plot(nodes, self.C[i_node:i_node+len(nodes),r], ".-",label=f'Mode {r+1}')
                    axs[i_var].set_xlabel(f"${var1_replaced}$")
                    # axs[i_var].set_ylabel(rf'$C_{{\mathrm{{{var1_replaced}}}}}$ (healthy)')
                    axs[i_var].set_ylabel(rf'{var1_replaced} modes (healthy)')
                    axs[i_var].legend()
                    i_var  += 1
                    i_node += len(nodes)

                    if num_classes > 1:
                        nodes = self.nodes_dict[f"{var1}-asthma"]
                        for r in range(self.R):
                            axs[i_var].plot(nodes, self.C[i_node:i_node+len(nodes),r], ".-",label=f'Mode {r+1}')
                        axs[i_var].set_xlabel(f"${var1_replaced}$")
                        # axs[i_var].set_ylabel(rf'$C_{{\mathrm{{{var1_replaced}}}}}$ (asthma)')
                        axs[i_var].set_ylabel(rf'{var1_replaced} modes (asthma)')
                        axs[i_var].legend()
                        i_var  += 1
                        i_node += len(nodes)

                    if num_classes > 2:
                        nodes = self.nodes_dict[f"{var1}-copd"]
                        for r in range(self.R):
                            axs[i_var].plot(nodes, self.C[i_node:i_node+len(nodes),r], ".-",label=f'Mode {r+1}')
                        axs[i_var].set_xlabel(f"${var1_replaced}$")
                        # axs[i_var].set_ylabel(rf'$C_{{\mathrm{{{var1_replaced}}}}}$ (COPD)')
                        axs[i_var].set_ylabel(rf'{var1_replaced} modes (COPD)')
                        axs[i_var].legend()
                        i_var  += 1
                        i_node += len(nodes)

            elif var in ["Position", "Groupe", "Sex"]:
                for r in range(self.R):
                    axs[i_var].plot(self.nodes_dict[var], self.C[:,r], ".", label=f'Mode {r+1}')
                axs[i_var].set_xlabel(f"${var}$")
                # axs[i_var].set_ylabel(rf'$C_{{\mathrm{{{var}}}}}$')
                axs[i_var].set_ylabel(rf'${var}$ modes')
                axs[i_var].legend()
                axs[i_var].set_xticks(self.nodes_dict[var])
                i_var += 1

            else:
                var_replaced = var.replace("p_transported_concatenated_norm", r"$p_{\text{norm}}$")
                var_replaced = var_replaced.replace("DeltaV_V", r"$\Delta V / V_{\rm{exhal}}$")
                var_replaced = var_replaced.replace("p_end_inh_mean", r"p_{\rm{i}}")
                # var_replaced = var_replaced.replace('_', r'\_')
                nodes = self.nodes_dict[var]
                for r in range(self.R):
                    axs[i_var].plot(nodes, self.C[i_node:i_node+len(nodes),r], ".-", label=f'Mode {r+1}')
                axs[i_var].set_xlabel(f"${var_replaced}$")
                # axs[i_var].set_ylabel(rf'$C_{{\mathrm{{{var_replaced}}}}}$')
                axs[i_var].set_ylabel(rf'${var_replaced}$ modes')
                axs[i_var].legend()
                i_var  += 1
                i_node += len(nodes)

        # Delete unused subplots
        for idx in range(i_var, len(axs)):
            fig.delaxes(axs[idx])

        plt.tight_layout()
        plt.savefig(filename, bbox_inches='tight')
        plt.show()