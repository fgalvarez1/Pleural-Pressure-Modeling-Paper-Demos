import numpy as np
import dolfin
import meshio

class MeshDisplacementReader:
    def __init__(self, acquisition, volunteer, region, n_phases):
        self.resultsPath = f"./Results_{acquisition}/{volunteer}"
        self.region      = region
        self.n_phases    = n_phases
        
        self.mesh_exhal_file = f"{self.resultsPath}/Mesh/Mesh_{region}_00.xdmf"
        self.mesh_exhal      = self.read_mesh(self.mesh_exhal_file)

        self.load_U_dict_array_exhal_to_t()



    def read_mesh(self, file, U_move=None, save_name=None):
        mesh = dolfin.Mesh()
        with dolfin.XDMFFile(file) as xdmf:
            xdmf.read(mesh)

        if U_move is not None:
            for U_move_i in U_move:
                dolfin.ALE.move(mesh, U_move_i)
        
        if save_name is not None:
            with dolfin.XDMFFile(save_name) as xdmf_file:
                xdmf_file.write(mesh)

        return mesh



    def load_U_dict_array_exhal_to_t(self):
        folder_mesh_disp = f"{self.resultsPath}/Fields_from_image_analysis/mesh_{self.region}.xdmf"

        U_dict_array_exhal_to_t = {} # dictionary with keys as phase numbers and values as displacement arrays
        with meshio.xdmf.TimeSeriesReader(folder_mesh_disp) as reader:
            points, cells = reader.read_points_cells()
            assert reader.num_steps == self.n_phases, f"Expected {self.n_phases} phases, but found {reader.num_steps} in {folder_mesh_disp}"
            
            for k in range(reader.num_steps):
                t, point_data, cell_data = reader.read_data(k)
                assert t % 1 == 0
                U_dict_array_exhal_to_t[int(t)] = point_data['deformation_fields']
        
        self.U_dict_array_exhal_to_t = U_dict_array_exhal_to_t



    def U_array_to_Function(self, mesh, U_array):
        """
        Convert displacement array from image analysis to a dolfin Function
        """
        # Create a VectorFunctionSpace on mesh
        V = dolfin.VectorFunctionSpace(mesh, "CG", 1)
        U_exhal_to_inhal = dolfin.Function(V)

        dof_coords = V.tabulate_dof_coordinates()
        for point_index in range(dof_coords.shape[0]//3):
            dof_index = 3 * point_index

            matches        = np.all(dof_coords[dof_index] == mesh.coordinates(), axis = 1)
            ref_mesh_index = np.where(matches)
            assert len(ref_mesh_index[0]) == 1, "dof coordinates not found in mesh nodes"

            U_exhal_to_inhal.vector()[dof_index]     = U_array[ref_mesh_index, 0][0][0]
            U_exhal_to_inhal.vector()[dof_index + 1] = U_array[ref_mesh_index, 1][0][0]
            U_exhal_to_inhal.vector()[dof_index + 2] = U_array[ref_mesh_index, 2][0][0]

        return U_exhal_to_inhal
    


    def get_U_exhal_to_t(self, mesh, phase):
        """
        Get the displacement field for a given phase
        """
        U_array = self.U_dict_array_exhal_to_t[phase]
        return self.U_array_to_Function(mesh, U_array)