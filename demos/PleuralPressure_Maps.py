import numpy as np

import dolfin
import dolfin_mech as dmech

from Reader_MeshDisplacement import MeshDisplacementReader
from Reader_xdmf import xdmfReader

class PleuralPressureMaps:
    def __init__(self, acquisition, volunteer, region, alpha_, gravity_, pe, params, mesh_unloaded, U_exhal_to_unloaded, MDR: MeshDisplacementReader,  Phis0_unloaded_array=None, assume_pf_0=True):
        self.resultsPath         = f"./Results_{acquisition}/{volunteer}"
        self.pleuralpressurePath = f"{self.resultsPath}/Pleural_pressure_estimation"
        self.region              = region
        self.alpha_              = alpha_
        self.gravity_            = gravity_
        self.pe                  = pe
        self.params              = params
        self.mesh_unloaded       = mesh_unloaded
        self.U_exhal_to_unloaded = U_exhal_to_unloaded
        self.MDR                 = MDR

        if Phis0_unloaded_array is not None:
            self.Phis0_unloaded_array = Phis0_unloaded_array
        else:
            # Read Phis0 from results in Unloaded configuration
            Unloaded_file             = f"{self.pleuralpressurePath}/Unloaded_{self.region}_alpha{self.alpha_}_gravity{self.gravity_}_pe{self.pe}.xdmf"
            unloaded_reader           = xdmfReader(xdmf_file=Unloaded_file)
            self.Phis0_unloaded_array = unloaded_reader.get_celldata_from_xdmf(field_name="Phis0")

        fs_DG0                     = dolfin.FunctionSpace(self.mesh_unloaded, "DG", 0)
        self.Phis0_fun             = dolfin.Function(fs_DG0)
        self.Phis0_fun.vector()[:] = self.Phis0_unloaded_array.reshape(-1) # Working as expected, you can save it in .xdmf to check it

        # Phi_s0_fun.vector().get_local()
        # dolfin.XDMFFile("Phis0_in_unloaded_configuration.xdmf").write_checkpoint(Phis0_fun, "Phis0", 0)

        ### Define boundary of unloaded mesh
        boundaries = dolfin.MeshFunction("size_t", self.mesh_unloaded, self.mesh_unloaded.topology().dim() - 1)
        boundaries.set_all(0)
        self.ds = dolfin.Measure('ds', subdomain_data = boundaries)

        ### Assumption on fluid pressure pf
        self.assume_pf_0 = assume_pf_0 # If False, the porosity from the image data is used to compute the pleural pressure
        
                 

    def compute_stress(self, kinematics, Phis0, Phis, params):
        # material.P: 1st Piola Kirchhoff stress tensor
        # material.Sigma: 2nd Piola Kirchhoff stress tensor
        # material.sigma: Cauchy stress tensor

        ## Psi_skel: solid deformation
        solid_material = dmech.WskelLungElasticMaterial(
            kinematics=kinematics,
            parameters=params)

        material_scaling = "linear"

        material = dmech.PorousElasticMaterial(
            solid_material=solid_material,
            scaling=material_scaling,
            Phis0=Phis0)

        Sigma_skel = material.Sigma

        if Phis is not None: # we are assuming pf != 0
            ## Psi_bulk: Response of the lungs to the pressure
            solid_material = dmech.WbulkLungElasticMaterial(
                Phis=Phis,
                Phis0=Phis0,
                parameters=params)

            material = dmech.PorousElasticMaterial(
                solid_material=solid_material,
                scaling=material_scaling,
                Phis0=Phis0)
            
            # Sigma_bulk = material.dWbulkdPhis * kinematics.I # TODO Check this
            # p_f   = -1. * Sigma_bulk # Using the porosity from the images gives anomalous peak tractions in some elements

        elif Phis is None: # we are assuming pf = 0
            p_f   = 0.
        
        phis_estimated = Phis0 * kinematics.J**(-1) # TODO Check this
        
        Sigma = Sigma_skel - p_f * kinematics.J * kinematics.C_inv
        P     = kinematics.F * Sigma
        sigma = kinematics.J**(-1) * P * kinematics.F.T
    
        return P, sigma, phis_estimated



    def save_boundary_traction(self):
        phis_from_data_lst = []
        phis_estimated_lst = []
        phis_proj_lst      = []
        p_pl_lst           = []

        self.bmesh_unloaded = dolfin.BoundaryMesh(self.mesh_unloaded, "exterior")

        hdf5_file = dolfin.HDF5File(self.bmesh_unloaded.mpi_comm(), # self.mesh_unloaded
                                    f"{self.pleuralpressurePath}/data_{self.region}_alpha{self.alpha_}_gravity{self.gravity_}_pe{self.pe}.h5", "w")

        with dolfin.XDMFFile(f"{self.pleuralpressurePath}/boundary_traction_{self.region}_alpha{self.alpha_}_gravity{self.gravity_}_pe{self.pe}.xdmf") as xdmf:
            xdmf.parameters["flush_output"] = True
            xdmf.parameters["functions_share_mesh"] = True
            xdmf.parameters["rewrite_function_mesh"] = False

            for phase_i in self.MDR.U_dict_array_exhal_to_t:
                ## Load displacements
                U_exhal_to_t          = self.MDR.get_U_exhal_to_t(mesh=self.MDR.mesh_exhal, phase=phase_i)
                mesh_deformed_phase_i = self.MDR.read_mesh(self.MDR.mesh_exhal_file, U_move = [U_exhal_to_t])

                V = dolfin.VectorFunctionSpace(self.mesh_unloaded, "CG", 1)
                U_unloaded_to_t = dolfin.Function(V)
                # U_unloaded_to_inhal.vector()[:] = U_unloaded_to_exhal.vector()[:] + U_exhal_to_t.vector()[:]
                U_unloaded_to_t.vector()[:] = -self.U_exhal_to_unloaded.vector()[:] + U_exhal_to_t.vector()[:]

                kinematics_unloaded_to_t = dmech.Kinematics(U = U_unloaded_to_t)

                fs_DG0        = dolfin.FunctionSpace(self.mesh_unloaded, "DG", 0)

                if self.assume_pf_0 == False:
                    ## Load Phis
                    porosity_phase_i_file = f"{self.resultsPath}/Porosities/projected_porosity_{self.region}_{phase_i:02d}.xml"
                    phis_mf               = dolfin.MeshFunction("double", mesh_deformed_phase_i, porosity_phase_i_file)

                    phis_unloaded = dolfin.Function(fs_DG0)

                    phis_unloaded.vector()[:] = phis_mf.array() # Working as expected, you can save it in .xdmf to check it
                    phis_from_data_lst += [phis_mf.array()]
                    # phis_unloaded.vector()[:] = phis_fun.vector()[:]

                    # To check that the porosity phis (obtained from the image data and projected in the deformed mesh at each phase_i) is correctly defined in the unloaded configuration:
                    # dolfin.XDMFFile("phis_in_unloaded_configuration.xdmf").write_checkpoint(phis_unloaded, "phis", 0)

                    Phis = phis_unloaded * kinematics_unloaded_to_t.J # Currently not being used

                elif self.assume_pf_0 == True: # If we assume pf = 0, we don't need to load the porosity from the image data
                    Phis = None

                ## Compute stress
                P, sigma, phis_estimated = self.compute_stress(kinematics=kinematics_unloaded_to_t, Phis0=self.Phis0_fun, Phis=Phis, params=self.params)

                # Tested: The results are very similar to using project_on_boundary
                p_scalar, t_bound, t_bound_norm, t_tan, t_tan_norm = self.project_traction_on_boundary_3d(self.mesh_unloaded,
                                                                                                          self.bmesh_unloaded,
                                                                                                          sigma_fn=sigma,
                                                                                                          F_fn=kinematics_unloaded_to_t.F)

                ### PREVIOUS ###
                # # Working in reference configuration
                # N       = dolfin.FacetNormal(self.mesh_unloaded)
                # FinvT_N = dolfin.dot(dolfin.inv(kinematics_unloaded_to_t.F).T, N)
                # n       = FinvT_N / dolfin.sqrt(dolfin.dot(FinvT_N, FinvT_N))

                # # Working in current configuration
                # # n = dolfin.FacetNormal(mesh_inhal)

                # t = dolfin.dot(sigma, n)
                # # T = dolfin.dot(P, N)

                # # Traction on boundary
                # vfs_DG0 = dolfin.VectorFunctionSpace(self.mesh_unloaded, "DG", 0) 
                # t_bound = self.project_on_boundary(fn=t, fn_name="t", fs=vfs_DG0, ds=self.ds)

                # # Pressure (scalar) on boundary, negative is outward
                # p_scalar = self.project_on_boundary(fn=-1.*dolfin.inner(t, n), fn_name="p_pl", fs=fs_DG0, ds=self.ds)
                ### PREVIOUS ###

                U_unloaded_to_t.rename("U", "U")

                V_b = dolfin.VectorFunctionSpace(self.bmesh_unloaded, "CG", 1)
                # dolfin.parameters["allow_extrapolation"] = True
                U_b = dolfin.interpolate(U_unloaded_to_t, V_b)
                U_b.rename("U_boundary", "U_boundary")


                # phis_test_fn   = dolfin.Function(DG0)
                phis_estimated_proj = dolfin.project(phis_estimated, fs_DG0)
                phis_estimated_proj.rename("phis_estimated", "phis_estimated")
                phis_proj_lst += [phis_estimated_proj]

                phis_estimated_lst += [phis_estimated_proj.vector().get_local()]
                p_pl_lst           += [p_scalar.vector().get_local()]

                # Center of mass
                fs_DG0_phase_i = dolfin.FunctionSpace(mesh_deformed_phase_i, "DG", 0)
                phis_phase_i   = dolfin.Function(fs_DG0_phase_i)
                
                phis_phase_i.vector()[:] = phis_estimated_proj.vector().get_local()
                
                # xdmf_test = dolfin.XDMFFile("test_phis.xdmf")
                # xdmf_test.write(phis_phase_i, phase_i)
                # xdmf_test.close()
                
                dv = dolfin.Measure("dx", domain=mesh_deformed_phase_i)
                m  = dolfin.assemble(phis_phase_i * dv)
                x0 = dolfin.assemble(phis_phase_i * dolfin.SpatialCoordinate(mesh_deformed_phase_i)[0] * dv)/m
                y0 = dolfin.assemble(phis_phase_i * dolfin.SpatialCoordinate(mesh_deformed_phase_i)[1] * dv)/m
                z0 = dolfin.assemble(phis_phase_i * dolfin.SpatialCoordinate(mesh_deformed_phase_i)[2] * dv)/m
                
                fs_R = dolfin.VectorFunctionSpace(self.bmesh_unloaded, "R", 0)
                center_of_mass = dolfin.Function(fs_R) # center of mass of deformed mesh, the real function is defined in the unloaded mesh only to save it in the same file
                center_of_mass.vector()[:] = np.array([x0, y0, z0])
                center_of_mass.rename("center_of_mass", "center_of_mass")
                hdf5_file.write(center_of_mass, "/center_of_mass", phase_i)

                xdmf.write(U_b, phase_i)
                xdmf.write(t_bound, phase_i)
                # xdmf.write(t_bound_norm, phase_i)
                xdmf.write(t_tan, phase_i)
                # xdmf.write(t_tan_norm, phase_i)
                xdmf.write(p_scalar, phase_i)
                
                hdf5_file.write(p_scalar, "/p_pl", phase_i)
                hdf5_file.write(U_b, "/U_b", phase_i)
            
            hdf5_file.close()

        # with dolfin.XDMFFile(f"{self.pleuralpressurePath}/phis_{self.region}_alpha{self.alpha_}_gravity{self.gravity_}_pe{self.pe}.xdmf") as xdmf:
        #     xdmf.parameters["flush_output"] = True
        #     xdmf.parameters["functions_share_mesh"] = True
        #     xdmf.parameters["rewrite_function_mesh"] = False

        #     for phase_i in self.MDR.U_dict_array_exhal_to_t:
        #         xdmf.write(phis_proj_lst[phase_i], phase_i)

        if self.assume_pf_0 == False:
            self.phis_from_data_lst = phis_from_data_lst
        elif self.assume_pf_0 == True:
            self.phis_from_data_lst = None

        self.phis_estimated_lst = phis_estimated_lst
        self.p_pl_lst           = p_pl_lst
        return phis_from_data_lst, phis_estimated_lst, p_pl_lst



    def project_on_boundary(self, fn, fn_name, fs, ds):
        # https://fenicsproject.discourse.group/t/how-to-plot-normal-unit-vector-of-faces-in-a-2d-mesh/3912
        u = dolfin.TrialFunction(fs)
        v = dolfin.TestFunction(fs)
        a = dolfin.inner(u,v)*ds
        l = dolfin.inner(fn, v)*ds
        A = dolfin.assemble(a, keep_diagonal=True)
        L = dolfin.assemble(l)

        A.ident_zeros()
        nh = dolfin.Function(fs, name=fn_name)
        dolfin.solve(A, nh.vector(), L)
        return nh
    


    def project_traction_on_boundary_3d(self, mesh, bmesh, sigma_fn, F_fn, test_normals=False):
        """
        Computes the normal traction (t . n) on the deformed boundary of a 3D mesh.
        Projects results onto a DG0 boundary mesh to handle corners and discontinuities correctly.
        """
        
        t_fs = dolfin.TensorFunctionSpace(mesh, "DG", 0)
        
        F_func     = dolfin.project(F_fn, t_fs)
        sigma_func = dolfin.project(sigma_fn, t_fs)
        
        # Reshape to (Num_Cells, 3, 3) for direct matrix algebra
        F_array     = F_func.vector().get_local().reshape((-1, 3, 3))
        sigma_array = sigma_func.vector().get_local().reshape((-1, 3, 3))
        
        # Map: Boundary Cell Index -> Global Volume Facet Index
        map_b_to_v = bmesh.entity_map(mesh.topology().dim() - 1)
        
        # Output Function: Scalar Traction (DG0 on Boundary)
        fs_b            = dolfin.FunctionSpace(bmesh, "DG", 0)
        traction_n_func = dolfin.Function(fs_b, name="p_pl")
        tn_values       = traction_n_func.vector().get_local()

        traction_norm_func = dolfin.Function(fs_b, name="t_norm")
        t_norm_values      = traction_norm_func.vector().get_local()
        
        traction_tan_norm_func = dolfin.Function(fs_b, name="t_tan_norm")
        tt_norm_values         = traction_tan_norm_func.vector().get_local()

        # Traction vector
        vfs_b         = dolfin.VectorFunctionSpace(bmesh, "DG", 0)
        traction_func = dolfin.Function(vfs_b, name="t")
        t_values      = traction_func.vector().get_local()

        traction_tan_func = dolfin.Function(vfs_b, name="t_tan")
        tt_values         = traction_tan_func.vector().get_local()
        

        if test_normals:
            vec_fs_b = dolfin.VectorFunctionSpace(bmesh, "DG", 0)
            n_func = dolfin.Function(vec_fs_b, name="N")
            n_values_flat = n_func.vector().get_local()
            n_values_view = n_values_flat.reshape((-1, 3))


        # Initialize connectivity: Facet (dim 2) -> Cell (dim 3)
        # This is required to find the parent cell of a boundary facet
        mesh.init(2, 3)

        for cell_b in dolfin.cells(bmesh):
            b_index = cell_b.index()
            
            # Global Volume Facet
            global_facet_index = map_b_to_v[b_index]
            facet = dolfin.Facet(mesh, global_facet_index)
            
            # Get Parent Cell (The volume element owning this face)
            # entities(3) returns indices of connected tetrahedra. 
            # For a boundary facet, there is exactly one parent.
            parent_cell_index = facet.entities(3)[0]
            parent_cell = dolfin.Cell(mesh, parent_cell_index)
            
            # Get Arbitrary Normal
            # facet.normal() is robust but direction is arbitrary
            n_dolfin = facet.normal()
            N = np.array([n_dolfin.x(), n_dolfin.y(), n_dolfin.z()])
            
            # Enforce Outward Orientation (Centroid Check)
            c_facet = facet.midpoint()
            c_cell = parent_cell.midpoint()
            
            # Vector from cell center to face center (Always points OUT)
            outward_vec = np.array([c_facet.x() - c_cell.x(), 
                                    c_facet.y() - c_cell.y(), 
                                    c_facet.z() - c_cell.z()])
            
            # If normal opposes the outward vector, flip it
            if np.dot(N, outward_vec) < 0:
                N = -N
                
            if test_normals:
                n_values_view[b_index, :] = N

            
            F_mat     = F_array[parent_cell_index]
            sigma_mat = sigma_array[parent_cell_index]

            # # Nanson's Formula: n = inv(F).T * N
            # try:
            #     F_inv = np.linalg.inv(F_mat)
            # except np.linalg.LinAlgError:
            #     # Fallback for degenerate elements (avoid crash)
            #     F_inv = np.eye(3)

            # Working in reference configuration
            F_inv = np.linalg.inv(F_mat)
            FinvT = F_inv.T
            FinvT_N = np.dot(FinvT, N)
            
            # Normalize deformed normal
            norm_val   = np.linalg.norm(FinvT_N)
            n_deformed = FinvT_N / norm_val
            
            # Traction vector
            t_vec       = np.dot(sigma_mat, n_deformed)
            t_norm      = np.linalg.norm(t_vec)
            t_n_norm    = np.dot(t_vec, n_deformed)
            t_tg        = t_vec - t_n_norm * n_deformed
            t_tg_norm   = np.linalg.norm(t_tg)
            
            t_values[3*b_index + 0] = t_vec[0]
            t_values[3*b_index + 1] = t_vec[1]
            t_values[3*b_index + 2] = t_vec[2]

            t_norm_values[b_index] = t_norm

            tn_values[b_index] = -1. * t_n_norm # so negative values point outward

            tt_values[3*b_index + 0] = t_tg[0]
            tt_values[3*b_index + 1] = t_tg[1]
            tt_values[3*b_index + 2] = t_tg[2]

            tt_norm_values[b_index] = t_tg_norm

        # Save Traction
        traction_func.vector().set_local(t_values)
        traction_norm_func.vector().set_local(t_norm_values)
        traction_n_func.vector().set_local(tn_values)
        traction_tan_func.vector().set_local(tt_values)
        traction_tan_norm_func.vector().set_local(tt_norm_values)
        
        if test_normals:
            n_func.vector().set_local(n_values_view.flatten())
            xf_n = dolfin.XDMFFile("test_normals.xdmf")
            xf_n.write(n_func)
            xf_n.close()
            print(f"Normals saved to test_normals.xdmf")
        
        return traction_n_func, traction_func, traction_norm_func, traction_tan_func, traction_tan_norm_func