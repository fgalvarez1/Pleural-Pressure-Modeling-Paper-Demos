import meshio
import xml.etree.ElementTree as ET

class xdmfReader:
    def __init__(self, xdmf_file):
        """
        Read an XDMF file (saved as a TimeSeries), extract and save the last time step, and read a field from it.
        It can be adapted to extract any field from any time step.
        """
        self.xdmf_file       = xdmf_file # File with results from simulation (TimeSeries)
        self.output_filename = xdmf_file.replace(".xdmf", "_final.xdmf")
        
        self.write_last_time_of_xdmf()



    def write_last_time_of_xdmf(self):
        # output_prefix = "Unloaded0.160"
        # output_folder = "split_xdmfs"
        # os.makedirs(output_folder, exist_ok=True)

        tree = ET.parse(self.xdmf_file)
        root = tree.getroot()

        # Get namespace
        ns = {'xi': 'http://www.w3.org/2001/XInclude'}

        # Find temporal Grid
        temporal_grid = root.find(".//Grid[@CollectionType='Temporal']")
        
        # Select last time step
        i    = len(temporal_grid.findall("Grid"))-1
        grid = temporal_grid.findall("Grid")[i]

        # Read all time steps
        # for i, grid in enumerate(temporal_grid.findall("Grid")):

        time_val = grid.find("Time").attrib["Value"]

        # Remove the <Time> tag
        time_tag = grid.find("Time")
        if time_tag is not None:
            grid.remove(time_tag)

        # Create new XDMF root
        new_root = ET.Element("Xdmf", Version="2.0")
        domain = ET.SubElement(new_root, "Domain")
        domain.append(grid)

        # Save the modified file
        new_tree = ET.ElementTree(new_root)
        # output_file = os.path.join(output_folder, f"{output_prefix}_t{time_val}.xdmf")
        new_tree.write(self.output_filename, xml_declaration=True, encoding="utf-8")



    def get_celldata_from_xdmf(self, field_name):
        # Read fields
        read_file    = meshio.read(self.output_filename) # Eventually, this could read any specific time step
        field_values = read_file.cell_data[field_name][0]
        # u_last   = read_file.point_data["u"]

        # mesh = dolfin.Mesh()
        # with dolfin.XDMFFile(xdmf_file) as xdmf:
        #     xdmf.read(mesh)

        # DG0   = dolfin.FunctionSpace(mesh, "DG", 0)
        # Phis0 = dolfin.Function(DG0)
        # Phis0.vector()[:] = Phis0_last.reshape(-1)

        # V = dolfin.VectorFunctionSpace(mesh, "CG", 1)
        # u = dolfin.Function(V)
        # u.vector()[:] = u_last.reshape(-1)

        return field_values