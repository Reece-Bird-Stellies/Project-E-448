from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation
from SQDMetal.PALACE.Capacitance_Simulation import PALACE_Capacitance_Simulation
from SQDMetal.Utilities.Materials import MaterialInterface
from SQDMetal.Utilities.QiskitShapelyRenderer import QiskitShapelyRenderer 
import numpy as np

def make_palace_cap_sim(design, design_name, palace_cap_config):

    path_to_output = "simulations/palace/capacitance/"

    # Define simulation options as variables for clarity and alignment
    mesh_refinement                 = palace_cap_config["mesh_refinement"]
    dielectric_material             = palace_cap_config["dielectric_material"]
    solver_order                    = palace_cap_config["solver_order"]
    solver_tolerance                = palace_cap_config["solver_tolerance"]
    solver_max_iterations           = palace_cap_config["solver_max_iterations"]
    mesh_min_size                   = palace_cap_config["mesh_min_size"]
    mesh_max_size                   = palace_cap_config["mesh_max_size"]
    mesh_sampling                   = palace_cap_config["mesh_sampling"]
    fillet_resolution               = palace_cap_config["fillet_resolution"]
    number_of_cores                 = palace_cap_config["number_of_cores"]              # Doesnt do anything here
    fine_mesh_min_size_components   = palace_cap_config["fine_mesh_min_size_components"]
    fine_mesh_max_size_components   = palace_cap_config["fine_mesh_max_size_components"]
    solutions_to_save               = palace_cap_config["solutions_to_save"]
    palace_dir                      = ""

    user_defined_options            = {
        "mesh_refinement":     mesh_refinement,
        "dielectric_material": dielectric_material,
        "solver_order":        solver_order,
        "solver_tol":          solver_tolerance,
        "solver_maxits":       solver_max_iterations,
        "mesh_min":            mesh_min_size,
        "mesh_max":            mesh_max_size,
        "mesh_sampling":       mesh_sampling,
        "fillet_resolution":   fillet_resolution,
        "num_cpus":            number_of_cores,
        "solns_to_save":       solutions_to_save,
        "palace_dir":          palace_dir
    }

    #Create the Palace Eigenmode simulation
    cap_sim = PALACE_Capacitance_Simulation(name                    = design_name,             #name of simulation
                                            metal_design            = design,                  #feed in qiskit metal design
                                            sim_parent_directory    = path_to_output,          #choose directory where mesh file, config file and HPC batch file will be saved
                                            mode                    = 'simPC',                 #choose simulation mode 'HPC' or 'simPC'                                          
                                            meshing                 = 'GMSH',                  #choose meshing 'GMSH' or 'COMSOL'
                                            user_options            = user_defined_options,    #provide options chosen above
                                            view_design_gmsh_gui    = False,                   #view design in GMSH gui 
                                            create_files            = True,
                                            ) 

    # Get all geometries from the design  
    qmpl         = QiskitShapelyRenderer(None, design, None)  
    gsdf        = qmpl.get_net_coordinates(resolution=4)  
    
    # Extract unique layer IDs  
    layer_ids   = np.unique(gsdf['layer'])  
    print(f"Layers in design: {layer_ids}")

    cap_sim.add_ground_plane(threshold=1e-10)

    if design_name == "transmon":
        fine_mesh_components = ['transmon', 'junction']
    if design_name == "transmon_no_jj":
        fine_mesh_components = ['transmon']
    else:
        fine_mesh_components = [design_name]

    for layer in layer_ids:
        cap_sim.add_metallic(layer, threshold=1e-10, fuse_threshold=1e-10)

    cap_sim.fine_mesh_around_comp_boundaries(fine_mesh_components, min_size=fine_mesh_min_size_components, max_size=fine_mesh_max_size_components) 

    cap_sim.prepare_simulation()
    cap_sim.display_conductor_indices()

def _setup_full_design(eigen_sim, design, ref_design_junction_inductance, assumed_junction_capacitance, mesh_sampling, fine_mesh_min_size_junction, fine_mesh_max_size_junction):
    fine_mesh_components = ['resonator', 'transmon', 'feedline', 'LP1', 'LP2']
    eigen_sim.create_port_CPW_on_Launcher('LP1', 20e-6)
    eigen_sim.create_port_CPW_on_Launcher('LP2', 20e-6)
    eigen_sim.create_port_JosephsonJunction('junction', L_J=ref_design_junction_inductance, C_J=assumed_junction_capacitance)
    bounds = design.components["junction"].qgeometry_bounds()
    bounds = bounds * 1e-3
    eigen_sim.fine_mesh_in_rectangle(bounds[0], bounds[1], bounds[2], bounds[3],
                                    mesh_sampling  = mesh_sampling, 
                                    min_size       = fine_mesh_min_size_junction, 
                                    max_size       = fine_mesh_max_size_junction
                                    )
    eigen_sim.set_port_impedance(port_ind=1, impedance_R=45.6, impedance_L=0, impedance_C=0)
    eigen_sim.set_port_impedance(port_ind=2, impedance_R=0, impedance_L=1e-15, impedance_C=0)
    return fine_mesh_components

def _setup_transmon_design(eigen_sim, design, ref_design_junction_inductance, assumed_junction_capacitance, mesh_sampling, fine_mesh_min_size_junction, fine_mesh_max_size_junction):
    fine_mesh_components = ['transmon']
    eigen_sim.create_port_CPW_on_Route('transmon', pin_name='readout', len_launch=20e-6)
    eigen_sim.create_port_JosephsonJunction('junction', L_J=ref_design_junction_inductance, C_J=assumed_junction_capacitance)
    bounds = design.components["junction"].qgeometry_bounds()
    bounds = bounds * 1e-3
    eigen_sim.fine_mesh_in_rectangle(bounds[0], bounds[1], bounds[2], bounds[3],
                                    mesh_sampling  = mesh_sampling, 
                                    min_size       = fine_mesh_min_size_junction, 
                                    max_size       = fine_mesh_max_size_junction
                                    )
    return fine_mesh_components

def _setup_resonator_design(eigen_sim):
    fine_mesh_components = ['resonator']
    eigen_sim.create_port_CPW_on_Route('resonator', pin_name='start', len_launch=20e-6)
    eigen_sim.create_port_CPW_on_Route('resonator', pin_name='end', len_launch=20e-6)
    eigen_sim.set_port_impedance(port_ind=1, impedance_R=50, impedance_L=0, impedance_C=0)
    eigen_sim.set_port_impedance(port_ind=2, impedance_R=0, impedance_L=1e-15, impedance_C=0)
    return fine_mesh_components

def make_palace_eigenmode_sim(design, design_name, palace_eigenmode_config):
    path_to_output = "simulations/palace/eigenmode/"

    # Extract Eigenmode Simulation Options from config
    mesh_refinement                 = palace_eigenmode_config["mesh_refinement"]
    dielectric_material             = palace_eigenmode_config["dielectric_material"]
    solver_order                    = palace_eigenmode_config["solver_order"]
    solver_tolerance                = palace_eigenmode_config["solver_tolerance"]
    solver_max_iterations           = palace_eigenmode_config["solver_max_iterations"]
    mesh_min_size                   = palace_eigenmode_config["mesh_min_size"]
    mesh_max_size                   = palace_eigenmode_config["mesh_max_size"]
    mesh_sampling                   = palace_eigenmode_config["mesh_sampling"]
    fillet_resolution               = palace_eigenmode_config["fillet_resolution"]
    number_of_cores                 = palace_eigenmode_config["number_of_cores"]
    fine_mesh_min_size_components   = palace_eigenmode_config["fine_mesh_min_size_components"]
    fine_mesh_max_size_components   = palace_eigenmode_config["fine_mesh_max_size_components"]
    fine_mesh_min_size_junction     = palace_eigenmode_config["fine_mesh_min_size_junction"]
    fine_mesh_max_size_junction     = palace_eigenmode_config["fine_mesh_max_size_junction"]
    starting_frequency              = palace_eigenmode_config["starting_frequency"]
    number_of_frequencies           = palace_eigenmode_config["number_of_frequencies"]
    solutions_to_save               = palace_eigenmode_config["solutions_to_save"]
    assumed_junction_capacitance    = palace_eigenmode_config["assumed_junction_capacitance"]
    ref_design_junction_inductance  = palace_eigenmode_config["lj"]
    path_to_palace                  = ""
    user_defined_options = {
        "mesh_refinement":     mesh_refinement,
        "dielectric_material": dielectric_material,
        "solver_order":        solver_order,
        "solver_tol":          solver_tolerance,
        "solver_maxits":       solver_max_iterations,
        "mesh_min":            mesh_min_size,
        "mesh_max":            mesh_max_size,
        "mesh_sampling":       mesh_sampling,
        "fillet_resolution":   fillet_resolution,
        "num_cpus":            number_of_cores,    # Doesnt do anything here
        "starting_freq":       starting_frequency,
        "number_of_freqs":     number_of_frequencies,
        "solns_to_save":       solutions_to_save,
        "palace_dir":          path_to_palace
    }

    # Create the Palace Eigenmode simulation
    eigen_sim = PALACE_Eigenmode_Simulation(name                    = design_name,
                                            metal_design            = design,
                                            sim_parent_directory    = path_to_output,
                                            mode                    = 'simPC',
                                            meshing                 = 'GMSH',
                                            user_options            = user_defined_options,
                                            view_design_gmsh_gui    = False,
                                            create_files            = True
                                            )

    # Get all geometries from the design  
    qmpl         = QiskitShapelyRenderer(None, design, None)  
    gsdf        = qmpl.get_net_coordinates(resolution=4)  
    
    # Extract unique layer IDs  
    layer_ids   = np.unique(gsdf['layer'])  
    print(f"Layers in design: {layer_ids}")
    
    for layer in layer_ids:
        eigen_sim.add_metallic(layer, threshold=1e-10, fuse_threshold=1e-10)
    eigen_sim.add_ground_plane(threshold=1e-10)
    eigen_sim.setup_EPR_interfaces(metal_air=MaterialInterface('Aluminium-Vacuum'),
                                substrate_air=MaterialInterface('Silicon-Vacuum'),
                                substrate_metal=MaterialInterface('Silicon-Aluminium'))

    # Call helper functions based on design_name
    if design_name == "full":
        fine_mesh_components = _setup_full_design(eigen_sim, design, ref_design_junction_inductance, assumed_junction_capacitance, mesh_sampling, fine_mesh_min_size_junction, fine_mesh_max_size_junction)
    elif design_name == "transmon":
        fine_mesh_components = _setup_transmon_design(eigen_sim, design, ref_design_junction_inductance, assumed_junction_capacitance, mesh_sampling, fine_mesh_min_size_junction, fine_mesh_max_size_junction)
    elif design_name == "resonator":
        fine_mesh_components = _setup_resonator_design(eigen_sim)

    eigen_sim.fine_mesh_around_comp_boundaries(fine_mesh_components,
                                               min_size=fine_mesh_min_size_components, 
                                               max_size=fine_mesh_max_size_components
                                               )

    eigen_sim.prepare_simulation()


