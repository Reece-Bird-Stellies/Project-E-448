from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation
from SQDMetal.PALACE.Capacitance_Simulation import PALACE_Capacitance_Simulation
from SQDMetal.Utilities.Materials import MaterialInterface

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
    number_of_cores                 = palace_cap_config["number_of_cores"]
    fine_mesh_min_size_components   = palace_cap_config["fine_mesh_min_size_components"]
    fine_mesh_min_size_components   = palace_cap_config["fine_mesh_min_size_components"]
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
        "palace_dir":          palace_dir
    }

    #Creat the Palace Eigenmode simulation
    cap_sim = PALACE_Capacitance_Simulation(name                    = design_name,             #name of simulation
                                            metal_design            = design,                  #feed in qiskit metal design
                                            sim_parent_directory    = path_to_output,          #choose directory where mesh file, config file and HPC batch file will be saved
                                            mode                    = 'simPC',                 #choose simulation mode 'HPC' or 'simPC'                                          
                                            meshing                 = 'GMSH',                  #choose meshing 'GMSH' or 'COMSOL'
                                            user_options            = user_defined_options,    #provide options chosen above
                                            view_design_gmsh_gui    = False,                   #view design in GMSH gui 
                                            create_files            = True,
                                            ) 

    cap_sim.add_metallic(1, threshold=1e-10, fuse_threshold=1e-10)
    cap_sim.add_ground_plane(threshold=1e-10)

    cap_sim.fine_mesh_around_comp_boundaries([design_name], min_size=fine_mesh_min_size_components, max_size=fine_mesh_min_size_components) 

    cap_sim.prepare_simulation()
    cap_sim.display_conductor_indices()

def make_palace_cap_sim_transmon_no_jj_gds(transmon_no_jj_design, palace_cap_config):
    return