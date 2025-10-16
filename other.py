


qubit_frequency_ghz                          = design_specs["qubit_frequency_ghz"]
cavity_frequency_ghz                         = design_specs["cavity_frequency_ghz"]
anharmonicity_mhz                            = design_specs["anharmonicity_mhz"]
coupling_g_mhz                               = design_specs["coupling_g_mhz"]
kappa_khz                                    = design_specs["kappa_khz"]
resonator_type                               = design_specs["resonator_type"]


# ========= SIMULATIONS =========
simulations                                  = data["simulations"]

# ----- INDUCTEX -----
inductex                                     = simulations["inductex"]

ixi                                          = inductex["ixi"]
mask_segment_size                            = ixi["mask_segment_size"]
mask_scaling_factor                          = ixi["mask_scaling_factor"]

layer_definition_file                        = inductex["layer_definition_file"]
segment_size_um                              = layer_definition_file["segment_size_um"]
substrate_choice                             = layer_definition_file["substrate_choice"]
thickness_um                                 = layer_definition_file["thickness_um"]
temperature_k                                = layer_definition_file["temperature_k"]
critical_temperature                         = layer_definition_file["critical_temperature"]
h_filaments_metal                            = layer_definition_file["h_filaments_metal"]
metal_choice                                 = layer_definition_file["metal_choice"]
decimation_distance                          = layer_definition_file["decimation_distance"]

# ----- PALACE -----
palace                                       = simulations["palace"]

capacitance                                  = palace["capacitance"]
mesh_refinement_cap                          = capacitance["mesh_refinement"]
dielectric_material_cap                      = capacitance["dielectric_material"]
solver_order_cap                             = capacitance["solver_order"]
solver_tolerance_cap                         = capacitance["solver_tolerance"]
solver_max_iterations_cap                    = capacitance["solver_max_iterations"]
mesh_min_size_cap                            = capacitance["mesh_min_size"]
mesh_max_size_cap                            = capacitance["mesh_max_size"]
mesh_sampling_cap                            = capacitance["mesh_sampling"]
fillet_resolution_cap                        = capacitance["fillet_resolution"]
number_of_cores_cap                          = capacitance["number_of_cores"]
fine_mesh_min_size_components_cap            = capacitance["fine_mesh_min_size_components"]
fine_mesh_max_size_components_cap            = capacitance["fine_mesh_max_size_components"]

eigenmode                                    = palace["eigenmode"]
mesh_refinement_eig                          = eigenmode["mesh_refinement"]
dielectric_material_eig                      = eigenmode["dielectric_material"]
solver_order_eig                             = eigenmode["solver_order"]
solver_tolerance_eig                         = eigenmode["solver_tolerance"]
solver_max_iterations_eig                    = eigenmode["solver_max_iterations"]
mesh_min_size_eig                            = eigenmode["mesh_min_size"]
mesh_max_size_eig                            = eigenmode["mesh_max_size"]
mesh_sampling_eig                            = eigenmode["mesh_sampling"]
fillet_resolution_eig                        = eigenmode["fillet_resolution"]
number_of_cores_eig                          = eigenmode["number_of_cores"]
fine_mesh_min_size_components_eig            = eigenmode["fine_mesh_min_size_components"]
fine_mesh_max_size_components_eig            = eigenmode["fine_mesh_max_size_components"]
fine_mesh_min_size_junction                  = eigenmode["fine_mesh_min_size_junction"]
fine_mesh_max_size_junction                  = eigenmode["fine_mesh_max_size_junction"]
starting_frequency_ghz                       = eigenmode["starting_frequency_ghz"]
number_of_frequencies                        = eigenmode["number_of_frequencies"]
solutions_to_save                            = eigenmode["solutions_to_save"]
assumed_junction_capacitance                 = eigenmode["assumed_junction_capacitance"]


# ========= QUANTUM ANALYSIS =========
quantum_analysis                             = data["quantum_analysis"]

qubit_truncation                             = quantum_analysis["qubit_truncation"]
cavity_truncation                            = quantum_analysis["cavity_truncation"]
charge_truncation                            = quantum_analysis["charge_truncation"]


# ========= MISCELLANEOUS =========
miscellaneous                                = data["miscellaneous"]

sims                                         = miscellaneous["sims"]
extra_sims                                   = sims["extra_sims"]

# ----- INDUCTEX -----
inductex_misc                                = sims["inductex"]
cap_transmon                                 = inductex_misc["cap_transmon"]
cap_resonator                                = inductex_misc["cap_resonator"]
cap_feedline                                 = inductex_misc["cap_feedline"]
induc_transmon                               = inductex_misc["induc_transmon"]
induc_resonator                              = inductex_misc["induc_resonator"]
induc_feedline                               = inductex_misc["induc_feedline"]
s_parm_transmon                              = inductex_misc["s_parm_transmon"]
s_parm_resonator                             = inductex_misc["s_parm_resonator"]
s_parm_whole                                 = inductex_misc["s_parm_whole"]
efield_whole                                 = inductex_misc["efield_whole"]
hfield_whole                                 = inductex_misc["hfield_whole"]

# ----- PALACE -----
palace_misc                                  = sims["palace"]
cap_transmon_palace                          = palace_misc["cap_transmon"]
cap_resonator_palace                         = palace_misc["cap_resonator"]
cap_feedline_palace                          = palace_misc["cap_feedline"]
eigenmode_transmon                           = palace_misc["eigenmode_transmon"]












# ========= PRINT VARIABLE TYPES =========

print("========== DESIGN SPECS =========")
print("qubit_frequency_ghz                :", type(qubit_frequency_ghz))
print("cavity_frequency_ghz               :", type(cavity_frequency_ghz))
print("anharmonicity_mhz                  :", type(anharmonicity_mhz))
print("coupling_g_mhz                     :", type(coupling_g_mhz))
print("kappa_khz                          :", type(kappa_khz))
print("resonator_type                     :", type(resonator_type))


print("\n========== SIMULATIONS / INDUCTEX / IXI =========")
print("mask_segment_size                  :", type(mask_segment_size))
print("mask_scaling_factor                :", type(mask_scaling_factor))

print("\n========== SIMULATIONS / INDUCTEX / LAYER DEFINITION FILE =========")
print("segment_size_um                    :", type(segment_size_um))
print("substrate_choice                   :", type(substrate_choice))
print("thickness_um                       :", type(thickness_um))
print("temperature_k                      :", type(temperature_k))
print("critical_temperature               :", type(critical_temperature))
print("h_filaments_metal                  :", type(h_filaments_metal))
print("metal_choice                       :", type(metal_choice))
print("decimation_distance                :", type(decimation_distance))


print("\n========== SIMULATIONS / PALACE / CAPACITANCE =========")
print("mesh_refinement_cap                :", type(mesh_refinement_cap))
print("dielectric_material_cap            :", type(dielectric_material_cap))
print("solver_order_cap                   :", type(solver_order_cap))
print("solver_tolerance_cap               :", type(solver_tolerance_cap))
print("solver_max_iterations_cap          :", type(solver_max_iterations_cap))
print("mesh_min_size_cap                  :", type(mesh_min_size_cap))
print("mesh_max_size_cap                  :", type(mesh_max_size_cap))
print("mesh_sampling_cap                  :", type(mesh_sampling_cap))
print("fillet_resolution_cap              :", type(fillet_resolution_cap))
print("number_of_cores_cap                :", type(number_of_cores_cap))
print("fine_mesh_min_size_components_cap  :", type(fine_mesh_min_size_components_cap))
print("fine_mesh_max_size_components_cap  :", type(fine_mesh_max_size_components_cap))


print("\n========== SIMULATIONS / PALACE / EIGENMODE =========")
print("mesh_refinement_eig                :", type(mesh_refinement_eig))
print("dielectric_material_eig            :", type(dielectric_material_eig))
print("solver_order_eig                   :", type(solver_order_eig))
print("solver_tolerance_eig               :", type(solver_tolerance_eig))
print("solver_max_iterations_eig          :", type(solver_max_iterations_eig))
print("mesh_min_size_eig                  :", type(mesh_min_size_eig))
print("mesh_max_size_eig                  :", type(mesh_max_size_eig))
print("mesh_sampling_eig                  :", type(mesh_sampling_eig))
print("fillet_resolution_eig              :", type(fillet_resolution_eig))
print("number_of_cores_eig                :", type(number_of_cores_eig))
print("fine_mesh_min_size_components_eig  :", type(fine_mesh_min_size_components_eig))
print("fine_mesh_max_size_components_eig  :", type(fine_mesh_max_size_components_eig))
print("fine_mesh_min_size_junction        :", type(fine_mesh_min_size_junction))
print("fine_mesh_max_size_junction        :", type(fine_mesh_max_size_junction))
print("starting_frequency_ghz             :", type(starting_frequency_ghz))
print("number_of_frequencies              :", type(number_of_frequencies))
print("solutions_to_save                  :", type(solutions_to_save))
print("assumed_junction_capacitance       :", type(assumed_junction_capacitance))


print("\n========== QUANTUM ANALYSIS =========")
print("qubit_truncation                   :", type(qubit_truncation))
print("cavity_truncation                  :", type(cavity_truncation))
print("charge_truncation                  :", type(charge_truncation))


print("\n========== MISCELLANEOUS / SIMS =========")
print("extra_sims                         :", type(extra_sims))

print("\n----- INDUCTEX -----")
print("cap_transmon                       :", type(cap_transmon))
print("cap_resonator                      :", type(cap_resonator))
print("cap_feedline                       :", type(cap_feedline))
print("induc_transmon                     :", type(induc_transmon))
print("induc_resonator                    :", type(induc_resonator))
print("induc_feedline                     :", type(induc_feedline))
print("s_parm_transmon                    :", type(s_parm_transmon))
print("s_parm_resonator                   :", type(s_parm_resonator))
print("s_parm_whole                       :", type(s_parm_whole))
print("efield_whole                       :", type(efield_whole))
print("hfield_whole                       :", type(hfield_whole))

print("\n----- PALACE -----")
print("cap_transmon_palace                :", type(cap_transmon_palace))
print("cap_resonator_palace               :", type(cap_resonator_palace))
print("cap_feedline_palace                :", type(cap_feedline_palace))
print("eigenmode_transmon                 :", type(eigenmode_transmon))
