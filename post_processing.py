import numpy as np
from scipy.constants import h, e
e, h
def post_processing(ref_design_specs, ref_cap_data, results_inductex, results_palace):
    mode_alt_cavity = "mode_2"
    Fc_hfss         = ref_design_specs["cavity_frequency_ghz"] * 1e9  # Convert GHz to Hz 
    Lr_inductex     = results_inductex["inductance"]["resonator"]["L1"] # Inductance not given so use inductex result
    Ec              = ref_cap_data["EC"] * 1e9 * h  # Convert GHz to Hz and then to J
    Cs              = ref_cap_data["cross_to_cross"] * 1e-15  # Convert fF to F
    Cg              = ref_cap_data["cross_to_claw"] * 1e-15  # Convert fF to F
    Fc_palace       = results_palace["eigenmode"]["full"]["eigen_mode_frequencies"][mode_alt_cavity] * 1e9                         # Changed this to use frequency from full
    Cr_palace       = results_palace["capacitance"]["resonator"]["C2-C2"]  # PALACE cannot do inductance sims so use eigenfrequencies to get Lr

    cj_hfss         = e**2 / (2*Ec) - Cg - Cs
    cr_hfss         = 1 / ((2*np.pi*Fc_hfss)**2 * Lr_inductex)
    lr_palace       = 1 / ((2*np.pi*Fc_palace)**2 * Cr_palace)

    return {
        "cj_hfss": cj_hfss,
        "cr_hfss": cr_hfss,
        "lr_palace": lr_palace
    }