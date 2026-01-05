import numpy as np
from scipy.constants import h, e 

def post_processing(ref_design_specs, ref_cap_data, results_inductex, results_palace):
    mode_alt_cavity = "mode_2"
    try: 
        Fc_hfss         = ref_design_specs["cavity_frequency_ghz"] * 1e9            # Convert GHz to Hz 
        Lr_inductex     = results_inductex["inductance"]["resonator"]["L1"]         # Inductance not given so use inductex result
        Ec              = ref_cap_data["EC"] * 1e9 * h                              # Convert GHz to Hz and then to J
        Cs              = ref_cap_data["cross_to_cross"] * 1e-15                    # Convert fF to F
        Cg              = ref_cap_data["cross_to_claw"] * 1e-15                     # Convert fF to F
        cj_hfss         = e**2 / (2*Ec) - Cg - Cs
        cr_hfss         = 1 / ((2*np.pi*Fc_hfss)**2 * Lr_inductex)
    except TypeError:
        print("InductEx results or reference data are incomplete. Ensure that inductance simulations have been run and reference design data is available.")
        Fc_hfss         = None
        Lr_inductex     = None
        Ec              = None
        Cs              = None
        Cg              = None
        cj_hfss         = None
        cr_hfss         = None

    try:
        Fc_palace       = results_palace["eigenmode"]["full"]["eigen_mode_frequencies"][mode_alt_cavity] * 1e9                         # Changed this to use frequency from full
        Cr_palace       = results_palace["capacitance"]["resonator"]["C2-C2"]                                                          # PALACE cannot do inductance sims so use eigenfrequencies to get Lr
        lr_palace       = 1 / ((2*np.pi*Fc_palace)**2 * Cr_palace)
    except TypeError:
        print("PALACE results are incomplete. Ensure that eigenmode simulations have been run and capacitance data is available.")
        Fc_palace       = None
        Cr_palace       = None
        lr_palace       = None
    
    return {
        "cj_hfss": cj_hfss,
        "cr_hfss": cr_hfss,
        "lr_palace": lr_palace
    }