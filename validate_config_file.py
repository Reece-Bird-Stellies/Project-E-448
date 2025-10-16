import json, os

# ---------- Default configuration ----------
DEFAULTS = {
    "design_specs": {
        "qubit_frequency_ghz": 3.7,
        "cavity_frequency_ghz": 6.98,
        "anharmonicity_mhz": -210,
        "coupling_g_mhz": 100,
        "kappa_khz": 405,
        "resonator_type": "quarter"
    },
    "simulations": {
        "inductex": {
            "ixi": {
                "iteration_count": 500,
                "mask_segment_size": 10, 
                "mask_scaling_factor": 1.1,
                "sparam_start_freq": 2e9,
                "sparam_stop_freq": 8e9,
                "sparam_steps": 200,
                "sparam_p1_resistance": 50,
                "sparam_p2_resistance": 50
                },
            "layer_definition_file": {
                "segment_size_um": 20,
                "substrate_choice": "silicon",
                "thickness_um": 0.07,
                "lambda0_um": 0.016,
                "temperature_k": 0.01,
                "critical_temperature": 1.2,
                "h_filaments_metal": 3,
                "decimation_distance": 1
            }
        },
        "palace": {
            "capacitance": {
                "mesh_refinement": 0,
                "dielectric_material": "silicon",
                "solver_order": 1,
                "solver_tolerance": 1e-8,
                "solver_max_iterations": 500,
                "mesh_min_size": 1e-5,
                "mesh_max_size": 1.2e-4,
                "mesh_sampling": 150,
                "fillet_resolution": 12,
                "number_of_cores": 4,
                "fine_mesh_min_size_components": 1e-5,
                "fine_mesh_max_size_components": 1.2e-4
            },
            "eigenmode": {
                "mesh_refinement": 0,
                "dielectric_material": "silicon",
                "solver_order": 1,
                "solver_tolerance": 1e-5,
                "solver_max_iterations": 300,
                "mesh_min_size": 1e-5,
                "mesh_max_size": 1.2e-4,
                "mesh_sampling": 130,
                "fillet_resolution": 12,
                "number_of_cores": 4,
                "fine_mesh_min_size_components": 1e-5,
                "fine_mesh_max_size_components": 1.2e-4,
                "fine_mesh_min_size_junction": 2e-6,
                "fine_mesh_max_size_junction": 6e-5,
                "starting_frequency": 2e9,
                "number_of_frequencies": 2,
                "solutions_to_save": 0,
                "assumed_junction_capacitance": 0
            }
        }
    },
    "quantum_analysis": {
        "qubit_truncation": 30,
        "cavity_truncation": 30,
        "charge_truncation": 30
    },
    "miscellaneous": {
        "sims": {
            "base_chip_scaling_factor": 1.5,
            "extra_sims": True,
            "inductex": {
                "cap_transmon": True,
                "cap_resonator": True,
                "cap_feedline": True,
                "induc_transmon": True,
                "induc_resonator": True,
                "induc_feedline": True,
                "s_parm_transmon": True,
                "s_parm_resonator": True,
                "s_parm_whole": True,
                "efield_whole": True,
                "hfield_whole": True
            },
            "palace": {
                "cap_transmon": True,
                "cap_resonator": True,
                "cap_feedline": True,
                "eigenmode_transmon": True
            }
        }
    }
}

# ---------- Validation keys ----------
REQUIRED_PATHS = [
    "design_specs.qubit_frequency_ghz",
    "design_specs.cavity_frequency_ghz",
    "design_specs.anharmonicity_mhz",
    "design_specs.coupling_g_mhz",
    "design_specs.kappa_khz",
    "design_specs.resonator_type",
    "simulations.inductex.ixi.mask_segment_size",
    "simulations.inductex.layer_definition_file.segment_size_um",
    "simulations.palace.capacitance.mesh_min_size",
    "simulations.palace.eigenmode.starting_frequency",
    "quantum_analysis.qubit_truncation",
    "miscellaneous.sims.extra_sims",
]

# ---------- Helper functions ----------
def _get(d, path):
    """Safely retrieve nested dict values by dot path."""
    cur = d
    for p in path.split("."):
        if not isinstance(cur, dict) or p not in cur:
            return None
        cur = cur[p]
    return cur

def _is_valid(d):
    """Check that all required paths exist."""
    return all(_get(d, p) is not None for p in REQUIRED_PATHS)

def ensure_config(path):
    """
    - If file missing: create it with defaults.
    - If present but invalid: print message and return None.
    - If valid: print success message and return loaded config.
    """
    # Case 1: Missing file → create
    if not os.path.exists(path):
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        with open(path, "w") as f:
            json.dump(DEFAULTS, f, indent=2)
        print(f"[config] Created new default config at '{path}'")
        return DEFAULTS

    # Case 2: File exists → try reading
    try:
        with open(path, "r") as f:
            data = json.load(f)
    except Exception as e:
        print(f"[config] Could not read '{path}': {e}")
        print("[config] Leaving file unchanged.")
        return None

    # Case 3: Validate structure
    if _is_valid(data):
        print(f"[config] Validation successful: '{path}' is valid.")
        return data

    # Case 4: Invalid → report missing keys
    print(f"[config] Invalid structure detected in '{path}'. Missing keys:")
    for p in REQUIRED_PATHS:
        if _get(data, p) is None:
            print(f"  - {p}")
    print("[config] Leaving file unchanged.")
    return None
