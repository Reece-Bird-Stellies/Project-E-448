import json, os
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
                "sparams_start_freq": 2e9,
                "sparams_stop_freq": 8e9,
                "sparams_steps": 200,
                "sparams_p1_resistance": 50,
                "sparams_p2_resistance": 50
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
                "solver_order": 2,
                "solver_tolerance": 1e-8,
                "solver_max_iterations": 500,
                "mesh_min_size": 1e-5,
                "mesh_max_size": 1.2e-4,
                "mesh_sampling": 150,
                "fillet_resolution": 12,
                "number_of_cores": 4,
                "fine_mesh_min_size_components": 1e-5,
                "fine_mesh_max_size_components": 1.2e-4,
                "solutions_to_save": 0,
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
                "fine_mesh_min_size_junction": 6e-6,
                "fine_mesh_max_size_junction": 6e-5,
                "starting_frequency": 2e9,
                "number_of_frequencies": 2,
                "solutions_to_save": 0,
                "assumed_junction_capacitance": 0,
                "port_1_resistance": 50,
                "port_1_inductance": 0,
                "port_1_capacitance": 0,
                "port_2_resistance": 1e-12,
                "port_2_inductance": 0,
                "port_2_capacitance": 0
             }
        }
    },
    "quantum_analysis": {
        "lom": {
            "qubit_truncation": 20,
            "cavity_truncation": 20,
            "charge_truncation": 20
        },
        "epr": {
            "mode_truncation": 20
        }
    },
    "miscellaneous": {
        "sims": {
            "base_chip_scaling_factor": 1.5,
            "extra_sims": True,
            "inductex": {
                "cap_transmon": True,
                "cap_transmon_no_jj": True,
                "cap_resonator": True,
                "cap_feedline": True,
                "induc_resonator": True,
                "s_parm_transmon": True,
                "s_parm_transmon_no_jj": True,
                "s_parm_resonator": True,
                "efield_whole": True,
                "hfield_whole": True
            },
            "palace": {
                "cap_transmon": True,
                "cap_transmon_no_jj": True,
                "cap_resonator": True,
                "cap_feedline": True,
                "eigenmode_full": True,
                "eigenmode_full_no_jj": True,
                "eigenmode_resonator": True
            }
        }
    },
    "lom": {
        "mesh_size": 1e-6,
        "solver_tolerance": 1e-9,
        "max_iterations": 1000,
        "material": "niobium",
        "temperature_k": 4.2
    },
    "epr": {
        "solver_order": 2,
        "solver_tolerance": 1e-6,
        "max_iterations": 500,
        "mesh_min_size": 1e-5,
        "mesh_max_size": 1e-4,
        "number_of_modes": 5,
        "save_fields": True
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
    "simulations.inductex.ixi.iteration_count",
    "simulations.inductex.ixi.mask_segment_size",
    "simulations.inductex.layer_definition_file.segment_size_um",
    "simulations.palace.capacitance.mesh_min_size",
    "simulations.palace.eigenmode.starting_frequency",
    "quantum_analysis.lom.qubit_truncation",
    "quantum_analysis.lom.cavity_truncation",
    "quantum_analysis.lom.charge_truncation",
    "quantum_analysis.epr.mode_truncation",
    "miscellaneous.sims.extra_sims",
    "miscellaneous.sims.inductex.cap_transmon",
    "miscellaneous.sims.palace.eigenmode_full",
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

def _set(d, path, value):
    """Safely set nested dict values by dot path."""
    keys = path.split(".")
    cur = d
    for key in keys[:-1]:
        if key not in cur:
            cur[key] = {}
        cur = cur[key]
    cur[keys[-1]] = value

def _compare_dicts(loaded, defaults, path=""):
    """
    Recursively compare loaded config with defaults.
    Print differences found.
    """
    if isinstance(defaults, dict):
        for key, default_value in defaults.items():
            new_path = f"{path}.{key}" if path else key
            
            if key not in loaded:
                continue
            
            loaded_value = loaded[key]
            
            if isinstance(default_value, dict):
                _compare_dicts(loaded_value, default_value, new_path)
            else:
                if loaded_value != default_value:
                    # Format the path for display
                    display_path = new_path.replace(".", '"]["')
                    print(f'[INFO] {display_path} is: {loaded_value} (Default: {default_value})')

def _is_valid(d):
    """Check that all required paths exist."""
    return all(_get(d, p) is not None for p in REQUIRED_PATHS)

def ensure_config(path):
    """
    - If file missing: create it with defaults.
    - If present but invalid: print message and return None.
    - If valid: print success message, compare with defaults, and return loaded config.
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
        
        # Compare with defaults and report differences
        _compare_dicts(data, DEFAULTS)
        
        return data

    # Case 4: Invalid → report missing keys
    print(f"[config] Invalid structure detected in '{path}'. Missing keys:")
    for p in REQUIRED_PATHS:
        if _get(data, p) is None:
            print(f"  - {p}")
    print("[config] Leaving file unchanged.")
    return None