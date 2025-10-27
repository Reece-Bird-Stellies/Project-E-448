from squadds import SQuADDS_DB, Analyzer;
import os
import json
import pandas as pd
import numpy as np
from scipy.constants import e, hbar, h

def _load_best_design_from_json(file_path="data/design_reference_data_all.json"):
    """
    Load the saved best_design from a JSON file and return it.
    If pandas is available, return a pandas.Series, otherwise return the raw dict/list.
    """
    with open(file_path, "r") as f:
        content = json.load(f)

    # If the JSON wraps the object under a key, try common names
    if isinstance(content, dict):
        if "best_design" in content:
            obj = content["best_design"]
        elif "data" in content and "best_design" in content["data"]:
            obj = content["data"]["best_design"]
        else:
            obj = content
    else:
        obj = content

    try:
        # If it's a mapping, convert to Series
        if isinstance(obj, dict):
            return pd.Series(obj)
        # If it's a list of [key, value] pairs
        if isinstance(obj, list) and all(isinstance(e, (list, tuple)) and len(e) == 2 for e in obj):
            return pd.Series({k: v for k, v in obj})
    except Exception:
        # pandas not available or conversion failed; fall back to raw object
        pass

    return obj

def fetch_design(design_specs, num_top=1):
    """
    Uses only design_specs to query SQuADDS, extracts the top design,
    and returns: best_design, ref_design_specs, ref_capacitance_specs,
                 qubit_geo, resonator_geo, feedline_geo, LJs
    """

    # ========= EXPLICIT VARIABLE DECLARATION =========
    qubit_frequency_ghz                          = design_specs["qubit_frequency_ghz"]
    cavity_frequency_ghz                         = design_specs["cavity_frequency_ghz"]
    anharmonicity_mhz                            = design_specs["anharmonicity_mhz"]
    coupling_g_mhz                               = design_specs["coupling_g_mhz"]
    kappa_khz                                    = design_specs["kappa_khz"]
    resonator_type                               = design_specs["resonator_type"]

    # ========= INITIALIZE DATABASE AND ANALYZER =========
    db                                           = SQuADDS_DB()
    db.unselect_all()
    db.select_system(["cavity_claw", "qubit"])
    db.select_qubit("TransmonCross")
    db.select_cavity_claw("RouteMeander")
    db.select_resonator_type(resonator_type)
    df                                           = db.create_system_df()
    analyzer                                     = Analyzer(db)

    # ========= BUILD TARGET PARAMETERS =========
    target_params                                = {
        "qubit_frequency_GHz":  qubit_frequency_ghz,
        "cavity_frequency_GHz": cavity_frequency_ghz,
        "anharmonicity_MHz":    anharmonicity_mhz,
        "kappa_kHz":            kappa_khz,
        "g_MHz":                coupling_g_mhz,
        "resonator_type":       resonator_type,
    }

    # ========= RUN ANALYSIS =========
    results                                      = analyzer.find_closest(
        target_params=target_params,
        num_top=num_top,
        metric="Euclidean"
    )

    best_design                                  = results.iloc[0] 

    # ========= EXTRACT NESTED OPTIONS (MATCH JSON KEYS) =========     
    qubit_options                                = best_design["design_options_qubit"]     #CHECK TO MAKE THERE ARE NO MISSING VALUES MAY NEED TO TRANSFER TO ANOTEHR 
    resonator_options                            = best_design["design_options_cavity_claw"]["cpw_opts"]
    feedline_options                             = best_design["design_options_cavity_claw"]["cplr_opts"]  

    # ========= REF OBJECTS =========
    ref_design_specs                             = {
        "qubit_frequency_ghz":  best_design["qubit_frequency_GHz"],
        "cavity_frequency_ghz": best_design["cavity_frequency_GHz"],
        "anharmonicity_mhz":    best_design["anharmonicity_MHz"],
        "coupling_g_mhz":       best_design["g_MHz"],
        "kappa_khz":            best_design["kappa_kHz"],
        "resonator_type":       best_design["resonator_type"],
    }

    ref_capacitance_specs                        = {
        "EC":                   best_design["EC"],            # GHz
        "claw_to_claw":         best_design["claw_to_claw"],  # fF
        "claw_to_ground":       best_design["claw_to_ground"],
        "cross_to_claw":        best_design["cross_to_claw"],
        "cross_to_cross":       best_design["cross_to_cross"],
        "cross_to_ground":      best_design["cross_to_ground"],
        "ground_to_ground":     best_design["ground_to_ground"],
    }

    qubit_geo                                    = {
        "cross_gap":            qubit_options["cross_gap"],
        "cross_length":         qubit_options["cross_length"],
        "cross_width":          qubit_options["cross_width"],
        "claw_cpw_length":      qubit_options["connection_pads"]["readout"]["claw_cpw_length"],
        "claw_cpw_width":       qubit_options["connection_pads"]["readout"]["claw_cpw_width"],
        "claw_gap":             qubit_options["connection_pads"]["readout"]["claw_gap"],
        "claw_length":          qubit_options["connection_pads"]["readout"]["claw_length"],
        "claw_width":           qubit_options["connection_pads"]["readout"]["claw_width"],
        "connector_location":   qubit_options["connection_pads"]["readout"]["connector_location"],
        "connector_type":       qubit_options["connection_pads"]["readout"]["connector_type"],
        "ground_spacing":       qubit_options["connection_pads"]["readout"]["ground_spacing"],
    }

    resonator_geo                                = {
        "fillet":                   resonator_options["fillet"],
        "total_length":             resonator_options["total_length"],
        "trace_gap":                resonator_options["trace_gap"],
        "trace_width":              resonator_options["trace_width"],
        "end_straight":             resonator_options["lead"]["end_straight"],
        "start_jogged_extension":   resonator_options["lead"]["start_jogged_extension"],
        "start_straight":           resonator_options["lead"]["start_straight"],
        "asymmetry":                resonator_options["meander"]["asymmetry"],
        "spacing":                  resonator_options["meander"]["spacing"],
    }

    feedline_geo                                 = {
        "cap_distance":         feedline_options["cap_distance"],
        "cap_gap":              feedline_options["cap_gap"],
        "cap_gap_ground":       feedline_options["cap_gap_ground"],
        "cap_width":            feedline_options["cap_width"],
        "coupling_length":      feedline_options["coupling_length"],
        "coupling_space":       feedline_options["coupling_space"],
        "down_length":          feedline_options["down_length"],
        "finger_count":         feedline_options["finger_count"],
        "finger_length":        feedline_options["finger_length"],
        "open_termination":     feedline_options["open_termination"],
        "orientation":          feedline_options["orientation"],
        "prime_gap":            feedline_options["prime_gap"],
        "prime_width":          feedline_options["prime_width"],
        "second_gap":           feedline_options["second_gap"],
        "second_width":         feedline_options["second_width"],
    }

    Ej      = best_design["EJ"] * 1e9   # GHz
    Ej      = Ej* h                     # Convert Hz to Joules
    phi_0   = hbar / (2*e)
    LJs     = ((phi_0**2) / Ej)
    return best_design, ref_design_specs, ref_capacitance_specs, qubit_geo, resonator_geo, feedline_geo, LJs

def fetch_design_fast():
    """
    Uses only design_specs to query SQuADDS, extracts the top design,
    and returns: best_design, ref_design_specs, ref_capacitance_specs,
                 qubit_geo, resonator_geo, feedline_geo, LJs
    """
    best_design                                  = _load_best_design_from_json()

    # ========= EXTRACT NESTED OPTIONS (MATCH JSON KEYS) =========     
    qubit_options                                = best_design["design_options_qubit"]     #CHECK TO MAKE THERE ARE NO MISSING VALUES MAY NEED TO TRANSFER TO ANOTEHR 
    resonator_options                            = best_design["design_options_cavity_claw"]["cpw_opts"]
    feedline_options                             = best_design["design_options_cavity_claw"]["cplr_opts"]  

    # ========= REF OBJECTS =========
    ref_design_specs                             = {
        "qubit_frequency_ghz":  best_design["qubit_frequency_GHz"],
        "cavity_frequency_ghz": best_design["cavity_frequency_GHz"],
        "anharmonicity_mhz":    best_design["anharmonicity_MHz"],
        "coupling_g_mhz":       best_design["g_MHz"],
        "kappa_khz":            best_design["kappa_kHz"],
        "resonator_type":       best_design["resonator_type"],
    }

    ref_capacitance_specs                        = {
        "EC":                   best_design["EC"],            # GHz
        "claw_to_claw":         best_design["claw_to_claw"],  # fF
        "claw_to_ground":       best_design["claw_to_ground"],
        "cross_to_claw":        best_design["cross_to_claw"],
        "cross_to_cross":       best_design["cross_to_cross"],
        "cross_to_ground":      best_design["cross_to_ground"],
        "ground_to_ground":     best_design["ground_to_ground"],
    }

    qubit_geo                                    = {
        "cross_gap":            qubit_options["cross_gap"],
        "cross_length":         qubit_options["cross_length"],
        "cross_width":          qubit_options["cross_width"],
        "claw_cpw_length":      qubit_options["connection_pads"]["readout"]["claw_cpw_length"],
        "claw_cpw_width":       qubit_options["connection_pads"]["readout"]["claw_cpw_width"],
        "claw_gap":             qubit_options["connection_pads"]["readout"]["claw_gap"],
        "claw_length":          qubit_options["connection_pads"]["readout"]["claw_length"],
        "claw_width":           qubit_options["connection_pads"]["readout"]["claw_width"],
        "connector_location":   qubit_options["connection_pads"]["readout"]["connector_location"],
        "connector_type":       qubit_options["connection_pads"]["readout"]["connector_type"],
        "ground_spacing":       qubit_options["connection_pads"]["readout"]["ground_spacing"],
    }

    resonator_geo                                = {
        "fillet":                   resonator_options["fillet"],
        "total_length":             resonator_options["total_length"],
        "trace_gap":                resonator_options["trace_gap"],
        "trace_width":              resonator_options["trace_width"],
        "end_straight":             resonator_options["lead"]["end_straight"],
        "start_jogged_extension":   resonator_options["lead"]["start_jogged_extension"],
        "start_straight":           resonator_options["lead"]["start_straight"],
        "asymmetry":                resonator_options["meander"]["asymmetry"],
        "spacing":                  resonator_options["meander"]["spacing"],
    }

    feedline_geo                                 = {
        "cap_distance":         feedline_options["cap_distance"],
        "cap_gap":              feedline_options["cap_gap"],
        "cap_gap_ground":       feedline_options["cap_gap_ground"],
        "cap_width":            feedline_options["cap_width"],
        "coupling_length":      feedline_options["coupling_length"],
        "coupling_space":       feedline_options["coupling_space"],
        "down_length":          feedline_options["down_length"],
        "finger_count":         feedline_options["finger_count"],
        "finger_length":        feedline_options["finger_length"],
        "open_termination":     feedline_options["open_termination"],
        "orientation":          feedline_options["orientation"],
        "prime_gap":            feedline_options["prime_gap"],
        "prime_width":          feedline_options["prime_width"],
        "second_gap":           feedline_options["second_gap"],
        "second_width":         feedline_options["second_width"],
    }

    Ej      = best_design["EJ"] * 1e9  # GHz
    Ej      = Ej* h # Convert Hz to Joules
    phi_0   = hbar / (2*e)
    LJs     = ((phi_0**2) / Ej)
    return best_design, ref_design_specs, ref_capacitance_specs, qubit_geo, resonator_geo, feedline_geo, LJs

def _to_json_safe(obj):
    """Recursively convert pandas/NumPy objects to JSON-safe Python types."""
    if isinstance(obj, pd.DataFrame):
        return json.loads(obj.to_json(orient="records"))
    if isinstance(obj, pd.Series):
        return _to_json_safe(obj.to_dict())
    if isinstance(obj, np.ndarray):
        return [_to_json_safe(x) for x in obj.tolist()]
    if isinstance(obj, np.generic):         # numpy scalar
        return obj.item()
    if isinstance(obj, dict):
        return {str(k): _to_json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set)):
        return [_to_json_safe(v) for v in obj]
    return obj  # assume already JSON-serializable

def save_results_to_json(data, name="data/design_reference_data", save=True, orient="records", indent=2):
    """
    Save a pandas DataFrame/Series or a Python dict/list to JSON.
    """
    if not save:
        print("Save skipped (save=False).")
        return

    dirpath = os.path.dirname(name) or "."
    os.makedirs(dirpath, exist_ok=True)
    json_path = f"{name}.json"

    if isinstance(data, pd.DataFrame):
        data.to_json(json_path, orient=orient, indent=indent)
    else:
        safe = _to_json_safe(data)
        with open(json_path, "w") as f:
            json.dump(safe, f, indent=indent)

    print(f"✅ Saved to {json_path}")

def to_gds(design, filename="gds_files/unnamed_design.gds"):
    a_gds = design.renderers.gds

    # Get all unique layers in your design  
    layers = design.qgeometry.get_all_unique_layers('main')  
    
    # Create a dictionary with all layers set to False  
    layer_dict = {layer: False for layer in layers}  
    
    # Apply to both cheese and no_cheese  
    a_gds.options['cheese']['view_in_file']['main'] = layer_dict  
    a_gds.options['no_cheese']['view_in_file']['main'] = layer_dict
    a_gds.options['fabricate'] = True
    if os.path.exists(filename):
        os.remove(filename)
        
    a_gds.export_to_gds(f"{filename}")
    
def extract_um(value: str) -> float:
    """Extract a float value from a string like '30um'."""
    return float(value.replace('um', ''))

