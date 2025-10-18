import json
import os
import re

# InductEx Collection Functions
def _get_inductex_time(file_path):
    """
    Extracts the execution time from an InductEx sol.txt file.
    
    Args:
        file_path: Path to the sol.txt file
    
    Returns:
        Float of execution time in seconds, or None if not found
    """
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
            # Check the last few lines for the "Job finished" message
            for line in reversed(lines[-10:]):  # Check last 10 lines
                match = re.search(r'Job finished in ([\d.]+) seconds', line)
                if match:
                    return float(match.group(1))
        return None
    except (FileNotFoundError, IOError) as e:
        return None

def _get_inductex_capacitance_results(file_path):
    """
    Extracts capacitance results from an InductEx elements.txt file.
    
    Args:
        file_path: Path to the elements.txt file
    
    Returns:
        Dictionary of capacitance results (Maxwell only), or None if file not found or invalid
    """
    try:
        with open(file_path, 'r') as file:
            capacitance_results = {}
            for line in file:
                # Check if line contains Maxwell capacitance (ignore SPICE)
                if '(Maxwell)' in line:
                    parts = line.split()
                    if len(parts) >= 2:
                        #parts[0] is the capacitor name (e.g., "CGROUND-CCROSS(Maxwell)")
                        # parts[1] is the value
                        cap_name = parts[0].replace('(Maxwell)', '')  # Remove (Maxwell) suffix
                        cap_value = float(parts[1])
                        capacitance_results[cap_name] = cap_value
            
            return capacitance_results if capacitance_results else None
    except (FileNotFoundError, ValueError, IOError) as e:
        return None

def _get_inductex_inductance_results(file_path):
    """
    Extracts inductance results from an InductEx elements.txt file.
    
    Args:
        file_path: Path to the elements.txt file
    
    Returns:
        Dictionary of inductance results, or None if file not found or invalid
    """
    try:
        with open(file_path, 'r') as file:
            inductance_results = {}
            for line in file:
                parts = line.split()
                # Check if line starts with 'L' and has a value
                if len(parts) >= 2 and parts[0].startswith('L'):
                    inductor_name = parts[0]
                    inductor_value = float(parts[1])
                    inductance_results[inductor_name] = inductor_value
            
            return inductance_results if inductance_results else None
    except (FileNotFoundError, ValueError, IOError) as e:
        return None
    
def _get_inductex_sparams_results(file_path):
    """
    Extracts S-parameter results from an InductEx sparams.m file.
    
    Args:
        file_path: Path to the sparams.m file
    
    Returns:
        Dictionary of S-parameter results and frequencies, or None if file not found or invalid
    """
    try:
        with open(file_path, 'r') as file:
            sparams_results = {}
            for line in file:
                # Skip comments and empty lines
                if line.strip().startswith('%') or not line.strip():
                    continue
                
                # Check if line contains S-parameter or frequency data
                if '=' in line and ';' in line:
                    parts = line.split('=')
                    if len(parts) == 2:
                        param_name = parts[0].strip()
                        # Extract the values between [ and ]
                        values_str = parts[1].strip().rstrip(';')
                        if values_str.startswith('[') and values_str.endswith(']'):
                            values_str = values_str[1:-1]  # Remove brackets
                            # Convert comma-separated values to list of floats
                            values = [float(v.strip()) for v in values_str.split(',')]
                            sparams_results[param_name] = values
            
            return sparams_results if sparams_results else None
    except (FileNotFoundError, ValueError, IOError) as e:
        return None
    
def collect_inductex_results(design_names, inductex_simulation_types):
    """
    Collects InductEx simulation results including execution times and specific results
    based on the simulation type.
    
    Args:
        design_names: List of design names
        inductex_simulation_types: List of InductEx simulation types
    
    Returns:
        Dictionary with nested structure: {simulation_type: {design_name: execution_time}}
        for use with print_simulation_times, or full results if needed
    """
    sim_durations = {}
    base_path = os.getcwd()
    results = {}
    
    for sim_type in inductex_simulation_types:
        sim_durations[sim_type] = {}
        results[sim_type] = {}  # ADD THIS LINE

        for design_name in design_names:
            # Construct the base path for the simulation type and design
            path = os.path.join(
                base_path, "simulations", "inductex", sim_type, design_name, "output"
            )
            # Collect results based on simulation type
            if sim_type == "capacitance":
                results[sim_type][design_name] = _get_inductex_capacitance_results(
                    os.path.join(path, "elements.txt")
                )
            elif sim_type == "inductance":
                results[sim_type][design_name] = _get_inductex_inductance_results(
                    os.path.join(path, "elements.txt")
                )
            elif sim_type == "sparams":
                results[sim_type][design_name] = _get_inductex_sparams_results(
                    os.path.join(path, "sparams.m")
                )
            exec_time = _get_inductex_time(os.path.join(path, "sol.txt"))
            sim_durations[sim_type][design_name] = exec_time

    return results, sim_durations

# PALACE Collection Functions
def _get_palace_time(file_path):
    try:
        with open(file_path, 'r') as file:
            data = json.load(file)
            return data.get("ElapsedTime", {}).get("Durations", {}).get("Total", None)
    except (FileNotFoundError, json.JSONDecodeError) as e:
        return None
    
def _get_palace_capacitance_results(file_path):
    """
    Extracts capacitance results from a Palace terminal-C.txt file.
    
    Args:
        file_path: Path to the terminal-C.txt file
    
    Returns:
        Dictionary of capacitance results in the format {name: value}, or None if file not found or invalid
    """
    try:
        with open(file_path, 'r') as file:
            capacitance_results = {}
            lines = file.readlines()
            
            # Find the header line to determine number of columns
            header_line = None
            data_start_idx = 0
            
            for idx, line in enumerate(lines):
                if 'C[i][' in line:
                    header_line = line
                    data_start_idx = idx + 1
                    break
            
            if header_line is None:
                return None
            
            # Parse data lines
            for line in lines[data_start_idx:]:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split(',')
                if len(parts) < 2:
                    continue
                
                try:
                    # First column is the row index
                    row_idx = int(float(parts[0].strip()))
                    
                    # Process each capacitance value in the row
                    for col_idx, value_str in enumerate(parts[1:], start=1):
                        value = float(value_str.strip())
                        
                        # Create a name for this capacitance element
                        if row_idx == col_idx:
                            # Diagonal element: C[i][i]
                            cap_name = f"C{row_idx}-C{row_idx}"
                        else:
                            # Off-diagonal element: C[i][j] (only store for i < j to avoid duplicates)
                            if row_idx < col_idx:
                                cap_name = f"C{row_idx}-C{col_idx}"
                            else:
                                continue  # Skip lower triangle
                        
                        capacitance_results[cap_name] = value
                
                except (ValueError, IndexError):
                    continue
            
            return capacitance_results if capacitance_results else None
    except (FileNotFoundError, IOError) as e:
        return None

def _get_palace_eigenmode_results(file_path):
    """
    Extracts eigenmode and EPR results from Palace eig.txt and port-EPR.txt files.
    
    Args:
        file_path: Path to the directory containing eig.txt and port-EPR.txt files
    
    Returns:
        Dictionary containing eigenmode frequencies and EPR values, or None if files not found or invalid
    """
    eigenmode_results = {
        "eigen_mode_frequencies": {},
        "eigen_mode_epr": {}
    }

    # Get eigenmode frequencies from eig.csv
    eig_file = os.path.join(file_path, "eig.csv")
    try:
        with open(eig_file, 'r') as file:
            for line in file:
                line = line.strip()
                if not line or line.startswith('#') or 'm,' in line:
                    continue
                
                parts = line.split(',')
                if len(parts) >= 2:
                    try:
                        # First column is mode number, second is Re{f} in GHz
                        mode_num = int(float(parts[0].strip()))
                        freq_ghz = float(parts[1].strip())
                        
                        mode_name = f"mode_{mode_num}"
                        eigenmode_results["eigen_mode_frequencies"][mode_name] = freq_ghz
                    except (ValueError, IndexError):
                        continue
    except (FileNotFoundError, IOError):
        pass

    # Get EPR values from port-EPR.csv
    epr_file = os.path.join(file_path, "port-EPR.csv")
    try:
        with open(epr_file, 'r') as file:
            for line in file:
                line = line.strip()
                if not line or line.startswith('#') or 'm,' in line:
                    continue
                
                parts = line.split(',')
                if len(parts) >= 2:
                    try:
                        # First column is mode number, second is EPR value
                        mode_num = int(float(parts[0].strip()))
                        epr_value = float(parts[1].strip())
                        
                        mode_name = f"mode_{mode_num}"
                        eigenmode_results["eigen_mode_epr"][mode_name] = epr_value
                    except (ValueError, IndexError):
                        continue
    except (FileNotFoundError, IOError):
        pass
    
    # Return None if no data was collected
    if not eigenmode_results["eigen_mode_frequencies"] and not eigenmode_results["eigen_mode_epr"]:
        return None
    
    return eigenmode_results

def collect_palace_results(design_names, palace_simulation_types):
    """
    Collects total times for all palace simulations.
    
    Args:
        design_names: List of design names
        palace_simulation_types: List of simulation types
    
    Returns:
        Dictionary with nested structure: {simulation_type: {design_name: total_time}}
    """
    sim_durations = {}
    base_path = os.getcwd()
    results = {}
    
    for sim_type in palace_simulation_types:
        sim_durations[sim_type] = {}
        results[sim_type] = {}  # ADD THIS LINE

        for design_name in design_names:
            # Construct the base path for the simulation type and design
            path = os.path.join(
                base_path, "simulations", "palace", sim_type, design_name, "outputFiles"
            )
            # Collect results based on simulation type
            if sim_type == "capacitance":
                results[sim_type][design_name] = _get_palace_capacitance_results(os.path.join(path, "terminal-C.csv"))
            elif sim_type == "eigenmode":
                results[sim_type][design_name] = _get_palace_eigenmode_results(path)

            exec_time = _get_palace_time(os.path.join(path, "palace.json"))
            sim_durations[sim_type][design_name] = exec_time

    return results, sim_durations