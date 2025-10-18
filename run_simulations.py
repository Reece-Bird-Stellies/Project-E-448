import os
import subprocess

def run_inductex_simulation(design_name, sim_type, keep_open=False):
    """
    Run InductEx in a new CMD window and optionally keep it open after completion.
    
    Args:
        design_name: Name of the design
        sim_type: Type of simulation ('capacitance', 'inductance', or 'sparams')
        keep_open: If True, keeps terminal open after completion. Default False.
    """
    control_map = {
        "capacitance": f"control_{design_name}_cap.ixi",
        "inductance":  f"control_{design_name}_ind.ixi",
        "sparams":     f"control_{design_name}_sparams.ixi",
    }
    if sim_type not in control_map:
        raise ValueError(f"Unknown sim_type: {sim_type}")

    workdir = os.path.join(os.getcwd(), "simulations", "inductex", sim_type, design_name)
    if not os.path.isdir(workdir):
        raise FileNotFoundError(f"Simulation directory not found: {workdir}")

    control_file = control_map[sim_type]
    cmd_to_run = f'inductex8 {control_file}'
    
    print(f"[InductEx] Running in: {workdir}")
    
    # Build command based on keep_open parameter
    if keep_open:
        full_command = f"{cmd_to_run} && echo. && echo Simulation complete. Press any key to exit... && pause > nul"
    else:
        full_command = f"{cmd_to_run} && exit"
    
    process = subprocess.Popen(
        ["cmd", "/c", full_command],
        cwd=workdir,
        shell=False,
        creationflags=subprocess.CREATE_NEW_CONSOLE
    )
    
    # Wait for process and check return code
    process.wait()
    if process.returncode == 0:
        print(f"InductEx {sim_type} simulation of {design_name} was successful")
    else:
        print(f"InductEx {sim_type} simulation of {design_name} was unsuccessful")

def run_palace_simulation(design_name, simulation, number_of_cores, keep_open=False):
    """
    Opens an Ubuntu WSL terminal, converts the Windows path to Linux format,
    runs PALACE commands, and optionally keeps the terminal open.
    
    Args:
        design_name: Name of the design
        simulation: Simulation type
        number_of_cores: Number of cores to use
        keep_open: If True, keeps terminal open after completion. Default False.
    """
    cwd = os.getcwd()
    workdir = os.path.join(cwd, "simulations", "palace", simulation, design_name)

    if not os.path.isdir(workdir):
        raise FileNotFoundError(f"Simulation directory not found: {workdir}")

    # Convert Windows path -> WSL path (/mnt/c/…)
    drive, path_rest = os.path.splitdrive(workdir)
    linux_path = f"/mnt/{drive[0].lower()}{path_rest.replace(os.sep, '/')}"

    # Command sequence - with optional pause at end
    if keep_open:
        wsl_command = (
            f"source $HOME/spack/share/spack/setup-env.sh; "
            f"cd '{linux_path}'; "
            f"spack load palace@develop; "
            f"palace -np {number_of_cores} {design_name}.json; "
            f"echo ''; "
            f"echo 'Simulation complete. Press any key to exit...'; "
            f"read -n 1"
        )
    else:
        wsl_command = (
            f"source $HOME/spack/share/spack/setup-env.sh; "
            f"cd '{linux_path}'; "
            f"spack load palace@develop; "
            f"palace -np {number_of_cores} {design_name}.json"
        )
    
    print(f"[PALACE] Running in: {linux_path}")

    process = subprocess.Popen(
        ["wsl", "-d", "Ubuntu", "bash", "-c", wsl_command],
        creationflags=subprocess.CREATE_NEW_CONSOLE
    )
    
    # Wait for process and check return code
    process.wait()
    if process.returncode == 0:
        print(f"PALACE {simulation} simulation of {design_name} was successful")
    else:
        print(f"PALACE {simulation} simulation of {design_name} was unsuccessful")