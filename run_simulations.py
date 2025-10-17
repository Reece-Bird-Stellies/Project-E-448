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

# No Used but could be useful 
def run_inductex_simulation_in_terminal(design_name, sim_type):
    """
    Run InductEx for a given design and sim type.
    Streams live output into the Jupyter Notebook.
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

    cmd = ["inductex8", control_map[sim_type]]
    print(f"[InductEx] Running in: {workdir}")
    print(f"[InductEx] Command: {' '.join(cmd)}\n")

    process = subprocess.Popen(
        cmd,
        cwd=workdir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )

    # Stream output live
    for line in process.stdout:
        print(line, end="")

    process.wait()
    print(f"\n[InductEx] Process exited with code {process.returncode}")
    return process.returncode

def run_palace_simulation_in_terminal(design_name, sim_type, palace_eigenmode_config):
    """
    Run PALACE in WSL for a given design and sim type.
    Streams output directly into the Jupyter Notebook.
    """
    number_of_cores = palace_eigenmode_config["number_of_cores"]
    config_file = f"{design_name}.json"
    win_dir = os.path.join(os.getcwd(), "simulations", "palace", sim_type, design_name)
    if not os.path.isdir(win_dir):
        raise FileNotFoundError(f"Simulation directory not found: {win_dir}")

    if not os.path.isfile(os.path.join(win_dir, config_file)):
        raise FileNotFoundError(f"Configuration file not found: {os.path.join(win_dir, config_file)}")

    # Convert Windows path to WSL format
    drive, rest = os.path.splitdrive(win_dir)
    wsl_dir = "/mnt/" + drive[0].lower() + rest.replace("\\", "/")

    bash_cmd = f'spack load palace@develop && cd "{wsl_dir}" && palace -np {number_of_cores} {config_file}'
    cmd = ["wsl", "bash", "-ic", bash_cmd]

    print(f"[PALACE] Running in: {win_dir}")
    print(f"[PALACE] Command: {' '.join(cmd)}\n")

    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )

    # Stream output live
    for line in process.stdout:
        print(line, end="")

    process.wait()
    print(f"\n[PALACE] Process exited with code {process.returncode}")
    return process.returncode


    """
    Run InductEx in a new CMD window and close that CMD window when the command finishes.
    Uses: cmd /c start "" cmd /c "your_command && exit"
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
    print(f"[InductEx] Launching CMD in: {workdir}")
    print(f"[InductEx] Command: {cmd_to_run}")

    # Make sure cmd.exe is available
    if shutil.which("cmd") is None:
        raise EnvironmentError("cmd.exe not found in PATH (this must run on Windows).")

    # Use cmd to start a new window which runs cmd /c "command && exit"
    # Note: start "" sets the window title to empty string (prevents confusion with the command)
    try:
        subprocess.Popen(
            ["cmd", "/c", "start", "", "cmd", "/c", f"{cmd_to_run} && exit"],
            cwd=workdir,
            shell=False
        )
        print("[InductEx] Simulation started (new CMD window). Window will close on completion.")
    except Exception as e:
        print("[InductEx][ERROR] Failed to start CMD window:", e)
        # Offer fallback: run in the current process (blocking) so user sees errors
        print("[InductEx] Fallback: running command in-process (will block until completion)")
        subprocess.run(cmd_to_run, cwd=workdir, shell=True, check=True)