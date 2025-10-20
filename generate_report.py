def _calculate_percentage_diff(value, reference):
    """
    Calculate percentage difference between a value and its reference.
    
    Parameters:
    -----------
    value : float
        The value to compare
    reference : float
        The reference value to compare against
    
    Returns:
    --------
    str
        Formatted percentage difference string with sign (e.g., "+2.34%" or "-1.23%")
    """
    if reference == 0:
        return "N/A"
    
    diff = ((value - reference) / reference) * 100
    
    # Format with sign
    if diff > 0:
        return f"+{diff:.2f}%"
    else:
        return f"{diff:.2f}%"

def _format_value(value, threshold=1e-3):
    """Format value: use scientific notation if |value| < threshold, else fixed-point."""
    if abs(value) < threshold and value != 0:
        return f"{value:.4e}"
    else:
        return f"{value:.4f}"
# ============================================================================
# LOM MATRICES (as strings for display)
# ============================================================================
# Helper function to format Hamiltonian matrix for display
def _format_hamiltonian_matrix(H, size=4, digits=2):
    """
    Format (size-1)x(size-1) corner of matrix H for LaTeX display.
    Always show real and imaginary components using scientific notation
    with proper LaTeX exponent formatting (e.g., 1.23×10^{-5}).
    If imag = 0, only show '+0j'; if both = 0, show '0'.
    """
    if H is None or not hasattr(H, "shape") or len(H.shape) != 2:
        return "N/A"
    if H.shape[0] == 0 or H.shape[1] == 0:
        return "N/A"

    nrows = min(size - 1, H.shape[0])
    ncols = min(size - 1, H.shape[1])

    def fmt_latex_sci(x: float) -> str:
        """Format number as LaTeX-style scientific notation."""
        if x == 0:
            return "0"
        mantissa, exp = f"{x:.{digits}e}".split("e")
        exp = int(exp)
        return f"{mantissa}×10^{{{exp}}}"

    out = []
    for i in range(nrows):
        row = []
        for j in range(ncols):
            try:
                z = complex(H[i, j])
            except Exception:
                z = 0.0 + 0.0j

            re, im = float(z.real), float(z.imag)

            if re == 0.0 and im == 0.0:
                cell = "0"
            else:
                re_s = fmt_latex_sci(re)
                if im == 0.0:
                    cell = f"{re_s}+0j"
                else:
                    im_s = fmt_latex_sci(abs(im))
                    sign = "+" if im >= 0 else "-"
                    cell = f"{re_s}{sign}{im_s}j"
            row.append(cell)

        row.append("\\cdots")
        out.append(" & ".join(row) + " \\\\")

    tail = ["\\vdots"] * max(1, ncols) + ["\\ddots"]
    out.append(" & ".join(tail))
    return "\n".join(out)

def generate_quantum_report(Lj, 
                            quantum_config,
                            best_design, 
                            results_inductex, 
                            palace_results, 
                            H_lom_hfss, 
                            H_lom_inductex, 
                            H_lom_palace, 
                            H_epr_inductex, 
                            H_epr_palace, 
                            results_lom_hfss,
                            results_lom_inductex,
                            results_lom_palace,
                            results_epr_inductex,
                            results_epr_palace,
                            compute_other_values,
                            output_file="Simulation Results and Analysis.html"
                            ):
    """
    Generate an HTML report for quantum simulation results.
    All values are defined as variables at the start for easy customization.
    """





    # ============================================================================
    # METADATA
    # ============================================================================
    creator                          = best_design["uploader_qubit"]
    institution                      = best_design["institution_qubit"]
    squadds_design_number            = str(best_design["index_qc"])
    design_type                      = "Qubit-Cavity System"
    resonator_type                   = best_design["resonator_type"]
    
    # ============================================================================
    # ELECTROSTATICS - HFSS (Reference)
    # ============================================================================
    hfss_Lj                          = Lj* 1e9
    hfss_Cj                          = compute_other_values["cj_hfss"]
    hfss_Cclaw_claw                  = best_design["claw_to_claw"]
    hfss_Cclaw_ground                = best_design["claw_to_ground"]
    hfss_Ccross_claw                 = best_design["cross_to_claw"]  # Cg
    hfss_Ccross_cross                = best_design["cross_to_cross"]  # Cs
    hfss_Ccross_ground               = best_design["cross_to_ground"]
    hfss_Cground_ground              = best_design["ground_to_ground"]
    hfss_Lr                          = results_inductex["inductance"]["resonator"]["L1"] * 1e9
    hfss_Cr                          = compute_other_values["cr_hfss"] * 1e15
    hfss_Cfeedline                   = 0
    
    # ============================================================================
    # ELECTROSTATICS - InductEx
    # ============================================================================
    inductex_Lj                      = Lj * 1e9
    inductex_Cj                      = 0 * 1e15
    inductex_Cclaw_claw              = results_inductex["capacitance"]["transmon"]["CCLAW-CCLAW"] * 1e15
    inductex_Cclaw_ground            = abs(results_inductex["capacitance"]["transmon"]["CGROUND-CCLAW"]) * 1e15
    inductex_Ccross_claw             = abs(results_inductex["capacitance"]["transmon"]["CCROSS-CCLAW"]) * 1e15  # Cg
    inductex_Ccross_cross            = results_inductex["capacitance"]["transmon"]["CCROSS-CCROSS"] * 1e15  # Cs
    inductex_Ccross_ground           = abs(results_inductex["capacitance"]["transmon"]["CGROUND-CCROSS"]) * 1e15
    inductex_Cground_ground          = results_inductex["capacitance"]["transmon"]["CGROUND-CGROUND"] * 1e15
    inductex_Lr                      = results_inductex["inductance"]["resonator"]["L1"] * 1e9
    inductex_Cr                      = results_inductex["capacitance"]["resonator"]["CRESONATOR-CRESONATOR"] * 1e15
    inductex_Cfeedline               = results_inductex["capacitance"]["feedline"]["CFEEDLINE-CFEEDLINE"] * 1e15
    
    # ============================================================================
    # ELECTROSTATICS - PALACE
    # ============================================================================
    palace_Lj                        = Lj* 1e9
    palace_Cj                        = 0 * 1e15
    palace_Cclaw_claw                = palace_results["capacitance"]["transmon"]["C3-C3"] * 1e15
    palace_Cclaw_ground              = abs(palace_results["capacitance"]["transmon"]["C1-C3"]) * 1e15
    palace_Ccross_claw               = abs(palace_results["capacitance"]["transmon"]["C2-C3"]) * 1e15  # Cg
    palace_Ccross_cross              = palace_results["capacitance"]["transmon"]["C2-C2"] * 1e15  # Cs
    palace_Ccross_ground             = abs(palace_results["capacitance"]["transmon"]["C1-C2"]) * 1e15
    palace_Cground_ground            = palace_results["capacitance"]["transmon"]["C1-C1"] * 1e15
    palace_Lr                        = compute_other_values["lr_palace"]
    palace_Cr                        = palace_results["capacitance"]["resonator"]["C2-C2"] * 1e15
    palace_Cfeedline                 = palace_results["capacitance"]["feedline"]["C2-C2"] * 1e15

    
    # ============================================================================
    # ELECTROSTATICS - Percentage Differences
    # ============================================================================
    # InductEx vs HFSS
    diff_inductex_Lj                 = _calculate_percentage_diff(inductex_Lj, hfss_Lj)
    diff_inductex_Cj                 = _calculate_percentage_diff(inductex_Cj, hfss_Cj) if hfss_Cj != 0 else "N/A"
    diff_inductex_Cclaw_claw         = _calculate_percentage_diff(inductex_Cclaw_claw, hfss_Cclaw_claw)
    diff_inductex_Cclaw_ground       = _calculate_percentage_diff(inductex_Cclaw_ground, hfss_Cclaw_ground)
    diff_inductex_Ccross_claw        = _calculate_percentage_diff(inductex_Ccross_claw, hfss_Ccross_claw)
    diff_inductex_Ccross_cross       = _calculate_percentage_diff(inductex_Ccross_cross, hfss_Ccross_cross)
    diff_inductex_Ccross_ground      = _calculate_percentage_diff(inductex_Ccross_ground, hfss_Ccross_ground)
    diff_inductex_Cground_ground     = _calculate_percentage_diff(inductex_Cground_ground, hfss_Cground_ground)
    diff_inductex_Lr                 = _calculate_percentage_diff(inductex_Lr, hfss_Lr) if hfss_Lr != 0 else "N/A"
    diff_inductex_Cr                 = _calculate_percentage_diff(inductex_Cr, hfss_Cr) if hfss_Cr != 0 else "N/A"
    diff_inductex_Cfeedline          = _calculate_percentage_diff(inductex_Cfeedline, hfss_Cfeedline) if hfss_Cfeedline != 0 else "N/A"
    
    # PALACE vs HFSS
    diff_palace_Lj                   = _calculate_percentage_diff(palace_Lj, hfss_Lj)
    diff_palace_Cj                   = _calculate_percentage_diff(palace_Cj, hfss_Cj) if hfss_Cj != 0 else "N/A"
    diff_palace_Cclaw_claw           = _calculate_percentage_diff(palace_Cclaw_claw, hfss_Cclaw_claw)
    diff_palace_Cclaw_ground         = _calculate_percentage_diff(palace_Cclaw_ground, hfss_Cclaw_ground)
    diff_palace_Ccross_claw          = _calculate_percentage_diff(palace_Ccross_claw, hfss_Ccross_claw)
    diff_palace_Ccross_cross         = _calculate_percentage_diff(palace_Ccross_cross, hfss_Ccross_cross)
    diff_palace_Ccross_ground        = _calculate_percentage_diff(palace_Ccross_ground, hfss_Ccross_ground)
    diff_palace_Cground_ground       = _calculate_percentage_diff(palace_Cground_ground, hfss_Cground_ground)
    diff_palace_Lr                   = _calculate_percentage_diff(palace_Lr, hfss_Lr) if hfss_Lr != 0 else "N/A"
    diff_palace_Cr                   = _calculate_percentage_diff(palace_Cr, hfss_Cr) if hfss_Cr != 0 else "N/A"
    diff_palace_Cfeedline            = _calculate_percentage_diff(palace_Cfeedline, hfss_Cfeedline) if hfss_Cfeedline != 0 else "N/A"
    
    # Mean differences for Cg and Cs
    # Calculate absolute values for mean
    abs_diff_inductex_cg = abs(((inductex_Ccross_claw - hfss_Ccross_claw) / hfss_Ccross_claw) * 100)
    abs_diff_inductex_cs = abs(((inductex_Ccross_cross - hfss_Ccross_cross) / hfss_Ccross_cross) * 100)
    mean_diff_inductex_cg_cs = f"{(abs_diff_inductex_cg + abs_diff_inductex_cs) / 2:.2f}%"
    
    abs_diff_palace_cg = abs(((palace_Ccross_claw - hfss_Ccross_claw) / hfss_Ccross_claw) * 100)
    abs_diff_palace_cs = abs(((palace_Ccross_cross - hfss_Ccross_cross) / hfss_Ccross_cross) * 100)
    mean_diff_palace_cg_cs = f"{(abs_diff_palace_cg + abs_diff_palace_cs) / 2:.2f}%"
    
    # ============================================================================
    # EIGENMODE ANALYSIS
    # ============================================================================
    # HFSS
    hfss_eigen_Fq                    = best_design["qubit_frequency_GHz"]
    hfss_eigen_Fc                    = best_design["cavity_frequency_GHz"]
    hfss_eigen_Pq                    = 0  # Not available in best_design
    hfss_eigen_Pc                    = 0  # Not available in best_design
    
    # InductEx
    inductex_eigen_Fq                = 0  # Not available in results_inductex
    inductex_eigen_Fc                = 0  # Not available in results_inductex
    inductex_eigen_Pq                = 0  # Not available in results_inductex
    inductex_eigen_Pc                = 0  # Not available in results_inductex
    
    # PALACE
    palace_eigen_Fq                  = palace_results["eigenmode"]["full"]["eigen_mode_frequencies"]["mode_1"]
    palace_eigen_Fc                  = palace_results["eigenmode"]["resonator"]["eigen_mode_frequencies"]["mode_1"]
    palace_eigen_Pq                  = abs(palace_results["eigenmode"]["full"]["eigen_mode_epr"]["mode_1"])
    palace_eigen_Pc                  = abs(palace_results["eigenmode"]["resonator"]["eigen_mode_epr"]["mode_1"])
    
    # Eigenmode percentage differences
    # InductEx vs HFSS
    diff_eigen_inductex_Fq           = _calculate_percentage_diff(inductex_eigen_Fq, hfss_eigen_Fq)
    diff_eigen_inductex_Fc           = _calculate_percentage_diff(inductex_eigen_Fc, hfss_eigen_Fc)
    diff_eigen_inductex_Pq           = _calculate_percentage_diff(inductex_eigen_Pq, hfss_eigen_Pq)
    diff_eigen_inductex_Pc           = _calculate_percentage_diff(inductex_eigen_Pc, hfss_eigen_Pc)
    
    # PALACE vs HFSS
    diff_eigen_palace_Fq             = _calculate_percentage_diff(palace_eigen_Fq, hfss_eigen_Fq)
    diff_eigen_palace_Fc             = _calculate_percentage_diff(palace_eigen_Fc, hfss_eigen_Fc)
    diff_eigen_palace_Pq             = _calculate_percentage_diff(palace_eigen_Pq, hfss_eigen_Pq)
    diff_eigen_palace_Pc             = _calculate_percentage_diff(palace_eigen_Pc, hfss_eigen_Pc)
    
    # ============================================================================
    # LOM MATRICES (as strings for display)
    # ============================================================================
    lom_hfss_matrix                  = _format_hamiltonian_matrix(H_lom_hfss)
    lom_inductex_matrix              = _format_hamiltonian_matrix(H_lom_inductex)
    lom_palace_matrix                = _format_hamiltonian_matrix(H_lom_palace)
    
    # ============================================================================
    # LOM HAMILTONIAN PARAMETERS
    # ============================================================================
    # Design (from best_design - these are target values)
    design_Fq                        = quantum_config["design_specs"]["qubit_frequency_ghz"]
    design_Fc                        = quantum_config["design_specs"]["cavity_frequency_ghz"]
    design_alpha                     = quantum_config["design_specs"]["anharmonicity_mhz"]
    design_kappa                     = quantum_config["design_specs"]["kappa_khz"]
    design_g                         = quantum_config["design_specs"]["coupling_g_mhz"]

    # HFSS Reference (from best_design)
    hfss_ref_Fq                      = best_design["qubit_frequency_GHz"]
    hfss_ref_Fc                      = best_design["cavity_frequency_GHz"]
    hfss_ref_alpha                   = best_design["anharmonicity_MHz"]
    hfss_ref_kappa                   = best_design["kappa_kHz"]
    hfss_ref_g                       = best_design["g_MHz"]
    
    # HFSS LOM
    hfss_lom_Fq                      = results_lom_hfss["qubit_frequency_ghz"]
    hfss_lom_Fc                      = results_lom_hfss["cavity_frequency_ghz"]
    hfss_lom_alpha                   = results_lom_hfss["anharmonicity_mhz"]
    hfss_lom_kappa                   = results_lom_hfss["kappa_khz"]
    hfss_lom_g                       = results_lom_hfss["coupling_g_mhz"]
    
    # InductEx LOM
    inductex_lom_Fq                  = results_lom_inductex["qubit_frequency_ghz"]
    inductex_lom_Fc                  = results_lom_inductex["cavity_frequency_ghz"]
    inductex_lom_alpha               = results_lom_inductex["anharmonicity_mhz"]
    inductex_lom_kappa               = results_lom_inductex["kappa_khz"]
    inductex_lom_g                   = results_lom_inductex["coupling_g_mhz"]
    
    # PALACE LOM
    palace_lom_Fq                    = results_lom_palace["qubit_frequency_ghz"]
    palace_lom_Fc                    = results_lom_palace["cavity_frequency_ghz"]
    palace_lom_alpha                 = results_lom_palace["anharmonicity_mhz"]
    palace_lom_kappa                 = results_lom_palace["kappa_khz"]
    palace_lom_g                     = results_lom_palace["coupling_g_mhz"]
    
    # ============================================================================
    # LOM PERCENTAGE DIFFERENCES
    # ============================================================================
    # HFSS Reference vs Design
    diff_hfss_ref_Fq                 = _calculate_percentage_diff(hfss_ref_Fq, design_Fq)
    diff_hfss_ref_Fc                 = _calculate_percentage_diff(hfss_ref_Fc, design_Fc)
    diff_hfss_ref_alpha              = _calculate_percentage_diff(hfss_ref_alpha, design_alpha)
    diff_hfss_ref_kappa              = _calculate_percentage_diff(hfss_ref_kappa, design_kappa)
    diff_hfss_ref_g                  = _calculate_percentage_diff(hfss_ref_g, design_g)

    # HFSS LOM vs Design / vs HFSS Reference
    diff_hfss_lom_Fq_design          = _calculate_percentage_diff(hfss_lom_Fq, design_Fq)
    diff_hfss_lom_Fq_ref             = _calculate_percentage_diff(hfss_lom_Fq, hfss_ref_Fq)
    diff_hfss_lom_Fc_design          = _calculate_percentage_diff(hfss_lom_Fc, design_Fc)
    diff_hfss_lom_Fc_ref             = _calculate_percentage_diff(hfss_lom_Fc, hfss_ref_Fc)
    diff_hfss_lom_alpha_design       = _calculate_percentage_diff(hfss_lom_alpha, design_alpha)
    diff_hfss_lom_alpha_ref          = _calculate_percentage_diff(hfss_lom_alpha, hfss_ref_alpha)
    diff_hfss_lom_kappa_design       = _calculate_percentage_diff(hfss_lom_kappa, design_kappa)
    diff_hfss_lom_kappa_ref          = _calculate_percentage_diff(hfss_lom_kappa, hfss_ref_kappa)
    diff_hfss_lom_g_design           = _calculate_percentage_diff(hfss_lom_g, design_g)
    diff_hfss_lom_g_ref              = _calculate_percentage_diff(hfss_lom_g, hfss_ref_g)

    # InductEx LOM vs Design / vs HFSS Reference
    diff_inductex_lom_Fq_design      = _calculate_percentage_diff(inductex_lom_Fq, design_Fq)
    diff_inductex_lom_Fq_ref         = _calculate_percentage_diff(inductex_lom_Fq, hfss_ref_Fq)
    diff_inductex_lom_Fc_design      = _calculate_percentage_diff(inductex_lom_Fc, design_Fc)
    diff_inductex_lom_Fc_ref         = _calculate_percentage_diff(inductex_lom_Fc, hfss_ref_Fc)
    diff_inductex_lom_alpha_design   = _calculate_percentage_diff(inductex_lom_alpha, design_alpha)
    diff_inductex_lom_alpha_ref      = _calculate_percentage_diff(inductex_lom_alpha, hfss_ref_alpha)
    diff_inductex_lom_kappa_design   = _calculate_percentage_diff(inductex_lom_kappa, design_kappa)
    diff_inductex_lom_kappa_ref      = _calculate_percentage_diff(inductex_lom_kappa, hfss_ref_kappa)
    diff_inductex_lom_g_design       = _calculate_percentage_diff(inductex_lom_g, design_g)
    diff_inductex_lom_g_ref          = _calculate_percentage_diff(inductex_lom_g, hfss_ref_g)

    # PALACE LOM vs Design / vs HFSS Reference
    diff_palace_lom_Fq_design        = _calculate_percentage_diff(palace_lom_Fq, design_Fq)
    diff_palace_lom_Fq_ref           = _calculate_percentage_diff(palace_lom_Fq, hfss_ref_Fq)
    diff_palace_lom_Fc_design        = _calculate_percentage_diff(palace_lom_Fc, design_Fc)
    diff_palace_lom_Fc_ref           = _calculate_percentage_diff(palace_lom_Fc, hfss_ref_Fc)
    diff_palace_lom_alpha_design     = _calculate_percentage_diff(palace_lom_alpha, design_alpha)
    diff_palace_lom_alpha_ref        = _calculate_percentage_diff(palace_lom_alpha, hfss_ref_alpha)
    diff_palace_lom_kappa_design     = _calculate_percentage_diff(palace_lom_kappa, design_kappa)
    diff_palace_lom_kappa_ref        = _calculate_percentage_diff(palace_lom_kappa, hfss_ref_kappa)
    diff_palace_lom_g_design         = _calculate_percentage_diff(palace_lom_g, design_g)
    diff_palace_lom_g_ref            = _calculate_percentage_diff(palace_lom_g, hfss_ref_g)
    
    # ============================================================================
    # EPR MATRICES (as strings for display)
    # ============================================================================
    epr_inductex_matrix              = _format_hamiltonian_matrix(H_epr_inductex)
    epr_palace_matrix                = _format_hamiltonian_matrix(H_epr_palace)
    
    # ============================================================================
    # EPR HAMILTONIAN PARAMETERS
    # ============================================================================
    # InductEx EPR
    inductex_epr_Fq                  = results_epr_inductex["qubit_frequency_ghz"]
    inductex_epr_Fc                  = results_epr_inductex["cavity_frequency_ghz"]
    inductex_epr_alpha               = results_epr_inductex["anharmonicity_mhz"]
    inductex_epr_kappa               = results_epr_inductex["kappa_khz"]
    inductex_epr_g                   = results_epr_inductex["coupling_g_mhz"]
    
    # PALACE EPR
    palace_epr_Fq                    = results_epr_palace["qubit_frequency_ghz"]
    palace_epr_Fc                    = results_epr_palace["cavity_frequency_ghz"]
    palace_epr_alpha                 = results_epr_palace["anharmonicity_mhz"]
    palace_epr_kappa                 = results_epr_palace["kappa_khz"]
    palace_epr_g                     = results_epr_palace["coupling_g_mhz"]
    
    # ============================================================================
    # EPR PERCENTAGE DIFFERENCES
    # ============================================================================
    # InductEx EPR vs Design / vs HFSS Reference
    diff_inductex_epr_Fq_design      = _calculate_percentage_diff(inductex_epr_Fq, design_Fq)
    diff_inductex_epr_Fq_ref         = _calculate_percentage_diff(inductex_epr_Fq, hfss_ref_Fq)
    diff_inductex_epr_Fc_design      = _calculate_percentage_diff(inductex_epr_Fc, design_Fc)
    diff_inductex_epr_Fc_ref         = _calculate_percentage_diff(inductex_epr_Fc, hfss_ref_Fc)
    diff_inductex_epr_alpha_design   = _calculate_percentage_diff(inductex_epr_alpha, design_alpha)
    diff_inductex_epr_alpha_ref      = _calculate_percentage_diff(inductex_epr_alpha, hfss_ref_alpha)
    diff_inductex_epr_kappa_design   = _calculate_percentage_diff(inductex_epr_kappa, design_kappa)
    diff_inductex_epr_kappa_ref      = _calculate_percentage_diff(inductex_epr_kappa, hfss_ref_kappa)
    diff_inductex_epr_g_design       = _calculate_percentage_diff(inductex_epr_g, design_g)
    diff_inductex_epr_g_ref          = _calculate_percentage_diff(inductex_epr_g, hfss_ref_g)

    # PALACE EPR vs Design / vs HFSS Reference
    diff_palace_epr_Fq_design        = _calculate_percentage_diff(palace_epr_Fq, design_Fq)
    diff_palace_epr_Fq_ref           = _calculate_percentage_diff(palace_epr_Fq, hfss_ref_Fq)
    diff_palace_epr_Fc_design        = _calculate_percentage_diff(palace_epr_Fc, design_Fc)
    diff_palace_epr_Fc_ref           = _calculate_percentage_diff(palace_epr_Fc, hfss_ref_Fc)
    diff_palace_epr_alpha_design     = _calculate_percentage_diff(palace_epr_alpha, design_alpha)
    diff_palace_epr_alpha_ref        = _calculate_percentage_diff(palace_epr_alpha, hfss_ref_alpha)
    diff_palace_epr_kappa_design     = _calculate_percentage_diff(palace_epr_kappa, design_kappa)
    diff_palace_epr_kappa_ref        = _calculate_percentage_diff(palace_epr_kappa, hfss_ref_kappa)
    diff_palace_epr_g_design         = _calculate_percentage_diff(palace_epr_g, design_g)
    diff_palace_epr_g_ref            = _calculate_percentage_diff(palace_epr_g, hfss_ref_g)

    # ============================================================================
    # GENERATE HTML
    # ============================================================================
    
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Simulation Results and Analysis</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            max-width: 1200px;
            margin: 0 auto;
            padding: 40px 20px;
            background: #0d1117;
            color: #c9d1d9;
        }}
        .container {{
            background: #161b22;
            padding: 40px;
            border-radius: 12px;
            box-shadow: 0 8px 24px rgba(0,0,0,0.4);
            border: 1px solid #30363d;
        }}
        h1 {{
            color: #58a6ff;
            border-bottom: 2px solid #30363d;
            padding-bottom: 15px;
            margin-bottom: 30px;
            font-size: 2.5em;
            font-weight: 600;
        }}
        h2 {{
            color: #58a6ff;
            margin-top: 40px;
            margin-bottom: 20px;
            font-size: 1.8em;
            font-weight: 600;
        }}
        .metadata {{
            background: #0d1117;
            padding: 20px;
            border-radius: 8px;
            margin-bottom: 30px;
            border: 1px solid #30363d;
        }}
        .metadata p {{
            margin: 8px 0;
            font-size: 1.1em;
            color: #8b949e;
        }}
        .metadata strong {{
            color: #c9d1d9;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 25px 0;
            border-radius: 8px;
            overflow: hidden;
            background: #0d1117;
            border: 1px solid #30363d;
        }}
        th {{
            background: #21262d;
            color: #58a6ff;
            padding: 15px;
            text-align: left;
            font-weight: 600;
            font-size: 1.1em;
            border-bottom: 1px solid #30363d;
        }}
        td {{
            padding: 12px 15px;
            border-bottom: 1px solid #21262d;
            color: #c9d1d9;
        }}
        tr:hover {{
            background-color: #161b22;
        }}
        tr:last-child td {{
            border-bottom: none;
        }}
        .reference {{
            background-color: #0d1117;
            font-weight: 500;
        }}
        .diff {{
            color: #58a6ff;
            font-weight: 600;
        }}
        .description {{
            background: #0d1117;
            padding: 20px;
            border-radius: 8px;
            margin: 25px 0;
            line-height: 1.8;
            color: #8b949e;
            border: 1px solid #30363d;
        }}
        .description p {{
            margin: 10px 0;
        }}
        .formula {{
            background: #0d1117;
            padding: 20px;
            border-radius: 8px;
            margin: 15px 0;
            text-align: center;
            font-size: 1.2em;
            border: 1px solid #30363d;
            color: #c9d1d9;
        }}
        .matrix {{
            background: #0d1117;
            padding: 25px;
            border-radius: 8px;
            margin: 20px 0;
            border: 1px solid #30363d;
        }}
        .matrix h3 {{
            color: #58a6ff;
            margin-top: 0;
            margin-bottom: 15px;
            font-weight: 600;
        }}
        .matrix-content {{
            font-family: 'Courier New', monospace;
            font-size: 1.1em;
            line-height: 2;
            background: #161b22;
            padding: 25px;
            border-radius: 6px;
            overflow-x: auto;
            color: #c9d1d9;
            border: 1px solid #30363d;
        }}
        .code-block {{
            background: #161b22;
            padding: 15px;
            border-radius: 6px;
            font-family: 'Courier New', monospace;
            font-size: 0.9em;
            color: #c9d1d9;
            border: 1px solid #30363d;
            margin: 10px 0;
        }}
    </style>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-MML-AM_CHTML"></script>
</head>
<body>
    <div class="container">
        <h1>Simulation Results and Analysis</h1>
        
        <div class="metadata">
            <p><strong>Design Created by:</strong> {creator} from {institution}</p>
            <p><strong>SQUaDDS Design Number:</strong> {squadds_design_number}</p>
            <p><strong>Design Type:</strong> {design_type}</p>
            <p><strong>Resonator Type:</strong> {resonator_type}</p>
        </div>

        <h2>Electrostatics</h2>
        <table>
            <thead>
                <tr>
                    <th>Parameter</th>
                    <th class="reference">HFSS (Reference)</th>
                    <th>InductEx</th>
                    <th>Δ%</th>
                    <th>PALACE</th>
                    <th>Δ%</th>
                </tr>
            </thead>
            <tbody>
                <tr>
                    <td><strong>L<sub>j</sub> (nH)</strong></td>
                    <td class="reference">{hfss_Lj:.4f}</td>
                    <td>{inductex_Lj:.4f}</td>
                    <td class="diff">{diff_inductex_Lj}</td>
                    <td>{palace_Lj:.4f}</td>
                    <td class="diff">{diff_palace_Lj}</td>
                </tr>
                <tr>
                    <td><strong>C<sub>j</sub> (fF)</strong></td>
                    <td class="reference">{_format_value(hfss_Cj)}*</td>
                    <td>{_format_value(inductex_Cj)}</td>
                    <td class="diff">{diff_inductex_Cj}</td>
                    <td>{_format_value(palace_Cj)}</td>
                    <td class="diff">{diff_palace_Cj}</td>
                </tr>
                <tr>
                    <td><strong>C<sub>claw-claw</sub> (fF)</strong></td>
                    <td class="reference">{hfss_Cclaw_claw:.4f}</td>
                    <td>{inductex_Cclaw_claw:.4f}</td>
                    <td class="diff">{diff_inductex_Cclaw_claw}</td>
                    <td>{palace_Cclaw_claw:.4f}</td>
                    <td class="diff">{diff_palace_Cclaw_claw}</td>
                </tr>
                <tr>
                    <td><strong>C<sub>claw-ground</sub> (fF)</strong></td>
                    <td class="reference">{hfss_Cclaw_ground:.4f}</td>
                    <td>{inductex_Cclaw_ground:.4f}</td>
                    <td class="diff">{diff_inductex_Cclaw_ground}</td>
                    <td>{palace_Cclaw_ground:.4f}</td>
                    <td class="diff">{diff_palace_Cclaw_ground}</td>
                </tr>
                <tr>
                    <td><strong>C<sub>cross-claw</sub> (C<sub>g</sub>) (fF)</strong></td>
                    <td class="reference">{hfss_Ccross_claw:.4f}</td>
                    <td>{inductex_Ccross_claw:.4f}</td>
                    <td class="diff">{diff_inductex_Ccross_claw}</td>
                    <td>{palace_Ccross_claw:.4f}</td>
                    <td class="diff">{diff_palace_Ccross_claw}</td>
                </tr>
                <tr>
                    <td><strong>C<sub>cross-cross</sub> (C<sub>s</sub>) (fF)</strong></td>
                    <td class="reference">{hfss_Ccross_cross:.4f}</td>
                    <td>{inductex_Ccross_cross:.4f}</td>
                    <td class="diff">{diff_inductex_Ccross_cross}</td>
                    <td>{palace_Ccross_cross:.4f}</td>
                    <td class="diff">{diff_palace_Ccross_cross}</td>
                </tr>
                <tr>
                    <td><strong>C<sub>cross-ground</sub> (fF)</strong></td>
                    <td class="reference">{hfss_Ccross_ground:.4f}</td>
                    <td>{inductex_Ccross_ground:.4f}</td>
                    <td class="diff">{diff_inductex_Ccross_ground}</td>
                    <td>{palace_Ccross_ground:.4f}</td>
                    <td class="diff">{diff_palace_Ccross_ground}</td>
                </tr>
                <tr>
                    <td><strong>C<sub>ground-ground</sub> (fF)</strong></td>
                    <td class="reference">{hfss_Cground_ground:.4f}</td>
                    <td>{inductex_Cground_ground:.4f}</td>
                    <td class="diff">{diff_inductex_Cground_ground}</td>
                    <td>{palace_Cground_ground:.4f}</td>
                    <td class="diff">{diff_palace_Cground_ground}</td>
                </tr>
                <tr>
                    <td><strong>L<sub>r</sub> (nH)</strong></td>
                    <td class="reference">{hfss_Lr:.4f} (InductEx)</td>
                    <td>{inductex_Lr:.4f}</td>
                    <td class="diff">{diff_inductex_Lr}</td>
                    <td>{palace_Lr:.4f}***</td>
                    <td class="diff">{diff_palace_Lr}</td>
                </tr>
                <tr>
                    <td><strong>C<sub>r</sub> (fF)</strong></td>
                    <td class="reference">{hfss_Cr:.4f}**</td>
                    <td>{inductex_Cr:.4f}</td>
                    <td class="diff">{diff_inductex_Cr}</td>
                    <td>{palace_Cr:.4f}</td>
                    <td class="diff">{diff_palace_Cr}</td>
                </tr>
                <tr>
                    <td><strong>C<sub>feedline</sub> (fF)</strong></td>
                    <td class="reference">{hfss_Cfeedline:.4f}</td>
                    <td>{inductex_Cfeedline:.4f}</td>
                    <td class="diff">{diff_inductex_Cfeedline}</td>
                    <td>{palace_Cfeedline:.4f}</td>
                    <td class="diff">{diff_palace_Cfeedline}</td>
                </tr>
            </tbody>
        </table>

        <div class="description">
            <p><strong>*</strong> Extracted using E<sub>c</sub>, C<sub>g</sub> and C<sub>s</sub>:</p>
            <div class="formula">
                \\( C_j = \\frac{{e^2}}{{2E_c}} - C_g - C_s \\approx 0 \\)
            </div>
            
            <p><strong>**</strong> Extracted using InductEx resonator inductance value as reference and HFSS reference cavity frequency:</p>
            <div class="formula">
                \\( C_r = \\frac{{1}}{{(2\\pi F_c)^2 \\cdot L_r}} \\)
            </div>
            
            <p><strong>***</strong> Extracted using PALACE resonator capacitance and PALACE cavity frequency:</p>
            <div class="formula">
                \\( L_r = \\frac{{1}}{{(2\\pi F_c)^2 \\cdot C_r}} \\)
            </div>
            
            <p style="margin-top: 20px;"><strong>Mean % difference for C<sub>g</sub> and C<sub>s</sub> for InductEx:</strong> {mean_diff_inductex_cg_cs}</p>
            <p><strong>Mean % difference for C<sub>g</sub> and C<sub>s</sub> for PALACE:</strong> {mean_diff_palace_cg_cs}</p>
        </div>

        <h2>Eigenmode Analysis</h2>
        <table>
            <thead>
                <tr>
                    <th>Mode Parameter</th>
                    <th>HFSS</th>
                    <th>InductEx</th>
                    <th>Δ%</th>
                    <th>PALACE</th>
                    <th>Δ%</th>
                </tr>
            </thead>
            <tbody>
                <tr>
                    <td><strong>Qubit Frequency (F<sub>q</sub>) (GHz)</strong></td>
                    <td>{hfss_eigen_Fq:.4f}</td>
                    <td>{inductex_eigen_Fq:.4f}</td>
                    <td class="diff">{diff_eigen_inductex_Fq}</td>
                    <td>{palace_eigen_Fq:.4f}</td>
                    <td class="diff">{diff_eigen_palace_Fq}</td>
                </tr>
                <tr>
                    <td><strong>Cavity Frequency (F<sub>c</sub>) (GHz)</strong></td>
                    <td>{hfss_eigen_Fc:.4f}</td>
                    <td>{inductex_eigen_Fc:.4f}</td>
                    <td class="diff">{diff_eigen_inductex_Fc}</td>
                    <td>{palace_eigen_Fc:.4f}</td>
                    <td class="diff">{diff_eigen_palace_Fc}</td>
                </tr>
                <tr>
                    <td><strong>Qubit Participation Ratio (P<sub>q</sub>) (Unitless)</strong></td>
                    <td>{hfss_eigen_Pq:.6f}</td>
                    <td>{inductex_eigen_Pq:.6f}</td>
                    <td class="diff">{diff_eigen_inductex_Pq}</td>
                    <td>{palace_eigen_Pq:.6f}</td>
                    <td class="diff">{diff_eigen_palace_Pq}</td>
                </tr>
                <tr>
                    <td><strong>Cavity Participation Ratio (P<sub>c</sub>) (Unitless)</strong></td>
                    <td>{hfss_eigen_Pc:.6f}</td>
                    <td>{inductex_eigen_Pc:.6f}</td>
                    <td class="diff">{diff_eigen_inductex_Pc}</td>
                    <td>{palace_eigen_Pc:.6f}</td>
                    <td class="diff">{diff_eigen_palace_Pc}</td>
                </tr>
            </tbody>
        </table>

        <h2>Local Oscillator Model (LOM): Hamiltonian Extraction and Analysis</h2>
        
        <div class="matrix">
            <h3>H<sub>HFSS</sub> [10×10]</h3>
            <div class="matrix-content">
\\[
\\begin{{bmatrix}}
{lom_hfss_matrix}
\\end{{bmatrix}}
\\]
            </div>
        </div>

        <div class="matrix">
            <h3>H<sub>InductEx</sub> [10×10]</h3>
            <div class="matrix-content">
\\[
\\begin{{bmatrix}}
{lom_inductex_matrix}
\\end{{bmatrix}}
\\]
            </div>
        </div>

        <div class="matrix">
            <h3>H<sub>PALACE</sub> [10×10]</h3>
            <div class="matrix-content">
\\[
\\begin{{bmatrix}}
{lom_palace_matrix}
\\end{{bmatrix}}
\\]
            </div>
        </div>

        <table>
            <thead>
                <tr>
                    <th>Parameter</th>
                    <th>Design</th>
                    <th>Δ%</th>
                    <th>HFSS Reference</th>
                    <th>Δ%</th>
                    <th>HFSS<sub>LOM</sub></th>
                    <th>Δ%</th>
                    <th>InductEx<sub>LOM</sub></th>
                    <th>Δ%</th>
                    <th>PALACE<sub>LOM</sub></th>
                    <th>Δ%</th>
                </tr>
            </thead>
            <tbody>
                <tr>
                    <td><strong>Qubit Frequency (F<sub>q</sub>) (GHz)</strong></td>
                    <td>{design_Fq:.4f}</td>
                    <td class="diff">-</td>
                    <td>{hfss_ref_Fq:.4f}</td>
                    <td class="diff">{diff_hfss_ref_Fq}</td>
                    <td>{hfss_lom_Fq:.4f}</td>
                    <td class="diff">{diff_hfss_lom_Fq_design}<br>{diff_hfss_lom_Fq_ref}</td>
                    <td>{inductex_lom_Fq:.4f}</td>
                    <td class="diff">{diff_inductex_lom_Fq_design}<br>{diff_inductex_lom_Fq_ref}</td>
                    <td>{palace_lom_Fq:.4f}</td>
                    <td class="diff">{diff_palace_lom_Fq_design}<br>{diff_palace_lom_Fq_ref}</td>
                </tr>
                <tr>
                    <td><strong>Cavity Frequency (F<sub>c</sub>) (GHz)</strong></td>
                    <td>{design_Fc:.4f}</td>
                    <td class="diff">-</td>
                    <td>{hfss_ref_Fc:.4f}</td>
                    <td class="diff">{diff_hfss_ref_Fc}</td>
                    <td>{hfss_lom_Fc:.4f}</td>
                    <td class="diff">{diff_hfss_lom_Fc_design}<br>{diff_hfss_lom_Fc_ref}</td>
                    <td>{inductex_lom_Fc:.4f}</td>
                    <td class="diff">{diff_inductex_lom_Fc_design}<br>{diff_inductex_lom_Fc_ref}</td>
                    <td>{palace_lom_Fc:.4f}</td>
                    <td class="diff">{diff_palace_lom_Fc_design}<br>{diff_palace_lom_Fc_ref}</td>
                </tr>
                <tr>
                    <td><strong>Anharmonicity (α) (MHz)</strong></td>
                    <td>{design_alpha:.2f}</td>
                    <td class="diff">-</td>
                    <td>{hfss_ref_alpha:.2f}</td>
                    <td class="diff">{diff_hfss_ref_alpha}</td>
                    <td>{hfss_lom_alpha:.2f}</td>
                    <td class="diff">{diff_hfss_lom_alpha_design}<br>{diff_hfss_lom_alpha_ref}</td>
                    <td>{inductex_lom_alpha:.2f}</td>
                    <td class="diff">{diff_inductex_lom_alpha_design}<br>{diff_inductex_lom_alpha_ref}</td>
                    <td>{palace_lom_alpha:.2f}</td>
                    <td class="diff">{diff_palace_lom_alpha_design}<br>{diff_palace_lom_alpha_ref}</td>
                </tr>
                <tr>
                    <td><strong>Kappa (κ) (kHz)</strong></td>
                    <td>{design_kappa:.2f}</td>
                    <td class="diff">-</td>
                    <td>{hfss_ref_kappa:.2f}</td>
                    <td class="diff">{diff_hfss_ref_kappa}</td>
                    <td>{hfss_lom_kappa:.2f}</td>
                    <td class="diff">{diff_hfss_lom_kappa_design}<br>{diff_hfss_lom_kappa_ref}</td>
                    <td>{inductex_lom_kappa:.2f}</td>
                    <td class="diff">{diff_inductex_lom_kappa_design}<br>{diff_inductex_lom_kappa_ref}</td>
                    <td>{palace_lom_kappa:.2f}</td>
                    <td class="diff">{diff_palace_lom_kappa_design}<br>{diff_palace_lom_kappa_ref}</td>
                </tr>
                <tr>
                    <td><strong>Coupling Strength (g) (MHz)</strong></td>
                    <td>{design_g:.2f}</td>
                    <td class="diff">-</td>
                    <td>{hfss_ref_g:.2f}</td>
                    <td class="diff">{diff_hfss_ref_g}</td>
                    <td>{hfss_lom_g:.2f}</td>
                    <td class="diff">{diff_hfss_lom_g_design}<br>{diff_hfss_lom_g_ref}</td>
                    <td>{inductex_lom_g:.2f}</td>
                    <td class="diff">{diff_inductex_lom_g_design}<br>{diff_inductex_lom_g_ref}</td>
                    <td>{palace_lom_g:.2f}</td>
                    <td class="diff">{diff_palace_lom_g_design}<br>{diff_palace_lom_g_ref}</td>
                </tr>
            </tbody>
        </table>

        <h2>Energy Participation Ratio (EPR): Hamiltonian Extraction and Analysis</h2>
        
        <div class="description">
            <p><em>Note: Design does not supply participation ratios</em></p>
        </div>

        <div class="matrix">
            <h3>H<sub>InductEx</sub> [10×10]</h3>
            <div class="matrix-content">
\\[
\\begin{{bmatrix}}
{epr_inductex_matrix}
\\end{{bmatrix}}
\\]
            </div>
        </div>

        <div class="matrix">
            <h3>H<sub>PALACE</sub> [10×10]</h3>
            <div class="matrix-content">
\\[
\\begin{{bmatrix}}
{epr_palace_matrix}
\\end{{bmatrix}}
\\]
            </div>
        </div>

        <table>
            <thead>
                <tr>
                    <th>Parameter</th>
                    <th>Design</th>
                    <th>Δ%</th>
                    <th>HFSS Reference</th>
                    <th>Δ%</th>
                    <th>InductEx<sub>EPR</sub></th>
                    <th>Δ%</th>
                    <th>PALACE<sub>EPR</sub></th>
                    <th>Δ%</th>
                </tr>
            </thead>
            <tbody>
                <tr>
                    <td><strong>Qubit Frequency (F<sub>q</sub>) (GHz)</strong></td>
                    <td>{design_Fq:.4f}</td>
                    <td class="diff">-</td>
                    <td>{hfss_ref_Fq:.4f}</td>
                    <td class="diff">{diff_hfss_ref_Fq}</td>
                    <td>{inductex_epr_Fq:.4f}</td>
                    <td class="diff">{diff_inductex_epr_Fq_design}<br>{diff_inductex_epr_Fq_ref}</td>
                    <td>{palace_epr_Fq:.4f}</td>
                    <td class="diff">{diff_palace_epr_Fq_design}<br>{diff_palace_epr_Fq_ref}</td>
                </tr>
                <tr>
                    <td><strong>Cavity Frequency (F<sub>c</sub>) (GHz)</strong></td>
                    <td>{design_Fc:.4f}</td>
                    <td class="diff">-</td>
                    <td>{hfss_ref_Fc:.4f}</td>
                    <td class="diff">{diff_hfss_ref_Fc}</td>
                    <td>{inductex_epr_Fc:.4f}</td>
                    <td class="diff">{diff_inductex_epr_Fc_design}<br>{diff_inductex_epr_Fc_ref}</td>
                    <td>{palace_epr_Fc:.4f}</td>
                    <td class="diff">{diff_palace_epr_Fc_design}<br>{diff_palace_epr_Fc_ref}</td>
                </tr>
                <tr>
                    <td><strong>Anharmonicity (α) (MHz)</strong></td>
                    <td>{design_alpha:.2f}</td>
                    <td class="diff">-</td>
                    <td>{hfss_ref_alpha:.2f}</td>
                    <td class="diff">{diff_hfss_ref_alpha}</td>
                    <td>{inductex_epr_alpha:.2f}</td>
                    <td class="diff">{diff_inductex_epr_alpha_design}<br>{diff_inductex_epr_alpha_ref}</td>
                    <td>{palace_epr_alpha:.2f}</td>
                    <td class="diff">{diff_palace_epr_alpha_design}<br>{diff_palace_epr_alpha_ref}</td>
                </tr>
                <tr>
                    <td><strong>Kappa (κ) (kHz)</strong></td>
                    <td>{design_kappa:.2f}</td>
                    <td class="diff">-</td>
                    <td>{hfss_ref_kappa:.2f}</td>
                    <td class="diff">{diff_hfss_ref_kappa}</td>
                    <td>{inductex_epr_kappa:.2f}</td>
                    <td class="diff">{diff_inductex_epr_kappa_design}<br>{diff_inductex_epr_kappa_ref}</td>
                    <td>{palace_epr_kappa:.2f}</td>
                    <td class="diff">{diff_palace_epr_kappa_design}<br>{diff_palace_epr_kappa_ref}</td>
                </tr>
                <tr>
                    <td><strong>Coupling Strength (g) (MHz)</strong></td>
                    <td>{design_g:.2f}</td>
                    <td class="diff">-</td>
                    <td>{hfss_ref_g:.2f}</td>
                    <td class="diff">{diff_hfss_ref_g}</td>
                    <td>{inductex_epr_g:.2f}</td>
                    <td class="diff">{diff_inductex_epr_g_design}<br>{diff_inductex_epr_g_ref}</td>
                    <td>{palace_epr_g:.2f}</td>
                    <td class="diff">{diff_palace_epr_g_design}<br>{diff_palace_epr_g_ref}</td>
                </tr>
            </tbody>
        </table>
    </div>
</body>
</html>"""
    
    # Write to file
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html)
    
    print(f"Report generated successfully: {output_file}")
    return output_file
