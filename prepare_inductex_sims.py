def _get_mask_trace_string(mask_bounds):
    """
    Converts [x1, y1, x2, y2] (bottom-left to top-right)
    into a formatted XY polygon trace string for IXI files.
    """
    x1, y1, x2, y2 = mask_bounds

    points = [
        (x1, y1),  # bottom-left
        (x2, y1),  # bottom-right
        (x2, y2),  # top-right
        (x1, y2),  # top-left
        (x1, y1)   # close loop
    ]

    def fmt(v): return f"{v:11.5f}"

    return "\n        ".join(f"({fmt(x)}, {fmt(y)})" for x, y in points)

def _make_ixi_file(param_string, design_name_sim, mask_string, output_path):
    """
    Creates a .IXI file content string and writes it to disk.
    param_string : str
        Full PARAMSTRING section contents (already formatted)
    gds_name : str
        Name of the GDS file (e.g., "qubit.gds")
    mask_bounds : dict
        {
            "mesh_size": 20,
            "top_left": (-1236.9, 348.45),
            "bottom_right": (-10.95, 1051.4),
            "outfile": "cap.ixi"
        }
    """

    # --- Get formatted mask coordinates ---
    # --- Build file content ---
    ixi_contents = f"""
{param_string}

$STRUCT
    Name cap
    $GDS
        name {design_name_sim}.gds
    $END
$END

{mask_string}"""



    # --- Write to output path ---
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(ixi_contents.strip() + "\n")

    return output_path

def _make_ldf_file():
    pass

def _make_cir_file():
    pass

def _make_cap_ixi_file(design_name, mask_bounds, inductex_ixi_config):
    design_name_cap = design_name + "_cap"
    param_string = f"""$PARAMSTRING
    -l process_{design_name_cap}.gds
    -qt 
    -zoff
$END"""

    if mask_bounds is None:
        mask_string = ""
    else:
        mask_string = f"""$MASK
    MeshSize {inductex_ixi_config['mask_segment_size']}
    XY
        {_get_mask_trace_string(mask_bounds*inductex_ixi_config["mask_scaling_factor"])}
$END"""
    
    output_path = rf"simulations/inductex/capacitance/{design_name}/control_{design_name_cap}.ixi"
    
    return _make_ixi_file(param_string, design_name_cap, mask_string, output_path)

def _make_cap_ldf_file():
    pass

def make_inductex_cap_sim(design, design_name, create_mask, mask_component=None, inductex_ixi_config=None, inductex_ldf_config=None):
    if create_mask and mask_component:
        mask_bounds = design.components[mask_component].qgeometry_bounds()
        mask_bounds = mask_bounds * 1e3  # convert to um
    else:
        mask_bounds = None
    _make_cap_ixi_file(design_name, mask_bounds, inductex_ixi_config or {})
    _make_cap_ldf_file(design_name, inductex_ldf_config or {})