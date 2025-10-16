def _get_mask_trace_string(mask_bounds, scale_factor=1.0):
    """
    Converts [x1, y1, x2, y2] (bottom-left to top-right)
    into a formatted XY polygon trace string for IXI files.
    Optionally scales about the rectangle's centre.
    """

    x1, y1, x2, y2 = mask_bounds

    # --- scale relative to centre ---
    if scale_factor != 1.0:
        cx = (x1 + x2) / 2.0
        cy = (y1 + y2) / 2.0
        x1 = cx + scale_factor * (x1 - cx)
        y1 = cy + scale_factor * (y1 - cy)
        x2 = cx + scale_factor * (x2 - cx)
        y2 = cy + scale_factor * (y2 - cy)

    # --- define closed rectangle (BL → BR → TR → TL → BL) ---
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

def _make_ldf_file(param_string, layers_list, output_path):
    layers_string = "\n".join(layers_list)
    ldf_contents = f"""* -------------------------------------------
* PARAMETERS
* -------------------------------------------
$Parameters
    {param_string}
$End
* -------------------------------------------
* LAYERS
* -------------------------------------------

{layers_string}
"""
    # --- Write to output path ---
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(ldf_contents.strip() + "\n")

def _make_cir_file():
    pass

def _make_cap_run_file(design_name):

    run_contents = f"""inductex8 control_{design_name}_cap.ixi"""

    output_path = rf"simulations/inductex/capacitance/{design_name}/run_cap_sim.cmd"

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(run_contents.strip() + "\n")

def _make_cap_ixi_file(design_name, mask_bounds, inductex_ixi_config):
    design_name_cap = design_name + "_cap"
    param_string = f"""$PARAMSTRING
    -l process_{design_name_cap}.ldf
    -qt 
    -zoff
$END"""

    if mask_bounds is None:
        mask_string = ""
    else:
        mask_string = f"""$MASK
    MeshSize {inductex_ixi_config['mask_segment_size']}
    XY
        {_get_mask_trace_string(mask_bounds, inductex_ixi_config["mask_scaling_factor"])}
$END"""
    
    output_path = rf"simulations/inductex/capacitance/{design_name}/control_{design_name_cap}.ixi"

    _make_ixi_file(param_string, design_name_cap, mask_string, output_path)

def _make_cap_ldf_file(design_name, inductex_ldf_config):
    design_name_cap = design_name + "_cap"

    if inductex_ldf_config["substrate_choice"] == "silicon":
        epsilon = 11.45
    else: 
        epsilon = 9.8  # sapphire

    param_string = f"""Units 					= 1e-6
	SegmentSize 			= {inductex_ldf_config['segment_size_um']}
	GPOverhang 				= 5.0
	AbsMin 					= 0.025
	ProcessHasGroundPlane 	= False
	TextLayer 				= 100
	TermLayer 				= 110
	TerminalInRange 		= 0.2
	UnitL 					= 1
	HFilaments 				= {inductex_ldf_config["h_filaments_metal"]}
	PolyDecimation 			= True 
	DecimationDistance 		= {inductex_ldf_config["decimation_distance"]}
    Lambda0					= {inductex_ldf_config["lambda0_um"]}
	LambdaFit				= 4
	Temperature		 		= {inductex_ldf_config["temperature_k"]}
	TC						= {inductex_ldf_config["critical_temperature"]}"""

    substrate_layer_string = f"""* ------------Silicon Base Layer-------------
$Layer
    Number 					= 90
	Name 					= sil
	Thickness 				= 750
	EpsilonR 				= {epsilon}
	MuR 					= 1
	Order 					= 0
	Mask 					= -1
	Filmtype 				= I
$End"""

    metal_layer_string = f"""* -----------Main Metal Layer----------------
$Layer
    Number 					= 0
	Name 					= main
	Thickness 				= {inductex_ldf_config["thickness_um"]}
	Order 					= 1
	Mask 					= 1
	Filmtype 				= S
$End"""
    layers_list = [substrate_layer_string, metal_layer_string]
    output_path = rf"simulations/inductex/capacitance/{design_name}/process_{design_name_cap}.ldf"
    _make_ldf_file(param_string, layers_list, output_path)

def make_inductex_cap_sim(design, design_name, create_mask, mask_component=None, inductex_ixi_config=None, inductex_ldf_config=None):
    if create_mask and mask_component:
        mask_bounds = design.components[mask_component].qgeometry_bounds()
        mask_bounds = mask_bounds * 1e3  # convert to um
    else:
        mask_bounds = None
    _make_cap_ixi_file(design_name, mask_bounds, inductex_ixi_config or {})
    _make_cap_ldf_file(design_name, inductex_ldf_config or {})
    _make_cap_run_file(design_name)