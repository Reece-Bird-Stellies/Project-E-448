import os
import pya
import numpy as np

from helpers import extract_um

def _convert_layer_to_dtype0(ly, layer_number):
    """
    Converts all datatypes of a given layer to datatype 0.
    Moves all shapes (polygons, paths, boxes) from (layer, any dtype) → (layer, 0)
    and deletes the original datatypes.
    """
    dst = ly.layer(pya.LayerInfo(layer_number, 0))
    for li in list(ly.layer_indexes()):
        info = ly.get_info(li)
        if info.layer == layer_number and info.datatype != 0:
            for c in ly.each_cell():
                reg = pya.Region(c.begin_shapes_rec(li))
                if not reg.is_empty():
                    c.shapes(dst).insert(reg)
                    c.shapes(li).clear()
            ly.delete_layer(li)

def _boolean_layer_operation(ly, top, A, B, C, operation="AND"):
    """
    Perform a boolean operation (AND or OR) between two layers (A, B)
    and write the result to output layer C.

    Parameters:
        ly : pya.Layout
            Loaded KLayout layout.
        top : pya.Cell
            Top cell of the layout.
        A, B, C : tuple(int, int)
            (layer_number, datatype) for input/output layers.
        operation : str
            Either "AND" or "OR" (case-insensitive).
    """
    LA = ly.find_layer(pya.LayerInfo(A[0], A[1]))
    LB = ly.find_layer(pya.LayerInfo(B[0], B[1]))

    # Create Regions
    r = pya.Region()
    if LA >= 0:
        r = pya.Region(top.begin_shapes_rec(LA))
    if LB >= 0:
        if operation.upper() == "AND":
            r &= pya.Region(top.begin_shapes_rec(LB))
        elif operation.upper() == "OR":
            r |= pya.Region(top.begin_shapes_rec(LB))
        else:
            raise ValueError("operation must be 'AND' or 'OR'")

    # Create destination layer C and clear any existing shapes
    LC = ly.layer(pya.LayerInfo(C[0], C[1]))
    for c in ly.each_cell():
        c.shapes(LC).clear()

    # Insert the boolean result into output layer on the top cell
    top.shapes(LC).insert(r)

    print(f"Boolean {operation.upper()}({A}, {B}) → {C} complete.")

def _delete_layer(ly, layer_tuple):
    """
    Deletes a specific layer (and datatype) from the layout if it exists.

    Parameters:
        ly : pya.Layout
            The loaded KLayout layout object.
        layer_tuple : tuple(int, int)
            (layer_number, datatype)
    """
    layer_num, dtype = layer_tuple
    li = ly.find_layer(pya.LayerInfo(layer_num, dtype))
    if li >= 0:
        ly.delete_layer(li)
        print(f"Deleted layer {layer_num}, datatype {dtype}")
    else:
        print(f"Layer {layer_num}, datatype {dtype} not found — nothing to delete.")

def _make_path_element(ly, cell, layer_tuple, start_coord, end_coord, width=0.5):
    """
    Creates a PATH element between two coordinates on the specified layer.

    Parameters:
        ly : pya.Layout
            Active KLayout layout object.
        cell : pya.Cell
            Target cell to insert the path.
        layer_tuple : tuple(int, int)
            (layer_number, datatype) for the layer.
        start_coord : tuple(float, float)
            (x_start, y_start) in microns.
        end_coord : tuple(float, float)
            (x_end, y_end) in microns.
        width : float
            Path width in microns (default = 0.5 µm)
    """
    layer_index = ly.layer(pya.LayerInfo(layer_tuple[0], layer_tuple[1]))
    dbu = ly.dbu  # micron per database unit

    # Convert µm coordinates to database units (integer grid)
    p1 = pya.Point(int(start_coord[0] / dbu), int(start_coord[1] / dbu))
    p2 = pya.Point(int(end_coord[0] / dbu), int(end_coord[1] / dbu))

    path = pya.Path([p1, p2], int(width / dbu))
    cell.shapes(layer_index).insert(path)

def _make_text_label(ly, cell, layer_tuple, text_str, text_size=10.0, position=(0, 0)):
    """
    Adds a TEXT element on the specified layer, with the coordinates included in the label.

    Parameters:
        ly : pya.Layout
            Active KLayout layout object.
        cell : pya.Cell
            Target cell to insert the text.
        layer_tuple : tuple(int, int)
            (layer_number, datatype) for the layer.
        text_str : str
            The text content to insert.
        text_size : float
            Text size in microns (default = 10 µm).
        position : tuple(float, float)
            (x, y) position of the text in microns.
    """
    layer_index = ly.layer(pya.LayerInfo(layer_tuple[0], layer_tuple[1]))
    dbu = ly.dbu

    # Coordinates
    x_um, y_um = position
    x, y = int(x_um / dbu), int(y_um / dbu)

    # Append coordinates to label text
    full_text = f"{text_str}"

    # Create text object
    text = pya.Text(full_text, pya.Trans(pya.Point(x, y)))
    text.size = text_size / dbu

    # Insert into layout
    cell.shapes(layer_index).insert(text)

def _get_component_center_um(design, component_name):  
    """  
    Get the center coordinates of a component in micrometers.  
      
    Args:  
        design: QDesign object  
        component_name (str): Name of the component  
          
    Returns:  
        tuple: (center_x_um, center_y_um) in micrometers  
    """  
    # Get bounds in meters (minx, miny, maxx, maxy)  
    bounds = design.components[component_name].qgeometry_bounds()  
      
    # Calculate center in meters  
    center_x = (bounds[0] + bounds[2]) / 2  
    center_y = (bounds[1] + bounds[3]) / 2  
    if component_name == "transmon":
        center_y = design.components[component_name].p.pos_y
      
    # Convert from meters to micrometers  
    center_x_um = center_x * 1e3 
    center_y_um = center_y * 1e3 
      
    return center_x_um, center_y_um

def _get_pin_center_um(design, component_name, pin_name):  
    """  
    Get the center coordinates of a pin in micrometers.  
      
    Args:  
        design: QDesign object  
        component_name (str): Name of the component  
        pin_name (str): Name of the pin  
          
    Returns:  
        tuple: (center_x_um, center_y_um) in micrometers  
    """  
    # Get the component  
    component = design.components[component_name]  
      
    # Get the pin dictionary  
    pin_dict = component.pins[pin_name]  
      
    # Extract middle coordinates (stored in meters)  
    center_x = pin_dict['middle'][0]  
    center_y = pin_dict['middle'][1]  
      
    # Convert from meters to micrometers  
    center_x_um = center_x * 1e3   
    center_y_um = center_y * 1e3   
      
    return center_x_um, center_y_um

def _get_pin_left_um(design, component_name, pin_name, pin_width):
    """
    Get the left edge coordinates of a pin in micrometers.

    Args:
        design: QDesign object
        component_name (str): Name of the component
        pin_name (str): Name of the pin
        pin_width (float): Width of the pin in meters

    Returns:
        tuple: (left_x_um, left_y_um) in micrometers
    """
    pin_width = pin_width*1e-3
    component = design.components[component_name]
    pin_dict = component.pins[pin_name]
    center = np.array(pin_dict['middle'])
    normal = np.array(pin_dict['normal'])

    tangent = np.array([-normal[1], normal[0]])
    left_pos = center + tangent * pin_width * 0.5

    left_x_um = left_pos[0] * 1e3
    left_y_um = left_pos[1] * 1e3

    return left_x_um, left_y_um

def _get_pin_right_um(design, component_name, pin_name, pin_width):
    """
    Get the right edge coordinates of a pin in micrometers.

    Args:
        design: QDesign object
        component_name (str): Name of the component
        pin_name (str): Name of the pin
        pin_width (float): Width of the pin in meters

    Returns:
        tuple: (right_x_um, right_y_um) in micrometers
    """
    pin_width = pin_width*1e-3
    component = design.components[component_name]
    pin_dict = component.pins[pin_name]
    center = np.array(pin_dict['middle'])
    normal = np.array(pin_dict['normal'])

    tangent = np.array([-normal[1], normal[0]])
    right_pos = center - tangent * pin_width * 0.5

    right_x_um = right_pos[0] * 1e3
    right_y_um = right_pos[1] * 1e3

    return right_x_um, right_y_um

def process_full_gds():
    """
    Process the full.gds file by performing boolean operations and converting layers.
    The processed files are saved as output_file_multi_layer and output_file_single_layer.
    """
    input_file = "gds_files/raw/full.gds"
    output_file_multi_layer = "gds_files/processed_multi_layer/full_processed_multi.gds"
    output_file_single_layer = "gds_files/processed_single_layer/full_processed.gds"

    ly = pya.Layout()
    ly.read(input_file)
    top = ly.top_cell()

    L0_0 = (0,0)
    L1_0 = (1,0)
    L2_0 = (2,0)
    L3_0 = (3,0)
    L4_0 = (4,0)
    L4_10 = (4,10)
    L4_11 = (4,11)

    _boolean_layer_operation(ly, top, L1_0, L3_0, L0_0, operation="AND")
    _boolean_layer_operation(ly, top, L0_0, L4_0, L0_0, operation="AND")
    _boolean_layer_operation(ly, top, L4_10, L4_11, L4_10, operation="OR")
    _delete_layer(ly, L1_0)
    _delete_layer(ly, L3_0)
    _delete_layer(ly, L4_0)
    _delete_layer(ly, L4_11)

    _convert_layer_to_dtype0(ly, 1)
    _convert_layer_to_dtype0(ly, 2)
    _convert_layer_to_dtype0(ly, 3)
    _convert_layer_to_dtype0(ly, 4)

    if os.path.exists(output_file_multi_layer):
        os.remove(output_file_multi_layer)
    ly.write(output_file_multi_layer)
    print(f"Done. Wrote {output_file_multi_layer}")

    _boolean_layer_operation(ly, top, L0_0, L1_0, L0_0, operation="OR")
    _boolean_layer_operation(ly, top, L0_0, L2_0, L0_0, operation="OR")
    _boolean_layer_operation(ly, top, L0_0, L3_0, L0_0, operation="OR")
    _boolean_layer_operation(ly, top, L0_0, L4_0, L0_0, operation="OR")
    _delete_layer(ly, L1_0)
    _delete_layer(ly, L2_0)
    _delete_layer(ly, L3_0)
    _delete_layer(ly, L4_0)

    if os.path.exists(output_file_single_layer):
        os.remove(output_file_single_layer)
    ly.write(output_file_single_layer)
    print(f"Done. Wrote {output_file_single_layer}")

def process_full_no_jj_gds():
    input_file = "gds_files/processed_multi_layer/full_processed_multi.gds"
    output_file_multi_layer = "gds_files/processed_multi_layer/full_no_jj_processed_multi.gds"
    output_file_single_layer = "gds_files/processed_single_layer/full_no_jj_processed.gds"

    ly = pya.Layout()
    ly.read(input_file)
    top = ly.top_cell()
    
    L0_0 = (0,0)
    L1_0 = (1,0)
    L2_0 = (2,0)
    L3_0 = (3,0)
    L4_0 = (4,0)
    _delete_layer(ly, L2_0)
    if os.path.exists(output_file_multi_layer):
        os.remove(output_file_multi_layer)
    ly.write(output_file_multi_layer)
    print(f"Done. Wrote {output_file_multi_layer}")

    _boolean_layer_operation(ly, top, L0_0, L1_0, L0_0, operation="OR")
    _boolean_layer_operation(ly, top, L0_0, L3_0, L0_0, operation="OR")
    _boolean_layer_operation(ly, top, L0_0, L4_0, L0_0, operation="OR")
    _delete_layer(ly, L1_0)
    _delete_layer(ly, L3_0)
    _delete_layer(ly, L4_0)

    if os.path.exists(output_file_single_layer):
        os.remove(output_file_single_layer)
    ly.write(output_file_single_layer)
    print(f"Done. Wrote {output_file_single_layer}")

def process_junction_gds():
    """
    Process the junction.gds file by converting layers to datatype 0.
    The processed files are saved as output_file_multi_layer and output_file_single_layer.
    """
    input_file = "gds_files/raw/junction.gds"
    output_file_multi_layer = "gds_files/processed_multi_layer/junction_processed_multi.gds"
    output_file_single_layer = "gds_files/processed_single_layer/junction_processed.gds"

    ly = pya.Layout()
    ly.read(input_file)
    top = ly.top_cell()

    L0_0 = (0,0)
    L2_0 = (2,0)

    _convert_layer_to_dtype0(ly, 2)

    if os.path.exists(output_file_multi_layer):
        os.remove(output_file_multi_layer)
    ly.write(output_file_multi_layer)
    print(f"Done. Wrote {output_file_multi_layer}")

    _boolean_layer_operation(ly, top, L2_0, L2_0, L0_0, operation="OR")
    _delete_layer(ly, L2_0)

    if os.path.exists(output_file_single_layer):
        os.remove(output_file_single_layer)
    ly.write(output_file_single_layer)
    print(f"Done. Wrote {output_file_single_layer}")

def process_transmon_gds():
    """
    Process the transmon.gds file by converting layers to datatype 0.
    The processed files are saved as output_file_multi_layer and output_file_single_layer.
    """
    input_file = "gds_files/raw/transmon.gds"
    output_file_multi_layer = "gds_files/processed_multi_layer/transmon_processed_multi.gds"
    output_file_single_layer = "gds_files/processed_single_layer/transmon_processed.gds"

    ly = pya.Layout()
    ly.read(input_file)
    top = ly.top_cell()

    L0_0 = (0,0)
    L1_0 = (1,0)
    L2_0 = (2,0)
    L1_10 = (1,10)

    _boolean_layer_operation(ly, top, L1_0, L1_0, L0_0, operation="OR")
    _boolean_layer_operation(ly, top, L1_10, L1_10, L1_0, operation="OR")
    _delete_layer(ly, L1_10)

    _convert_layer_to_dtype0(ly, 1)
    _convert_layer_to_dtype0(ly, 2)

    if os.path.exists(output_file_multi_layer):
        os.remove(output_file_multi_layer)
    ly.write(output_file_multi_layer)
    print(f"Done. Wrote {output_file_multi_layer}")

    _boolean_layer_operation(ly, top, L0_0, L1_0, L0_0, operation="OR")
    _boolean_layer_operation(ly, top, L0_0, L2_0, L0_0, operation="OR")
    _delete_layer(ly, L1_0)
    _delete_layer(ly, L2_0)

    if os.path.exists(output_file_single_layer):
        os.remove(output_file_single_layer)
    ly.write(output_file_single_layer)
    print(f"Done. Wrote {output_file_single_layer}")

def process_transmon_no_jj_gds():
    input_file = "gds_files/processed_multi_layer/transmon_processed_multi.gds"
    output_file_multi_layer = "gds_files/processed_multi_layer/transmon_no_jj_processed_multi.gds"
    output_file_single_layer = "gds_files/processed_single_layer/transmon_no_jj_processed.gds"

    ly = pya.Layout()
    ly.read(input_file)
    top = ly.top_cell()
    
    L0_0 = (0,0)
    L1_0 = (1,0)
    L2_0 = (2,0)

    _delete_layer(ly, L2_0)

    if os.path.exists(output_file_multi_layer):
        os.remove(output_file_multi_layer)
    ly.write(output_file_multi_layer)
    print(f"Done. Wrote {output_file_multi_layer}")

    _boolean_layer_operation(ly, top, L0_0, L1_0, L0_0, operation="OR")
    _delete_layer(ly, L1_0)

    if os.path.exists(output_file_single_layer):
        os.remove(output_file_single_layer)
    ly.write(output_file_single_layer)
    print(f"Done. Wrote {output_file_single_layer}")

def process_resonator_gds():
    """
    Process the resonator.gds file by converting layers to datatype 0.
    The processed files are saved as output_file_multi_layer and output_file_single_layer.
    """
    input_file = "gds_files/raw/resonator.gds"
    output_file_multi_layer = "gds_files/processed_multi_layer/resonator_processed_multi.gds"
    output_file_single_layer = "gds_files/processed_single_layer/resonator_processed.gds"

    ly = pya.Layout()
    ly.read(input_file)
    top = ly.top_cell()

    L0_0 = (0,0)
    L1_0 = (1,0)
    L3_0 = (3,0)
    L3_11 = (3,11)

    _boolean_layer_operation(ly, top, L3_0, L3_0, L0_0, operation="OR")
    _delete_layer(ly, L3_0)
    _boolean_layer_operation(ly, top, L3_11, L3_11, L1_0, operation="OR")
    _delete_layer(ly, L3_11)

    _convert_layer_to_dtype0(ly,1)                              #Can probably remove these

    if os.path.exists(output_file_multi_layer):
        os.remove(output_file_multi_layer)
    ly.write(output_file_multi_layer)
    print(f"Done. Wrote {output_file_multi_layer}")

    _boolean_layer_operation(ly, top, L0_0, L1_0, L0_0, operation="OR")
    _delete_layer(ly, L1_0)

    if os.path.exists(output_file_single_layer):
        os.remove(output_file_single_layer)
    ly.write(output_file_single_layer)
    print(f"Done. Wrote {output_file_single_layer}")

def process_feedline_gds():
    """
    Process the feedline.gds file by converting layers to datatype 0.
    The processed files are saved as output_file_multi_layer and output_file_single_layer.
    """
    input_file = "gds_files/raw/feedline.gds"
    output_file_multi_layer = "gds_files/processed_multi_layer/feedline_processed_multi.gds"
    output_file_single_layer = "gds_files/processed_single_layer/feedline_processed.gds"

    ly = pya.Layout()
    ly.read(input_file)
    top = ly.top_cell()

    L0_0 = (0,0)
    L1_0 = (1,0)
    L4_0 = (4,0)
    L4_10 = (4,10)
    L4_11 = (4,11)

    _boolean_layer_operation(ly, top, L4_10, L4_11, L4_10, operation="OR")
    _delete_layer(ly, L4_11)
    _boolean_layer_operation(ly, top, L4_0, L4_0, L0_0, operation="OR")
    _delete_layer(ly, L4_0)
    _boolean_layer_operation(ly, top, L4_10, L4_10, L1_0, operation="OR")
    _delete_layer(ly, L4_10)

    if os.path.exists(output_file_multi_layer):
        os.remove(output_file_multi_layer)
    ly.write(output_file_multi_layer)
    print(f"Done. Wrote {output_file_multi_layer}")

    _convert_layer_to_dtype0(ly, 3)

    _boolean_layer_operation(ly, top, L0_0, L1_0, L0_0, operation="OR")
    _delete_layer(ly, L1_0)

    if os.path.exists(output_file_single_layer):
        os.remove(output_file_single_layer)
    ly.write(output_file_single_layer)
    print(f"Done. Wrote {output_file_single_layer}")

def make_inductex_cap_sim_full_gds(full_design, qubit_geo):
    return
    input_file = "gds_files/processed_single_layer/full_processed.gds"
    output_file = "simulations/inductex/capacitance/full/full.gds"

    ly = pya.Layout()
    ly.read(input_file)
    top = ly.top_cell()

    # Example: Add labels or paths as needed here

    if os.path.exists(output_file):
        os.remove(output_file)
    ly.write(output_file)
    print(f"Done. Wrote {output_file}")

def make_inductex_cap_sim_full_no_jj_gds(full_no_jj_design, qubit_geo):
    return
    input_file = "gds_files/processed_single_layer/full_no_jj_processed.gds"
    output_file = "simulations/inductex/capacitance/full_no_jj/full_no_jj.gds"

    ly = pya.Layout()
    ly.read(input_file)
    top = ly.top_cell()

    # Example: Add labels or paths as needed here

    if os.path.exists(output_file):
        os.remove(output_file)
    ly.write(output_file)
    print(f"Done. Wrote {output_file}")

def make_inductex_cap_sim_junction_gds(junction_design, qubit_geo):
    return
    input_file = "gds_files/processed_single_layer/junction_processed.gds"
    output_file = "simulations/inductex/capacitance/junction/junction.gds"

    ly = pya.Layout()
    ly.read(input_file)
    top = ly.top_cell()

    # Example: Add labels or paths as needed here

    if os.path.exists(output_file):
        os.remove(output_file)
    ly.write(output_file)
    print(f"Done. Wrote {output_file}")

def make_inductex_cap_sim_transmon_gds(transmon_design, qubit_geo):
    input_file = "gds_files/processed_single_layer/transmon_processed.gds"
    output_file = "simulations/inductex/capacitance/transmon/transmon_cap.gds"
    
    ly = pya.Layout()
    ly.read(input_file)
    top = ly.top_cell()

    L100_0 = (100, 0)

    (cx, cy)            = _get_component_center_um(transmon_design, "transmon")
    (px, py)            = _get_pin_center_um(transmon_design, "transmon", "readout")
    cross_length        = extract_um(qubit_geo["cross_length"])
    cross_gap           = extract_um(qubit_geo["cross_gap"])
    claw_cpw_length     = extract_um(qubit_geo["claw_cpw_length"])
    p1                  = (cx + cross_length + 2*cross_gap, cy)
    p2                  = (px, py - claw_cpw_length/2)
    p3                  = (cx, cy)

    # Order which appears here is how they are positioned in capacitance matrix
    _make_text_label(ly, top, L100_0, "Cground main", text_size=10.0, position=p1)
    _make_text_label(ly, top, L100_0, "Ccross main", text_size=10.0, position=p3)
    _make_text_label(ly, top, L100_0, "Cclaw main", text_size=10.0, position=p2)

    if os.path.exists(output_file):
        os.remove(output_file)
    ly.write(output_file)
    print(f"Done. Wrote {output_file}")

def make_inductex_cap_sim_transmon_no_jj_gds(transmon_no_jj_design, qubit_geo):
    input_file = "gds_files/processed_single_layer/transmon_no_jj_processed.gds"
    output_file = "simulations/inductex/capacitance/transmon_no_jj/transmon_no_jj_cap.gds"        # Might change all functions to take file paths, then can remove this
    
    ly = pya.Layout()
    ly.read(input_file)
    top = ly.top_cell()

    L100_0 = (100, 0)
    
    (cx, cy)            = _get_component_center_um(transmon_no_jj_design, "transmon")
    (px, py)            = _get_pin_center_um(transmon_no_jj_design, "transmon", "readout")
    cross_length        = extract_um(qubit_geo["cross_length"])
    cross_gap           = extract_um(qubit_geo["cross_gap"])
    claw_cpw_length     = extract_um(qubit_geo["claw_cpw_length"])
    p1                  = (cx + cross_length + 2*cross_gap, cy)
    p2                  = (px, py - claw_cpw_length/2)
    p3                  = (cx, cy)
    # Order which appears here is how they are positioned in capacitance matrix
    _make_text_label(ly, top, L100_0, "Cground main", text_size=10.0, position=p1)
    _make_text_label(ly, top, L100_0, "Ccross main", text_size=10.0, position=p3)
    _make_text_label(ly, top, L100_0, "Cclaw main", text_size=10.0, position=p2)

    if os.path.exists(output_file):
        os.remove(output_file)
    ly.write(output_file)
    print(f"Done. Wrote {output_file}")

def make_inductex_cap_sim_resonator_gds(resonator_design, resonator_geo):
    input_file = "gds_files/processed_single_layer/resonator_processed.gds"
    output_file = "simulations/inductex/capacitance/resonator/resonator_cap.gds"

    ly = pya.Layout()
    ly.read(input_file)
    top = ly.top_cell()

    L100_0 = (100, 0)
    
    (px, py)            = _get_pin_center_um(resonator_design, "resonator", "end")
    trace_gap           = extract_um(resonator_geo["trace_gap"])
    trace_width         = extract_um(resonator_geo["trace_width"])
    p1                  = (px + trace_width + 2*trace_gap, py)
    p2                  = (px, py + trace_width)

    _make_text_label(ly, top, L100_0, "Cground main", text_size=10.0, position=p1)
    _make_text_label(ly, top, L100_0, "Cresonator main", text_size=10.0, position=p2)

    if os.path.exists(output_file):
        os.remove(output_file)
    ly.write(output_file)
    print(f"Done. Wrote {output_file}")

def make_inductex_cap_sim_feedline_gds(feedline_design, feedline_geo):
    input_file = "gds_files/processed_single_layer/feedline_processed.gds"
    output_file = "simulations/inductex/capacitance/feedline/feedline_cap.gds"

    ly = pya.Layout()
    ly.read(input_file)
    top = ly.top_cell()

    L100_0 = (100, 0)

    (cx, cy)            = _get_component_center_um(feedline_design, "feedline")
    prime_width        = extract_um(feedline_geo["prime_width"])

    p1                  = (cx, cy + 3*prime_width)
    p2                  = (cx, cy)

    _make_text_label(ly, top, L100_0, "Cground main", text_size=10.0, position=p1)
    _make_text_label(ly, top, L100_0, "Cfeedline main", text_size=10.0, position=p2)

    if os.path.exists(output_file):
        os.remove(output_file)
    ly.write(output_file)
    print(f"Done. Wrote {output_file}")

def make_inductex_ind_sim_full_gds(full_design, qubit_geo):
    return
    input_file = "gds_files/processed_single_layer/full_processed.gds"
    output_file = "simulations/inductex/inductance/full/full.gds"

    ly = pya.Layout()
    ly.read(input_file)
    top = ly.top_cell()

    # Example: Add labels or paths as needed here

    if os.path.exists(output_file):
        os.remove(output_file)
    ly.write(output_file)
    print(f"Done. Wrote {output_file}")

def make_inductex_ind_sim_full_no_jj_gds(full_no_jj_design, qubit_geo):
    return
    input_file = "gds_files/processed_single_layer/full_no_jj_processed.gds"
    output_file = "simulations/inductex/inductance/full_no_jj/full_no_jj.gds"

    ly = pya.Layout()
    ly.read(input_file)
    top = ly.top_cell()

    # Example: Add labels or paths as needed here

    if os.path.exists(output_file):
        os.remove(output_file)
    ly.write(output_file)
    print(f"Done. Wrote {output_file}")

def make_inductex_ind_sim_junction_gds(junction_design, qubit_geo):
    return
    input_file = "gds_files/processed_single_layer/junction_processed.gds"
    output_file = "simulations/inductex/inductance/junction/junction.gds"

    ly = pya.Layout()
    ly.read(input_file)
    top = ly.top_cell()

    # Example: Add labels or paths as needed here

    if os.path.exists(output_file):
        os.remove(output_file)
    ly.write(output_file)
    print(f"Done. Wrote {output_file}")

def make_inductex_ind_sim_transmon_gds(transmon_design, qubit_geo):
    return
    input_file = "gds_files/processed_single_layer/transmon_processed.gds"
    output_file = "simulations/inductex/inductance/transmon/transmon.gds"

    ly = pya.Layout()
    ly.read(input_file)
    top = ly.top_cell()

    # Example: Add labels or paths as needed here

    if os.path.exists(output_file):
        os.remove(output_file)
    ly.write(output_file)
    print(f"Done. Wrote {output_file}")

def make_inductex_ind_sim_transmon_no_jj_gds(transmon_no_jj_design, qubit_geo):
    return
    input_file = "gds_files/processed_single_layer/transmon_no_jj_processed.gds"
    output_file = "simulations/inductex/inductance/transmon_no_jj/transmon_no_jj.gds"

    ly = pya.Layout()
    ly.read(input_file)
    top = ly.top_cell()

    # Example: Add labels or paths as needed here

    if os.path.exists(output_file):
        os.remove(output_file)
    ly.write(output_file)
    print(f"Done. Wrote {output_file}")

def make_inductex_ind_sim_resonator_gds(resonator_design, resonator_geo):
    input_file = "gds_files/processed_multi_layer/resonator_processed_multi.gds"
    output_file = "simulations/inductex/inductance/resonator/resonator_ind.gds"

    ly = pya.Layout()
    ly.read(input_file)
    top = ly.top_cell()

    L100_0 = (100, 0)
    L110_0 = (110, 0)

    trace_gap        = extract_um(resonator_geo["trace_gap"])
    trace_width      = extract_um(resonator_geo["trace_width"])
    (start_left_x, start_left_y)   = _get_pin_left_um(resonator_design, "resonator", "start", trace_width)
    (start_right_x, start_right_y) = _get_pin_right_um(resonator_design, "resonator", "start", trace_width)
    (start_mid_x, start_mid_y)     = _get_pin_center_um(resonator_design, "resonator", "start")
    (end_left_x, end_left_y)       = _get_pin_left_um(resonator_design, "resonator", "end", trace_width)
    (end_right_x, end_right_y)     = _get_pin_right_um(resonator_design, "resonator", "end", trace_width)
    (end_mid_x, end_mid_y)         = _get_pin_center_um(resonator_design, "resonator", "end")

    p1 = (start_mid_x, start_mid_y)
    p2 = (start_mid_x - trace_gap, start_mid_y)
    p3 = (end_mid_x, end_mid_y)
    p4 = (end_mid_x, end_mid_y - trace_gap)

    _make_text_label(ly, top, L100_0, "P1+ main", text_size=10.0, position=p1)
    _make_text_label(ly, top, L100_0, "P1- gd", text_size=10.0, position=p2)
    _make_text_label(ly, top, L100_0, "P2+ main", text_size=10.0, position=p3)
    _make_text_label(ly, top, L100_0, "P2- gd", text_size=10.0, position=p4)

    start_left      = (start_left_x, start_left_y)
    start_right     = (start_right_x, start_right_y)
    start_left_gap  = (start_left_x - trace_gap, start_left_y)
    start_right_gap = (start_right_x - trace_gap, start_right_y)
    end_left        = (end_left_x, end_left_y)
    end_right       = (end_right_x, end_right_y)
    end_left_gap    = (end_left_x, end_left_y - trace_gap)
    end_right_gap   = (end_right_x, end_right_y - trace_gap)

    _make_path_element(ly, top, L110_0, start_left, start_right, width=0.5)
    _make_path_element(ly, top, L110_0, start_left_gap, start_right_gap, width=0.5)
    _make_path_element(ly, top, L110_0, end_left, end_right, width=0.5)
    _make_path_element(ly, top, L110_0, end_left_gap, end_right_gap, width=0.5)

    if os.path.exists(output_file):
        os.remove(output_file)
    ly.write(output_file)
    print(f"Done. Wrote {output_file}")

def make_inductex_ind_sim_feedline_gds(feedline_design, feedline_geo):
    pass
    input_file = "gds_files/processed_multi_layer/feedline_processed_multi.gds"
    output_file = "simulations/inductex/inductance/feedline/feedline_ind.gds"

    ly = pya.Layout()
    ly.read(input_file)
    top = ly.top_cell()

    if os.path.exists(output_file):
        os.remove(output_file)
    ly.write(output_file)
    print(f"Done. Wrote {output_file}")

def make_inductex_sparam_sim_full_gds(full_design, qubit_geo):
    return
    input_file = "gds_files/processed_single_layer/full_processed.gds"
    output_file = "simulations/inductex/sparams/full/full.gds"

    ly = pya.Layout()
    ly.read(input_file)
    top = ly.top_cell()

    # Example: Add labels or paths as needed here

    if os.path.exists(output_file):
        os.remove(output_file)
    ly.write(output_file)
    print(f"Done. Wrote {output_file}")

def make_inductex_sparam_sim_full_no_jj_gds(full_no_jj_design, qubit_geo):
    return
    input_file = "gds_files/processed_single_layer/full_no_jj_processed.gds"
    output_file = "simulations/inductex/sparams/full_no_jj/full_no_jj.gds"

    ly = pya.Layout()
    ly.read(input_file)
    top = ly.top_cell()

    # Example: Add labels or paths as needed here

    if os.path.exists(output_file):
        os.remove(output_file)
    ly.write(output_file)
    print(f"Done. Wrote {output_file}")

def make_inductex_sparam_sim_junction_gds(junction_design, qubit_geo):
    return
    input_file = "gds_files/processed_single_layer/junction_processed.gds"
    output_file = "simulations/inductex/sparams/junction/junction.gds"

    ly = pya.Layout()
    ly.read(input_file)
    top = ly.top_cell()

    # Example: Add labels or paths as needed here

    if os.path.exists(output_file):
        os.remove(output_file)
    ly.write(output_file)
    print(f"Done. Wrote {output_file}")

def make_inductex_sparam_sim_transmon_gds(transmon_design, qubit_geo):

    input_file = "gds_files/processed_single_layer/transmon_processed.gds"
    output_file = "simulations/inductex/sparams/transmon/transmon_sparams.gds"

    ly = pya.Layout()
    ly.read(input_file)
    top = ly.top_cell()

    L100_0 = (100, 0)
    L110_0 = (110, 0)

    cross_width        = extract_um(qubit_geo["cross_width"])
    cross_length        = extract_um(qubit_geo["cross_length"])
    cross_gap           = extract_um(qubit_geo["cross_gap"])
    claw_gap            = extract_um(qubit_geo["claw_gap"])
    claw_width          = extract_um(qubit_geo["claw_width"])

    (cx, cy)            = _get_component_center_um(transmon_design, "transmon")
    (start_left_x, start_left_y)   = _get_pin_left_um(transmon_design, "transmon", "readout", claw_width)
    (start_right_x, start_right_y) = _get_pin_right_um(transmon_design, "transmon", "readout", claw_width)
    (start_mid_x, start_mid_y)     = _get_pin_center_um(transmon_design, "transmon", "readout")

    p1 = (start_mid_x, start_mid_y)
    p2 = (start_mid_x , start_mid_y + claw_gap)
    p3 = (cx + cross_length, cy)
    p4 = (cx + cross_length + cross_gap, cy)

    _make_text_label(ly, top, L100_0, "P1+ main", text_size=10.0, position=p1)
    _make_text_label(ly, top, L100_0, "P1- main", text_size=10.0, position=p2)
    _make_text_label(ly, top, L100_0, "P2+ main", text_size=10.0, position=p3)
    _make_text_label(ly, top, L100_0, "P2- main", text_size=10.0, position=p4)

    start_left      = (start_left_x, start_left_y)
    start_right     = (start_right_x, start_right_y)
    start_left_gap  = (start_left_x , start_left_y + claw_gap)
    start_right_gap = (start_right_x , start_right_y + claw_gap)
    end_left        = (cx + cross_length, cy - cross_width/2)
    end_right       = (cx + cross_length, cy + cross_width/2)
    end_left_gap    = (cx + cross_length + cross_gap, cy + cross_width/2)
    end_right_gap   = (cx + cross_length + cross_gap, cy - cross_width/2)

    _make_path_element(ly, top, L110_0, start_left, start_right, width=0.5)
    _make_path_element(ly, top, L110_0, start_left_gap, start_right_gap, width=0.5)
    _make_path_element(ly, top, L110_0, end_left, end_right, width=0.5)
    _make_path_element(ly, top, L110_0, end_left_gap, end_right_gap, width=0.5)


    if os.path.exists(output_file):
        os.remove(output_file)
    ly.write(output_file)
    print(f"Done. Wrote {output_file}")

def make_inductex_sparam_sim_transmon_no_jj_gds(transmon_no_jj_design, qubit_geo):
    input_file = "gds_files/processed_single_layer/transmon_no_jj_processed.gds"
    output_file = "simulations/inductex/sparams/transmon_no_jj/transmon_no_jj_sparams.gds"

    ly = pya.Layout()
    ly.read(input_file)
    top = ly.top_cell()

    L100_0 = (100, 0)
    L110_0 = (110, 0)

    cross_width        = extract_um(qubit_geo["cross_width"])
    cross_length        = extract_um(qubit_geo["cross_length"])
    cross_gap           = extract_um(qubit_geo["cross_gap"])
    claw_gap            = extract_um(qubit_geo["claw_gap"])
    claw_width          = extract_um(qubit_geo["claw_width"])

    (cx, cy)            = _get_component_center_um(transmon_no_jj_design, "transmon")
    (start_left_x, start_left_y)   = _get_pin_left_um(transmon_no_jj_design, "transmon", "readout", claw_width)
    (start_right_x, start_right_y) = _get_pin_right_um(transmon_no_jj_design, "transmon", "readout", claw_width)
    (start_mid_x, start_mid_y)     = _get_pin_center_um(transmon_no_jj_design, "transmon", "readout")

    p1 = (start_mid_x, start_mid_y)
    p2 = (start_mid_x , start_mid_y + claw_gap)
    p3 = (cx + cross_length, cy)
    p4 = (cx + cross_length + cross_gap, cy)

    _make_text_label(ly, top, L100_0, "P1+ main", text_size=10.0, position=p1)
    _make_text_label(ly, top, L100_0, "P1- main", text_size=10.0, position=p2)
    _make_text_label(ly, top, L100_0, "P2+ main", text_size=10.0, position=p3)
    _make_text_label(ly, top, L100_0, "P2- main", text_size=10.0, position=p4)

    start_left      = (start_left_x, start_left_y)
    start_right     = (start_right_x, start_right_y)
    start_left_gap  = (start_left_x , start_left_y + claw_gap)
    start_right_gap = (start_right_x , start_right_y + claw_gap)
    end_left        = (cx + cross_length, cy - cross_width/2)
    end_right       = (cx + cross_length, cy + cross_width/2)
    end_left_gap    = (cx + cross_length + cross_gap, cy + cross_width/2)
    end_right_gap   = (cx + cross_length + cross_gap, cy - cross_width/2)

    _make_path_element(ly, top, L110_0, start_left, start_right, width=0.5)
    _make_path_element(ly, top, L110_0, start_left_gap, start_right_gap, width=0.5)
    _make_path_element(ly, top, L110_0, end_left, end_right, width=0.5)
    _make_path_element(ly, top, L110_0, end_left_gap, end_right_gap, width=0.5)


    if os.path.exists(output_file):
        os.remove(output_file)
    ly.write(output_file)
    print(f"Done. Wrote {output_file}")

def make_inductex_sparam_sim_resonator_gds(resonator_design, resonator_geo):
    input_file = "gds_files/processed_single_layer/resonator_processed.gds"
    output_file = "simulations/inductex/sparams/resonator/resonator_sparams.gds"

    ly = pya.Layout()
    ly.read(input_file)
    top = ly.top_cell()

    L100_0 = (100, 0)
    L110_0 = (110, 0)

    trace_gap        = extract_um(resonator_geo["trace_gap"])
    trace_width      = extract_um(resonator_geo["trace_width"])
    (start_left_x, start_left_y)   = _get_pin_left_um(resonator_design, "resonator", "start", trace_width)
    (start_right_x, start_right_y) = _get_pin_right_um(resonator_design, "resonator", "start", trace_width)
    (start_mid_x, start_mid_y)     = _get_pin_center_um(resonator_design, "resonator", "start")
    (end_left_x, end_left_y)       = _get_pin_left_um(resonator_design, "resonator", "end", trace_width)
    (end_right_x, end_right_y)     = _get_pin_right_um(resonator_design, "resonator", "end", trace_width)
    (end_mid_x, end_mid_y)         = _get_pin_center_um(resonator_design, "resonator", "end")

    p1 = (start_mid_x, start_mid_y)
    p2 = (start_mid_x - trace_gap, start_mid_y)
    p3 = (end_mid_x, end_mid_y)
    p4 = (end_mid_x, end_mid_y - trace_gap)

    _make_text_label(ly, top, L100_0, "P1+ main", text_size=10.0, position=p1)
    _make_text_label(ly, top, L100_0, "P1- main", text_size=10.0, position=p2)
    _make_text_label(ly, top, L100_0, "P2+ main", text_size=10.0, position=p3)
    _make_text_label(ly, top, L100_0, "P2- main", text_size=10.0, position=p4)

    start_left      = (start_left_x, start_left_y)
    start_right     = (start_right_x, start_right_y)
    start_left_gap  = (start_left_x - trace_gap, start_left_y)
    start_right_gap = (start_right_x - trace_gap, start_right_y)
    end_left        = (end_left_x, end_left_y)
    end_right       = (end_right_x, end_right_y)
    end_left_gap    = (end_left_x, end_left_y - trace_gap)
    end_right_gap   = (end_right_x, end_right_y - trace_gap)

    _make_path_element(ly, top, L110_0, start_left, start_right, width=0.5)
    _make_path_element(ly, top, L110_0, start_left_gap, start_right_gap, width=0.5)
    _make_path_element(ly, top, L110_0, end_left, end_right, width=0.5)
    _make_path_element(ly, top, L110_0, end_left_gap, end_right_gap, width=0.5)


    if os.path.exists(output_file):
        os.remove(output_file)
    ly.write(output_file)
    print(f"Done. Wrote {output_file}")

def make_inductex_sparam_sim_feedline_gds(feedline_design, feedline_geo):
    return
    input_file = "gds_files/processed_single_layer/feedline_processed.gds"
    output_file = "simulations/inductex/sparams/feedline/feedline.gds"

    ly = pya.Layout()
    ly.read(input_file)
    top = ly.top_cell()

    # Example: Add labels or paths as needed here

    if os.path.exists(output_file):
        os.remove(output_file)
    ly.write(output_file)
    print(f"Done. Wrote {output_file}")
