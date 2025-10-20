from qiskit_metal import designs, Dict
from qiskit_metal.qlibrary.tlines.meandered import RouteMeander
from qiskit_metal.qlibrary.tlines.pathfinder import RoutePathfinder
from qiskit_metal.qlibrary.terminations.launchpad_wb import LaunchpadWirebond
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround
from SQDMetal.Comps.Xmon import Xmon
from SQDMetal.Comps.Junctions import JunctionDolanPinStretch
from qiskit_metal.qlibrary.qubits.transmon_cross import TransmonCross

from qiskit_metal.qlibrary.couplers.coupled_line_tee import CoupledLineTee
from qiskit_metal.qlibrary.terminations.short_to_ground import ShortToGround

from helpers import extract_um

def _center_and_scale_chip(design, component_name='transmon', scale_factor=1.3):
    # Get the components's bounding box (minx, miny, maxx, maxy)  
    comp_bounds = design.components[component_name].qgeometry_bounds()  
    
    # Calculate component's dimensions  
    comp_width = comp_bounds[2] - comp_bounds[0]  # maxx - minx  
    comp_height = comp_bounds[3] - comp_bounds[1]  # maxy - miny  
    
    # Calculate component's center  
    comp_center_x = (comp_bounds[0] + comp_bounds[2]) / 2  
    comp_center_y = (comp_bounds[1] + comp_bounds[3]) / 2  
    
    # Special handling for feedline to ensure enough height
    if component_name == 'feedline':
        comp_height = design.components[component_name].p.trace_width * 20 * scale_factor

    # Apply scaling factor  
    chip_width = comp_width * scale_factor  
    chip_height = comp_height * scale_factor  

    # Set chip size and center  
    design.chips.main.size.size_x = f'{chip_width}mm'  
    design.chips.main.size.size_y = f'{chip_height}mm'  
    design.chips.main.size.size_z = '-280um'                                                                #MAYBE CHANGE THIS
    design.chips.main.size.center_x = f'{comp_center_x}mm'  
    design.chips.main.size.center_y = f'{comp_center_y}mm' 

def make_full_design(best_design, qubit_geo, resonator_geo, feedline_geo):
    full_design = designs.DesignPlanar({}, overwrite_enabled=True);

    # Set up chip dimensions 
    full_design.chips.main.size.size_x = '4.6mm'
    full_design.chips.main.size.size_y = '2.4mm'
    full_design.chips.main.size.size_z = '-750um'                                                               #changed this
    full_design.chips.main.size.center_x = '0mm'
    full_design.chips.main.size.center_y = '-1mm'

    layer_transmon  = 1
    layer_junction  = 2
    layer_xmon      = 1
    layer_rout      = 3
    layer_otg       = 3
    layer_launch1   = 4
    layer_launch2   = 4
    layer_feed      = 4

    # Adjust qubit positions and orientation to make compact
    qubit_options                   = best_design["design_options_qubit"]
    qubit_options["orientation"]    = 0
    qubit_options["pos_x"]          = '0.375mm'
    qubit_options["pos_y"]          = '-1.3'
    qubit_options["layer"]          = f'{layer_transmon}'

    # Create TransmonCross (pin name: 'readout')
    TransmonCross(full_design, "transmon", options=qubit_options)

    # Components are tested seperately so some geometry is tuned for this
    # we chose this to adjust it for a complete system: DOESNT WORK, CAN USE IN REPORT FOR 
    # SOURCES OF ERROR, too much time needed to fix such a small issue
    # Q1.options['connection_pads']['readout']['cpw_width'] = resonator_geo["trace_width"]
    # Q1.options['connection_pads']['readout']['cpw_gap'] = '0um'  

    # Extract geometry parameters
    cross_width     = extract_um(qubit_geo['cross_width'])
    cross_length    = extract_um(qubit_geo['cross_length'])
    cross_gap       = extract_um(qubit_geo['cross_gap'])

    # Create dummy structure to place junction, will delete later
    xmon = Xmon(
        full_design,
        'xmon',
        options=Dict(
            pos_x='0.375mm',
            pos_y='-1.3',
            layer=layer_xmon,
            hBar_width=f"{cross_width}um",
            vBar_width=f"{cross_width}um",
            vBar_gap=f"{cross_gap}um",
            hBar_gap=f"{cross_gap}um",
            cross_width=f"{2 * cross_length}um",
            cross_height=f"{2 * cross_length}um",
            gap_up=qubit_geo['cross_gap'],
            gap_left=qubit_geo['cross_gap'],
            gap_right=qubit_geo['cross_gap'],
            gap_down=qubit_geo['cross_gap']
        )
    )

    # Junction
    JunctionDolanPinStretch(
        full_design,
        'junction',
        options=Dict(
            pin_inputs=Dict(
                start_pin=Dict(
                    component='xmon',
                    pin='right'
                )
            ),
            layer=layer_junction,
            dist_extend=qubit_geo['cross_gap'],
            finger_width='0.4um',
            t_pad_size='0.385um',
            squid_width='5.4um',
            prong_width='0.9um'
        )
    )
    # Launchpads
    # launch pads are placed arbitrarily as x,y positions where not given on blueprint
    x1 = '-2mm'  # Replace this with database length?
    y1 = '0mm'
    x2 = '2mm'
    launch_options1 = dict(
        chip='main',
        layer=layer_launch1,
        pos_x=x1,
        pos_y=y1,
        orientation='360',
        lead_length='30um',
        pad_height='103um',
        pad_width='103um',
        pad_gap='60um',
        trace_width=resonator_geo["trace_width"],
        trace_gap=resonator_geo["trace_gap"]
    )

    launch_options2 = dict(
        chip='main',
        layer=layer_launch2,
        pos_x=x2,
        pos_y=y1,
        orientation='180',
        lead_length='30um',
        pad_height='103um',
        pad_width='103um',
        pad_gap='60um',
        trace_width=resonator_geo["trace_width"],
        trace_gap=resonator_geo["trace_gap"]
    )

    LaunchpadWirebond(full_design, 'LP1', options=launch_options1)
    LaunchpadWirebond(full_design, 'LP2', options=launch_options2)

    # Use RoutePathfinder to connect the two launchpads
    # The refenence example did not follow the strict database geometries and 
    # opted to use a RoutePathfinder instead of the data component coupled_line_tee
    # works just fine so will keep it this way for nw
    RoutePathfinder(
        full_design, 'feedline', 
        options=dict(
            chip='main',
            layer=layer_feed,
            trace_width=resonator_geo["trace_width"],
            trace_gap=resonator_geo["trace_gap"],
            fillet='90um',
            hfss_wire_bonds=True,
            lead=dict(end_straight='0.1mm'),
            pin_inputs=Dict(
                start_pin=Dict(
                    component='LP1',
                    pin='tie'
                ),
                end_pin=Dict(
                    component='LP2',
                    pin='tie'
                )
            )
        )
    )

    # Will implement the CoupledLineTee in another version of the program
    # for now the ghetto feedline will be used
    """
    coupler = CoupledLineTee(  
        full_design,   
        'coupler_name',   
        options=dict(  
            pos_x='1.0mm',  # Set your position  
            pos_y='2.2mm',  # Set your position  
            orientation=feedline_geo.get("orientation", "0"),  
            coupling_length=feedline_geo.get("coupling_length", "350um"),  
            coupling_space=feedline_geo.get("coupling_space", "7.9um"),  
            down_length=feedline_geo.get("down_length", "50um"),  
            prime_width=feedline_geo.get("prime_width", "11.7um"),  
            prime_gap=feedline_geo.get("prime_gap", "5.1um"),  
            second_width=feedline_geo.get("second_width", "11.7um"),  
            second_gap=feedline_geo.get("second_gap", "5.1um"),  
            open_termination=feedline_geo.get("open_termination", False),  
            # Optional capacitor parameters (if needed)  
            cap_distance=feedline_geo.get("cap_distance", None),  
            cap_gap=feedline_geo.get("cap_gap", None),  
            cap_gap_ground=feedline_geo.get("cap_gap_ground", None),  
            cap_width=feedline_geo.get("cap_width", None),  
            finger_count=feedline_geo.get("finger_count", None),  
            finger_length=feedline_geo.get("finger_length", None)  
        )  
    )
    """

    # Open to ground placed at arbitary x,y position as not given in blueprint
    OpenToGround(
            full_design,
            'otg1',
            options=dict(
                    chip='main',
                    layer=layer_otg,
                    pos_x='-0.2mm',
                    pos_y='-40um',
                    orientation=180,
                    termination_gap = resonator_geo["trace_gap"],
            )
    )

    # Resonator and feedline gap width (W) and center conductor width (S)
    full_design.variables['cpw_width'] = resonator_geo["trace_width"]
    full_design.variables['cpw_gap'] = resonator_geo["trace_gap"]

    # The resonator is a similar story to the feedline but have updated it to be more 
    # compliant with the database geometry

    # Original resonator from tutorial
    """
    res1 = RouteMeander( 
        full_design,
        'resonator',
        Dict(
            trace_width='10um',
            trace_gap='6um',
            total_length='3.7mm',
            hfss_wire_bonds=False,
            fillet='99.9 um',
            lead=dict(
                start_straight='300um'
            ),
            pin_inputs=Dict(
                start_pin=Dict(
                    component='otg1',
                    pin='open'
                ),
                end_pin=Dict(
                    component='claw',
                    pin='a'
                )
            )
        )
    )
    """
    # Proper implementation of the resonator using database geometry
    """    
    res1 = RouteMeander(
        full_design, 'resonator', Dict(
            layer=layer_rout,
            trace_width=resonator_geo['trace_width'],
            trace_gap=resonator_geo['trace_gap'],
            total_length=resonator_geo['total_length'],
            fillet=resonator_geo['fillet'],
            hfss_wire_bonds=False,
            lead=Dict(
                start_straight=resonator_geo['start_straight'],
            ),
            meander=Dict(
                asymmetry=resonator_geo['asymmetry'],
                spacing=resonator_geo['spacing']
            ),
            pin_inputs=Dict(
                start_pin=Dict(component='otg1', pin='open'),
                end_pin=Dict(component='transmon', pin='readout')
            )
        )
    )
    """
    # Adjusted resonator from tutorial, still 'incorrect' asymmetry, spacing and start length but shouldn't effect result
    RouteMeander(
        full_design,
        'resonator',
        Dict(
            layer=layer_rout,
            trace_width=resonator_geo['trace_width'],
            trace_gap=resonator_geo['trace_gap'],
            total_length=resonator_geo['total_length'],
            hfss_wire_bonds=False,
            fillet=resonator_geo['fillet'],
            lead=dict(
                start_straight='100um',                                                                 # This is not strictly correct but won't effect results
            ),
            meander=Dict(
                spacing=resonator_geo['spacing'],   # Spacing between meander segments
                asymmetry=resonator_geo['asymmetry']    # Asymmetry offset for meander pattern
            ),
            pin_inputs=Dict(
                start_pin=Dict(
                    component='otg1',
                    pin='open'
                ),
                end_pin=Dict(
                    component='transmon',
                    pin='readout'
                )
            )
        )
    )

    # Use RouteMeander to fix the total length of the resonator
    """res1 = RouteMeander(
        full_design,
        'resonator',
        Dict(
            layer=layer_rout,
            trace_width='10um',
            trace_gap='6um',
            total_length='3.7mm',
            hfss_wire_bonds=False,
            fillet='99.9um',
            lead=Dict(start_straight='300um'),
            pin_inputs=Dict(
                start_pin=Dict(component='otg1', pin='open'),
                end_pin=Dict(component='transmon', pin='readout')
            )
        )
    )"""
    full_design.delete_component('xmon')
    return full_design

def make_full_no_jj_design(best_design, qubit_geo, resonator_geo, feedline_geo):
    full_no_jj_design = make_full_design(best_design, qubit_geo, resonator_geo, feedline_geo)
    full_no_jj_design.delete_component('junction')
    return full_no_jj_design

def make_junction_design(best_design, qubit_geo, resonator_geo, feedline_geo, scaling_factor=1.3):
    junction_design = make_full_design(best_design, qubit_geo, resonator_geo, feedline_geo)
    junction_design.delete_component('LP1')
    junction_design.delete_component('LP2')
    junction_design.delete_component('otg1')
    junction_design.delete_component('feedline')
    junction_design.delete_component('transmon')
    junction_design.delete_component('resonator')
    _center_and_scale_chip(junction_design, component_name='junction', scale_factor=scaling_factor)
    return junction_design

def make_transmon_design(best_design, qubit_geo, resonator_geo, feedline_geo, scaling_factor=1.3):
    transmon_design = make_full_design(best_design, qubit_geo, resonator_geo, feedline_geo)
    transmon_design.delete_component('LP1')
    transmon_design.delete_component('LP2')
    transmon_design.delete_component('otg1')
    transmon_design.delete_component('feedline')
    transmon_design.delete_component('resonator')
    _center_and_scale_chip(transmon_design, component_name='transmon', scale_factor=scaling_factor)
    return transmon_design

def make_transmon_no_jj_design(best_design, qubit_geo, resonator_geo, feedline_geo, scaling_factor=1.3):
    transmon_design = make_full_design(best_design, qubit_geo, resonator_geo, feedline_geo)
    transmon_design.delete_component('LP1')
    transmon_design.delete_component('LP2')
    transmon_design.delete_component('otg1')
    transmon_design.delete_component('feedline')
    transmon_design.delete_component('resonator')
    transmon_design.delete_component('junction')
    _center_and_scale_chip(transmon_design, component_name='transmon', scale_factor=scaling_factor)
    return transmon_design

def make_resonator_design(best_design, qubit_geo, resonator_geo, feedline_geo, scaling_factor=1.3):
    resonator_design = make_full_design(best_design, qubit_geo, resonator_geo, feedline_geo)
    pin_coords = resonator_design.components['transmon'].pins['readout']['middle']
    otg2 = OpenToGround(
            resonator_design,
            'otg2',
            options=dict(
                    chip='main',
                    layer=resonator_design.components['resonator'].options["layer"], # Ensure OTG doesnt make its own ground layer by putting it on the same layer as resonator
                    pos_x=f'{pin_coords[0]}mm',
                    pos_y=f'{pin_coords[1]}mm',
                    orientation=270,
                    termination_gap = resonator_geo["trace_gap"],
            )
    )
    resonator_design.delete_component('LP1')
    resonator_design.delete_component('LP2')
    resonator_design.delete_component('feedline')
    resonator_design.delete_component('transmon')
    resonator_design.delete_component('junction')
    _center_and_scale_chip(resonator_design, component_name='resonator', scale_factor=scaling_factor)
    return resonator_design

def make_feedline_design(best_design, qubit_geo, resonator_geo, feedline_geo, scaling_factor=1.3):
    feedline_design = make_full_design(best_design, qubit_geo, resonator_geo, feedline_geo)
    feedline_design.delete_component('otg1')
    feedline_design.delete_component('resonator')
    feedline_design.delete_component('transmon')
    feedline_design.delete_component('junction')
    _center_and_scale_chip(feedline_design, component_name='feedline', scale_factor=scaling_factor)
    return feedline_design



