# read in solution, which may be split into multiple files
@exnodes=<./Cantilever.part*.exnode>;
@exelems=<./Cantilever.part*.exelem>;
foreach $filename (@exnodes) {
    print "Reading $filename\n";
    gfx read node "$filename";
}
foreach $filename (@exelems) {
    print "Reading $filename\n";
    gfx read elem "$filename";
}

# define deformed geometry and pressure
gfx define field "deformed_geom" component Dependent.1 Dependent.2 Dependent.3
gfx define field "hydrostatic_pressure" component Dependent.4

gfx create window 1

gfx create spectrum pressure
gfx modify spectrum pressure clear overwrite_colour
gfx modify spectrum pressure linear reverse range 0.0 1.0 extend_above extend_below rainbow colour_range 0 1

# display deformed geometry
gfx define faces egroup "Region 1"
gfx modify g_element "Region 1" lines coordinate deformed_geom select_on material default selected_material default_selected
j
gfx modify g_element "Region 1" surfaces coordinate deformed_geom select_on material tissue selected_material default_selected data hydrostatic_pressure spectrum pressure render_shaded

gfx modify g_element "Region 1" node_points coordinate deformed_geom glyph sphere General size "2*2*2" centre 0,0,0 font default select_on material default selected_material default_selected

# display undeformed lines
gfx modify g_element "Region 1" lines select_on material green selected_material default_selected

gfx modify g_element "/" point  glyph axes general size "75*60*60" centre 0,0,0 font default select_on material default selected_material default_selected;

gfx modify spectrum pressure autorange
gfx create colour_bar spectrum pressure material default
gfx modify g_element "/" point NORMALISED_WINDOW_FIT_LEFT glyph colour_bar general size "0.5*0.5*0.2" centre -0.4,0.4,0.0 select_on material default selected_material default;

gfx edit scene
gfx modify window 1 set antialias 2
gfx modify window 1 view parallel eye_point 20 -200 20 interest_point 20 20 20 up_vector 0 0 1 view_angle 40 near_clipping_plane 1.5 far_clipping_plane 700 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1
