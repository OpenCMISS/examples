# read in description
gfx read node Cantilever.part0.exnode
gfx read element Cantilever.part0.exelem

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

gfx create axes length 20 material default
gfx draw axes

gfx modify spectrum pressure autorange

gfx edit scene
gfx modify window 1 set antialias 2
