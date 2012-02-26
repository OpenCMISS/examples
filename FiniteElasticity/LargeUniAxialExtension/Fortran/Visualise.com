# read in Region 1 description
gfx read node LargeUniaxialExtension.part0.exnode
gfx read element LargeUniaxialExtension.part0.exelem

# define deformed geometry
gfx define field "deformed_geom" component Dependent.1 Dependent.2 Dependent.3

# display deformed geometry
gfx define faces egroup "Region"
gfx modify g_element "Region" lines coordinate deformed_geom select_on material default selected_material default_selected
#gfx modify g_element "Region" surfaces coordinate deformed_geom select_on material tissue selected_material default_selected render_shaded
gfx modify g_element "Region" node_points coordinate deformed_geom glyph sphere General size "0.1*0.1*0.1" centre 0,0,0 font default select_on material default selected_material default_selected

# display undeformed lines
gfx modify g_element "Region" lines select_on material green selected_material default_selected

gfx create window 1

gfx create axes length 5 material default
gfx draw axes

gfx edit scene
gfx modify window 1 set antialias 2
