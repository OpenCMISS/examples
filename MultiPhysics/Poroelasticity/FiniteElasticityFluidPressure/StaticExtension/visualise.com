system("sed -i 's%unknown field variable type%rectangular cartesian%' StaticExtension.part0.exnode");
system("sed -i 's%unknown field variable type%rectangular cartesian%' StaticExtension.part0.exelem");

# read in description
gfx read node StaticExtension.part0.exnode
gfx read element StaticExtension.part0.exelem

# define deformed geometry and pressure fields
gfx define field "deformed_geom" component U.1 U.2 U.3
gfx define field "fluid_pressure" component V.1

gfx create window 1

gfx create spectrum fluid
gfx modify spectrum fluid clear overwrite_colour
gfx modify spectrum fluid linear reverse range 0.0 1.0 extend_above extend_below rainbow colour_range 0 1

# display deformed geometry
gfx define faces egroup "Region 1"
gfx modify g_element "Region 1" surfaces coordinate deformed_geom select_on material tissue selected_material default_selected data fluid_pressure spectrum fluid render_shaded
gfx modify g_element "Region 1" lines coordinate deformed_geom line_width 2 select_on material default selected_material default_selected

gfx modify g_element "Region 1" node_points coordinate deformed_geom glyph sphere general size "0.5*0.5*0.5" centre 0,0,0 font default select_on material default selected_material default_selected

# display undeformed nodes
gfx modify g_element "Region 1" lines select_on material green line_width 2 selected_material default_selected

gfx create axes length 12.0 material default
gfx draw axes

gfx modify spectrum fluid autorange

gfx edit scene
gfx modify window 1 set antialias 8
