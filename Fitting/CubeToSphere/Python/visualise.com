gfx read node UndeformedGeometry.part0.exnode time 0
gfx read node DeformedGeometry1.part0.exnode time 1
#gfx read node DeformedGeometry2.part0.exnode time 2

gfx read elem UndeformedGeometry.part0.exelem

gfx def faces egroup FittingRegion
gfx modify g_element "/" lines coordinate Coordinate tessellation default LOCAL native_discretization NONE select_on material default selected_material default_selected;
gfx modify g_element "/" surfaces coordinate Coordinate exterior tessellation default LOCAL native_discretization NONE select_on material muscle selected_material default_selected render_shaded;
gfx modify g_element "/" node_points coordinate Coordinate LOCAL glyph sphere general size "0.1*0.1*0.1" centre 0,0,0 font default select_on material default selected_material default_selected;

gfx read data DataPoints.part0.exdata
#gfx read data DataPoints.part1.exdata
gfx modify g_element DataPoints data_points coordinate data_coordinates LOCAL glyph point general size "0.1*0.1*0.1" centre 0,0,0 font default select_on material blue selected_material default_selected;
gfx cre wi
gfx cre time_editor
gfx edit scene
