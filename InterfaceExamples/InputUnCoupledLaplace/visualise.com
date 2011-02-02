gfx read node CoupledLaplace_1.part0.exnode 
gfx read element CoupledLaplace_1.part0.exelem 
gfx define faces egroup Region1

gfx change_identifier no +100000
gfx change_identifier ele +100000


gfx read node CoupledLaplace_2.part0.exnode
gfx read element CoupledLaplace_2.part0.exelem
gfx define faces egroup Region2

gfx change_identifier no +100000
gfx change_identifier ele +100000

gfx read node CoupledLaplace_Interface.part0.exnode
gfx read element CoupledLaplace_Interface.part0.exelem
gfx define faces egroup Interface

gfx create window 1
gfx modify window 1 view interest_point 2.0,0.5,0.0 eye_point 2.0,0.5,10.0 up_vector 0.0,1.0,0.0
gfx modify spectrum default clear overwrite_colour;
gfx modify spectrum default linear reverse range -1.0 1.0 extend_above extend_below rainbow colour_range 0 1 component 1
gfx modify spectrum default linear reverse range -1.0 1.0 extend_above extend_below banded number_of_bands 20 band_ratio 0.06 component 1
gfx modify g_element Region1 node_points glyph sphere general size "0.1*0.1*0.1" centre 0,0,0 
gfx modify g_element Region2 node_points glyph sphere general size "0.1*0.1*0.1" centre 0,0,0 

gfx modify g_element Interface node_points glyph sphere general size "0.1*0.1*0.1" centre 0,0,0 

gfx modify g_element Region1 lines select_on material default selected_material default_selected
gfx modify g_element Region2 lines select_on material default selected_material default_selected

gfx modify g_element Interface lines select_on material default data general spectrum default selected_material default_selected
gfx modify g_element Interface element_points glyph point general size "1*1*1" centre 0,0,0 label general use_elements cell_centres discretization "2*2*2" native_discretization NONE select_on material default selected_material default_selected
gfx modify g_element Region1 surfaces select_on material default data general spectrum default selected_material default_selected render_shaded
gfx modify g_element Region2 surfaces select_on material default data general spectrum default selected_material default_selected render_shaded
