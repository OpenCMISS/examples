gfx read node LagrangeSimplexMeshExample.part0.exnode
gfx read elem LagrangeSimplexMeshExample.part0.exelem

gfx define faces egroup "Region 1"
gfx create window 1
gfx modify window 1 view interest_point 1.0,0.5,0.0 eye_point 1.0,0.5,5.0 up_vector 0.0,1.0,0.0
gfx modify g_element "Region 1" node_points glyph sphere general size "0.05*0.05*0.05" centre 0,0,0 
gfx modify g_element "Region 1" lines select_on material default selected_material default_selected

gfx draw axes
gfx edit scene

