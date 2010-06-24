gfx read node LinearTetrahedraSimplex.part0.exnode
gfx read element LinearTetrahedraSimplex.part0.exelem
gfx define faces egroup "Region 2"
gfx create window 1
gfx modify window 1 view interest_point 1.0,0.5,0.0 eye_point 1.0,0.5,5.0 up_vector 0.0,1.0,0.0
gfx modify g_element "Region 2" lines select_on 
#gfx modify g_element "Region 2" surfaces select_on 
gfx modify g_element "Region 2" node_points select_on material red label cmiss_number
gfx modify g_element "Region 2" element_points select_on material blue label cmiss_number use_elements cell_centres discretization "1*1*1"
gfx 
