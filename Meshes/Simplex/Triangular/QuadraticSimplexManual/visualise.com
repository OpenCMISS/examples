gfx read node SimplexMesh.part0.exnode
gfx read element SimplexMesh.part0.exelem generate
gfx create window 1
gfx modify window 1 view interest_point 1.0,0.5,0.0 eye_point 1.0,0.5,5.0 up_vector 0.0,1.0,0.0
gfx modify g_element "Region 2" lines select_on material default selected_material default_selected
gfx modify g_element "Region 2" surfaces select_on material default selected_material default_selected render_shaded
gfx 
