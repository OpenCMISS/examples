gfx read node QuadraticTriangleSimplex.part0.exnode
gfx read element QuadraticTriangleSimplex.part0.exelem generate
gfx create window 1
gfx modify g_element "Region 2" lines select_on material default selected_material default_selected
gfx modify g_element "Region 2" surfaces select_on material default selected_material default_selected render_shaded
gfx 
