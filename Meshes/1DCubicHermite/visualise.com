gfx read node 1DCubicHermite.part0.exnode
gfx read element 1DCubicHermite.part0.exelem
gfx define faces egroup 1DCubicHermiteRegion
gfx create window 1
gfx modify window 1 view interest_point 2.5,0.0,0.0 eye_point 2.5,0.0,15.0 up_vector 0.0,1.0,0.0
gfx modify g_element 1DCubicHermiteRegion lines invisible
gfx modify g_element 1DCubicHermiteRegion general default_coordinate UnitScaling element_discretization 100*100*100
gfx modify g_element 1DCubicHermiteRegion lines as Unit coordinate UnitScaling select_on material white
gfx modify g_element 1DCubicHermiteRegion node_points select_on label cmiss_number
gfx modify g_element 1DCubicHermiteRegion element_points select_on label cmiss_number use_elements cell_centres discretization "1*1*1"
gfx modify g_element 1DCubicHermiteRegion lines as ArithmeticMean coordinate ArithmeticMeanScaling select_on material red
gfx modify g_element 1DCubicHermiteRegion lines as GeometricMean coordinate GeometricMeanScaling select_on material blue
gfx modify g_element 1DCubicHermiteRegion lines as HarmonicMean coordinate HarmonicMeanScaling select_on material green
gfx 
