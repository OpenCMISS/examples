gfx read node Embedded1DLaplace.part0.exnode
gfx read elem Embedded1DLaplace.part0.exelem
gfx define faces egroup LaplaceRegion
gfx create window 1
gfx draw axes
gfx modify window 1 view interest_point 1.0,0.5,0.0 eye_point 1.0,0.5,5.0 up_vector 0.0,1.0,0.0
gfx modify spectrum default clear overwrite_colour;
gfx modify spectrum default linear reverse range 0.0 1.0 extend_above extend_below rainbow colour_range 0 1 component 1
gfx modify spectrum default linear reverse range 0.0 1.0 extend_above extend_below banded number_of_bands 20 band_ratio 0.06 component 1
gfx modify g_element LaplaceRegion node_points glyph sphere general size "0.1*0.1*0.1" centre 0,0,0 
gfx modify g_element LaplaceRegion lines select_on material default data Phi spectrum default selected_material default_selected
gfx 
