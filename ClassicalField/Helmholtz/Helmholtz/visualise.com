# read in description
gfx read node Helmholtz.part0.exnode
gfx read element Helmholtz.part0.exelem generate

gfx define field "Region 1/amplitude" component general.1

gfx create window 1

gfx create spectrum amplitude
gfx modify spectrum amplitude clear overwrite_colour
gfx modify spectrum amplitude linear reverse range -1 1 extend_above extend_below rainbow colour_range 0 1

# display geometry with amplitude
gfx modify g_element "Region 1" surfaces coordinate coordinates select_on material blue selected_material default_selected data amplitude spectrum amplitude render_shaded
gfx modify g_element "Region 1" lines coordinate coordinates select_on material default selected_material default_selected

gfx modify window 1 view interest_point 1.0,0.5,0.0 eye_point 1.0,0.5,5.0 up_vector 0.0,1.0,0.0
