# read in description
gfx read node MultipleMeshComponents.part0.exnode
gfx read element MultipleMeshComponents.part0.exelem

gfx create window 1

gfx create material surface_material ambient 0.2 0.2 0.8 diffuse 0.2 0.2 0.8 alpha 0.4

# display geometry
gfx define faces egroup "Region"
gfx modify g_element "Region" lines select_on material default selected_material default_selected
gfx modify g_element "Region" node_points glyph point label cmiss_number font default select_on material default selected_material default_selected
gfx modify g_element "Region" surfaces select_on material surface_material render_shaded

# display second field
gfx modify g_element "Region" surfaces select_on material tissue selected_material default_selected data "Extra Field" spectrum default render_shaded

gfx modify g_element "/" point  glyph axes general size "10*10*10" centre 0,0,0 font default select_on material default selected_material default_selected

gfx modify window 1 image view_all;
