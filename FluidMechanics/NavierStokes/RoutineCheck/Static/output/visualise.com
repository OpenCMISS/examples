gfx read node './output/STATICSOLUTION.exnode'
gfx read elem './output/STATICSOLUTION.exelem'
gfx create window 1

gfx cre spectrum flow
gfx modify spectrum flow clear overwrite_colour;
gfx modify spectrum flow linear reverse range 1 1.25 extend_above extend_below rainbow colour_range 0 1 component 2;

gfx def field vector_field coord rectangular_cartesian component general.1 general.2 general.3

gfx cre spectrum pressure
gfx modify spectrum pressure linear reverse range 5 15 extend_above extend_below rainbow colour_range 0 1 component 4;
gfx modify spectrum pressure overlay_colour;


gfx modify window 1 background colour 1 1 1

gfx define faces egroup OpenCMISS
gfx create material default_transparent normal_mode ambient 0 0 0 diffuse 0 0 0 emission 0 0 0 specular 0.3 0.3 0.3 alpha 0.5 shininess 0.2;
gfx modify g_element OpenCMISS surfaces select_on material default_transparent selected_material default_selected data general spectrum pressure

gfx modify g_element OpenCMISS element_points glyph arrow_solid general size "0*0.01*0.01" centre 0,0,0 font default orientation vector_field scale_factors "0.1*0*0" use_elements cell_centres discretization "3*3*3" native_discretization NONE select_on material gold selected_material default_selected;


gfx modify window 1 set order_independent_transparency 5;
gfx modify window 1 view parallel eye_point -0.303916 0.0379809 -3.87343 interest_point 0 0 0 up_vector 0.411286 0.911208 -0.0233354 view_angle 40 near_clipping_plane 0.0388552 far_clipping_plane 13.8855 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;
