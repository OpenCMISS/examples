gfx read node ExtracellularBidomain.part0.exnode;
gfx read elem ExtracellularBidomain.part0.exelem;

gfx define faces egroup ExtracellularBidomainRegion;

gfx create spectrum phi;
gfx create spectrum vm;
#gfx modify spectrum phi linear reverse range 0.0 14.0 extend_above extend_below rainbow colour_range 0 1 component 1;
gfx modify spectrum phi linear reverse range 0.0 66.0 extend_above extend_below rainbow colour_range 0 1 component 1;
gfx modify spectrum phi autorange;
#gfx modify spectrum phi linear reverse range -2.0 30.0 extend_above extend_below banded number_of_bands 20 band_ratio 0.06 component 1;
gfx modify spectrum vm linear reverse range -82.0 33.0 extend_above extend_below rainbow colour_range 0 1 component 1;
gfx modify spectrum vm autorange;
#gfx modify spectrum vm linear reverse range -70.0 30.0 extend_above extend_below banded number_of_bands 20 band_ratio 0.06 component 1;

gfx define tessellation default minimum_divisions "1" refinement_factors "4";
gfx define tessellation tess20 minimum_divisions "20" refinement_factors "20";

#use "Coordinate" instead of "Geometry" if generated mesh is used
gfx modify g_element "/" general clear;
gfx modify g_element "/" node_points subgroup ExtracellularBidomainRegion coordinate Geometry LOCAL glyph sphere general size "0.01*0.01*0.01" centre 0,0,0 font default select_on material black selected_material default_selected;
gfx modify g_element "/" lines subgroup ExtracellularBidomainRegion coordinate Geometry tessellation default LOCAL select_on material black selected_material default_selected;
gfx modify g_element "/" surfaces subgroup ExtracellularBidomainRegion coordinate Geometry tessellation tess20 LOCAL select_on material default data Phi spectrum phi selected_material default_selected render_shaded;
#gfx modify g_element "/" surfaces subgroup ExtracellularBidomainRegion coordinate Geometry tessellation tess20 LOCAL select_on material default data Vm spectrum vm selected_material default_selected render_shaded;
gfx modify g_element "/" point as axes coordinate Geometry LOCAL glyph axes general size "1*1*1" centre 0,0,0 font default select_on material default selected_material default;


gfx create window 1 double_buffer;
gfx modify window 1 image scene default light_model default;
gfx modify window 1 image add_light default;
gfx modify window 1 layout simple ortho_axes z -y eye_spacing 0.25 width 1452 height 543;
gfx modify window 1 set current_pane 1;
gfx modify window 1 background colour 1 1 1 texture none;
gfx modify window 1 view parallel eye_point 4.39635 4.01765 16.2323 interest_point 6.79636 0.314534 0 up_vector 0.0785339 0.974395 -0.21068 view_angle 28.5453 near_clipping_plane 0.168214 far_clipping_plane 60.114 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;
#gfx modify window 1 set transform_tool current_pane 1 std_view_angle 40 normal_lines no_antialias depth_of_field 0.0 fast_transparency blend_normal;



gfx create colour_bar axis 0 1 0 label_material black number_format [%+.2e] spectrum phi tick_direction 1 0 0;
#gfx create colour_bar axis 0 1 0 label_material black number_format [%+.2e] spectrum vm tick_direction 1 0 0;
gfx modify g_element "/" point coordinate Geometry NORMALISED_WINDOW_FIT_LEFT glyph colour_bar general size "0.5*0.5*0.5" centre 0.5,0,0 font default select_on material default selected_material default_selected;

gfx edit scene;





# get com commands:
# gfx list g_elem com
# gfx list win com

#gfx define field annotation string_constant "potential phi!" material black;
#gfx create colour_bar as legend_iso axis 0 0 -1 centre -1 1.5 0.5 number_format [%+.3e] label_material black spectrum STRESS
#gfx create annotation as text_iso1 material black position -1 1.5 1.5 text "von Mises"
#gfx create annotation as text_iso2 material black position -1 1.5 -0.5 text "[N/mm^2]"
#gfx draw legend_iso
#gfx draw text_iso1
#gfx draw text_iso2 

