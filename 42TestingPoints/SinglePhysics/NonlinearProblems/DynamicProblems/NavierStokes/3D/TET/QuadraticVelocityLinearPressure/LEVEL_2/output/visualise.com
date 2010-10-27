#Read in the sequence of nodal positions.
for $i (1..1)
  {
	 $filename = sprintf("./output/TIME_STEP_%04d.exnode", $i);
	 print "Reading $filename time $i\n";
	 gfx read node "$filename" time $i;
  }

#Read in the element description
gfx read elem ./output/TIME_STEP_0001.exelem;
gfx create window 1

gfx cre spectrum flow
gfx modify spectrum flow clear overwrite_colour;
gfx modify spectrum flow linear reverse range 1 1.25 extend_above extend_below rainbow colour_range 0 1 component 2;

gfx def field vector_field coord rectangular_cartesian component general.1 general.2 general.3

gfx cre spectrum pressure
gfx modify spectrum pressure linear reverse range -1 1 extend_above extend_below rainbow colour_range 0 1 component 4;

gfx modify g_element OpenCMISS node_points glyph arrow_solid general size "0.1*0.1*0.1" centre 0,0,0 select_on material default selected_material default_selected data vector_field orientation vector_field scale_factors "0.1*0.1*0.1" spectrum flow

gfx modify window 1 background colour 1 1 1

gfx define faces egroup OpenCMISS
gfx modify g_element OpenCMISS surfaces select_on material default selected_material default_selected data general spectrum pressure

#Set the timekeeper playing
gfx timekeeper default play speed 1 skip;
gfx create time_editor


gfx define field p_field coordinate_system rectangular_cartesian composite general.4;
gfx define field vector_field coordinate_system rectangular_cartesian composite general.1 general.2 general.3;
gfx define field mag_field coordinate_system rectangular_cartesian magnitude field vector_field;
gfx modify spectrum default clear overwrite_colour;
gfx modify spectrum default linear reverse range -34.946 50.493 extend_above extend_below rainbow colour_range 0 1 component 1;
gfx modify spectrum flow clear overwrite_colour;
gfx modify spectrum flow linear reverse range 1 1.25 extend_above extend_below rainbow colour_range 0 1 component 2;
gfx modify spectrum pressure clear overwrite_colour;
gfx modify spectrum pressure linear reverse range 0 1 extend_above extend_below rainbow colour_range 0 1 component 1;
gfx modify spectrum pressure linear reverse range -1 1 extend_above extend_below rainbow colour_range 0 1 component 4;
gfx create material black normal_mode ambient 0 0 0 diffuse 0 0 0 emission 0 0 0 specular 0.3 0.3 0.3 alpha 1 shininess 0.2;
gfx create material blue normal_mode ambient 0 0 0.5 diffuse 0 0 1 emission 0 0 0 specular 0.2 0.2 0.2 alpha 1 shininess 0.2;
gfx create material bone normal_mode ambient 0.7 0.7 0.6 diffuse 0.9 0.9 0.7 emission 0 0 0 specular 0.1 0.1 0.1 alpha 1 shininess 0.2;
gfx create material default normal_mode ambient 1 1 1 diffuse 1 1 1 emission 0 0 0 specular 0 0 0 alpha 1 shininess 0;
gfx create material default_selected normal_mode ambient 1 0.2 0 diffuse 1 0.2 0 emission 0 0 0 specular 0 0 0 alpha 1 shininess 0;
gfx create material gold normal_mode ambient 1 0.4 0 diffuse 1 0.7 0 emission 0 0 0 specular 0.5 0.5 0.5 alpha 1 shininess 0.3;
gfx create material gray50 normal_mode ambient 0.5 0.5 0.5 diffuse 0.5 0.5 0.5 emission 0.5 0.5 0.5 specular 0.5 0.5 0.5 alpha 1 shininess 0.2;
gfx create material green normal_mode ambient 0 0.5 0 diffuse 0 1 0 emission 0 0 0 specular 0.2 0.2 0.2 alpha 1 shininess 0.1;
gfx create material muscle normal_mode ambient 0.4 0.14 0.11 diffuse 0.5 0.12 0.1 emission 0 0 0 specular 0.3 0.5 0.5 alpha 1 shininess 0.2;
gfx create material red normal_mode ambient 0.5 0 0 diffuse 1 0 0 emission 0 0 0 specular 0.2 0.2 0.2 alpha 1 shininess 0.2;
gfx create material silver normal_mode ambient 0.4 0.4 0.4 diffuse 0.7 0.7 0.7 emission 0 0 0 specular 0.5 0.5 0.5 alpha 1 shininess 0.3;
gfx create material tissue normal_mode ambient 0.9 0.7 0.5 diffuse 0.9 0.7 0.5 emission 0 0 0 specular 0.2 0.2 0.3 alpha 1 shininess 0.2;
gfx create material transparent_gray50 normal_mode ambient 0.5 0.5 0.5 diffuse 0.5 0.5 0.5 emission 0.5 0.5 0.5 specular 0.5 0.5 0.5 alpha 0 shininess 0.2;
gfx create material white normal_mode ambient 1 1 1 diffuse 1 1 1 emission 0 0 0 specular 0 0 0 alpha 1 shininess 0;
gfx modify g_element OpenCMISS general clear circle_discretization 6 default_coordinate coordinates element_discretization "4*4*4" native_discretization none;
gfx modify g_element OpenCMISS lines select_on material default selected_material default_selected;
gfx modify g_element OpenCMISS node_points glyph arrow_solid general size "0.1*0.1*0.1" centre 0,0,0 font default orientation vector_field scale_factors "0.1*0.1*0.1" select_on material default data mag_field spectrum default selected_material default_selected;
gfx modify g_element OpenCMISS surfaces select_on material default data p_field spectrum default selected_material default_selected render_shaded;
gfx modify window 1 image scene default light_model default;
gfx modify window 1 image add_light default;
gfx modify window 1 layout simple ortho_axes z -y eye_spacing 0.25 width 512 height 512;
gfx modify window 1 set current_pane 1;
gfx modify window 1 background colour 1 1 1 texture none;
gfx modify window 1 view parallel eye_point -5.14426 -26.338 22.9826 interest_point 0 0 0 up_vector 0.282893 0.598661 0.749384 view_angle 40 near_clipping_plane 0.353321 far_clipping_plane 126.265 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;
gfx modify window 1 overlay scene none;
gfx modify window 1 set transform_tool current_pane 1 std_view_angle 40 normal_lines no_antialias depth_of_field 0.0 fast_transparency blend_normal;
