gfx re no TriQuadraticCube.part0.exnode
gfx re el TriQuadraticCube.part0.exelem generate
gfx def field deformed composite general.1 general.2 general.3

gfx define field pressure coordinate_system rectangular_cartesian composite general.4;
gfx define field deformed coordinate_system rectangular_cartesian composite general.1 general.2 general.3;
gfx modify spectrum default clear overwrite_colour;
gfx modify spectrum default linear reverse range -13.6372 -12.5593 extend_above extend_below rainbow colour_range 0 1 component 1;
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
gfx modify g_element "Region 1" general clear circle_discretization 6 default_coordinate coordinates element_discretization "20*20*20" native_discretization none;
gfx modify g_element "Region 1" lines select_on material white selected_material default_selected;
gfx modify g_element "Region 1" lines coordinate deformed select_on material green selected_material default_selected;
gfx modify g_element "Region 1" node_points glyph point general size "1*1*1" centre 0,0,0 font default label cmiss_number select_on material gold selected_material default_selected;
gfx modify g_element "Region 1" surfaces coordinate deformed select_on invisible material default data pressure spectrum default selected_material default_selected render_shaded;
gfx create window 1 double_buffer;
gfx modify window 1 image scene default light_model default;
gfx modify window 1 image add_light default;
gfx modify window 1 layout simple ortho_axes z -y eye_spacing 0.25 width 519 height 682;
gfx modify window 1 set current_pane 1;
gfx modify window 1 background colour 0 0 0 texture none;
gfx modify window 1 view parallel eye_point -0.828084 -2.02811 2.66528 interest_point 0.542705 0.587444 0.725341 up_vector 0.0677115 0.571217 0.818001 view_angle 40 near_clipping_plane 0.0353321 far_clipping_plane 12.6265 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;
gfx modify window 1 overlay scene none;
gfx modify window 1 set transform_tool current_pane 1 std_view_angle 40 normal_lines no_antialias depth_of_field 0.0 fast_transparency blend_normal;

gfx mod spect default autorange