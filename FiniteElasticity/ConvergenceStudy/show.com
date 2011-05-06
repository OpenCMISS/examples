gfx re no ConvergenceStudy.part0.exnode
gfx re no ConvergenceStudy.part1.exnode
gfx re el ConvergenceStudy.part0.exelem generate
gfx re el ConvergenceStudy.part1.exelem generate

# define fields: pressure, deformed
gfx define field pressure coordinate_system rectangular_cartesian composite general.4;
gfx define field deformed_coordinates coordinate_system rectangular_cartesian composite general.1 general.2 general.3;

# display commands
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
gfx create material bluey ambient 0 0.25 0.5 diffuse 0 0.4 1 specular 0.5 0.5 0.5 shininess 0.3
gfx create material jade ambient 0.1 0.56 0 diffuse 0.1 0.39 0 emission 0 0 0 specular 0.8 0.8 0.8 alpha 1 shininess 0.8

gfx modify g_element "Region 1" general clear circle_discretization 6 default_coordinate coordinates element_discretization "20*20*20" native_discretization none;
gfx modify g_element "Region 1" lines select_on material white selected_material default_selected;
gfx modify g_element "Region 1" lines coordinate deformed_coordinates select_on material green selected_material default_selected;
gfx modify g_element "Region 1" node_points glyph point general size "1*1*1" centre 0,0,0 font default label cmiss_number select_on material gold selected_material default_selected;
gfx modify g_element "Region 1" surfaces coordinate deformed_coordinates select_on invisible material default data pressure spectrum default selected_material default_selected render_shaded;

# window and scene commands
gfx create window 1 double_buffer;
gfx modify window 1 image scene default light_model default;
gfx modify window 1 image add_light default;
gfx modify window 1 layout simple ortho_axes z -y eye_spacing 0.25 width 519 height 682;
gfx modify window 1 set current_pane 1;
gfx modify window 1 background colour 0 0 0 texture none;
gfx modify window 1 view parallel eye_point -0.828084 -2.02811 2.66528 interest_point 0.542705 0.587444 0.725341 up_vector 0.0677115 0.571217 0.818001 view_angle 40 near_clipping_plane 0.0353321 far_clipping_plane 12.6265 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;
gfx modify window 1 overlay scene none;
gfx modify window 1 set transform_tool current_pane 1 std_view_angle 40 normal_lines no_antialias depth_of_field 0.0 fast_transparency blend_normal;

# autorange the spectrum
gfx mod spect default autorange

#######################################################

  ## show strains ##

# Define fields necessary to show the strains.
#
# Define F = dx/dX, where X are undeformed coordinates, x are deformed coordinates.
# The 9 values are returned in the order dx/dX dx/dY dx/dZ dy/X etc., ie. across
# the rows of a 3 x 3 matrix first to present a single vector.
# Then:
#   C = F_transpose . F     = right Cauchy-Green deformation tensor
#   E = 0.5*(C - Identity3) = Lagrangian / Greens [finite] strain tensor
#
gfx define field F gradient coordinate coordinates field deformed_coordinates
gfx define field F_transpose transpose source_number_of_rows 3 field F
gfx define field Identity3 composite 1 0 0 0 1 0 0 0 1
gfx define field C matrix_multiply number_of_rows 3 fields F_transpose F
gfx define field E2 add fields C Identity3 scale_factors 1 -1
gfx define field E scale field E2 scale_factors 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5

# Obtain the principal strains and their direction vectors by performing
# an eigenanalysis on E.
#
gfx define field principal_strains eigenvalues field E
gfx define field principal_strain_vectors coordinate_system rectangular_cartesian eigenvectors eigenvalues principal_strains
#
# The above vectors mark the orientation of the principal strains in the
# undeformed configuration, where the first 3 values in the 9-component field
# are the components of the first eigenvector in the principal_strains field,
# and so on. We convert to deformed directions by post-multiplying them by
# F_transpose. This scales the vectors by the stretches in the principal
# directions, which are the ratios of final length to initial length. The
# deformed principal strain directions remain orthogonal by definition.
#
gfx define field deformed_principal_strain_vectors coordinate_system rectangular_cartesian matrix_multiply number_of_rows 3 fields principal_strain_vectors F_transpose

# Extract individual eigenvalues and eigenvectors for display.
#
gfx define field principal_strain1 composite principal_strains.1
gfx define field principal_strain2 composite principal_strains.2
gfx define field principal_strain3 composite principal_strains.3
gfx define field principal_strain_vector1 coordinate_system rectangular_cartesian composite principal_strain_vectors.1 principal_strain_vectors.2 principal_strain_vectors.3
gfx define field principal_strain_vector2 coordinate_system rectangular_cartesian composite principal_strain_vectors.4 principal_strain_vectors.5 principal_strain_vectors.6
gfx define field principal_strain_vector3 coordinate_system rectangular_cartesian composite principal_strain_vectors.7 principal_strain_vectors.8 principal_strain_vectors.9
gfx define field deformed_principal_strain_vector1 coordinate_system rectangular_cartesian composite deformed_principal_strain_vectors.1 deformed_principal_strain_vectors.2 deformed_principal_strain_vectors.3
gfx define field deformed_principal_strain_vector2 coordinate_system rectangular_cartesian composite deformed_principal_strain_vectors.4 deformed_principal_strain_vectors.5 deformed_principal_strain_vectors.6
gfx define field deformed_principal_strain_vector3 coordinate_system rectangular_cartesian composite deformed_principal_strain_vectors.7 deformed_principal_strain_vectors.8 deformed_principal_strain_vectors.9
#
# Normalise the deformed principal strain vectors to remove stretches.
#
gfx define field norm_def_principal_strain_vector1 coordinate_system rectangular_cartesian normalise field deformed_principal_strain_vector1
gfx define field norm_def_principal_strain_vector2 coordinate_system rectangular_cartesian normalise field deformed_principal_strain_vector2
gfx define field norm_def_principal_strain_vector3 coordinate_system rectangular_cartesian normalise field deformed_principal_strain_vector3

gfx create spectrum strain
gfx modify spectrum strain clear overwrite_colour
gfx modify spectrum strain linear reverse range -1 0 extend_below red colour_range 1 1 ambient diffuse component 1
gfx modify spectrum strain linear reverse range 0 1 extend_above blue colour_range 1 1 ambient diffuse component 1
gfx modify spectrum strain linear reverse range 0 1 extend_above green colour_range 0.5 0.5 ambient diffuse component 1
#
# Make the deformed coordinates the default coordinates for the cube, and
# display lines and node points (to edit). At any time in the future you
# can override the coordinate field in any settings for the cube in the
# graphical element editor to show the undeformed state.
#
gfx modify g_element "Region 1" general clear circle_discretization 6 default_coordinate deformed_coordinates element_discretization "4*4*4" native_discretization none
gfx modify g_element "Region 1" lines select_on material default selected_material default_selected
gfx modify g_element "Region 1" node_points glyph sphere general size "0.05*0.05*0.05" centre 0,0,0 select_on material default selected_material default_selected
#
# Also show undeformed lines to see what you are deforming from.
#
gfx modify g_element "Region 1" lines coordinate coordinates select_on material jade selected_material default_selected
#
# Display eigenvalue/vector sets 1, 2 and 3. Each command ensures the strains
# are shown at the deformed coordinates and oriented with the normalised
# deformed principal strain vectors. The corresponding principal strain
# eigenvalue is supplied as the variable_scale field. This, in combination with
# the mirrow glyph (a cone in this place) allows negative values to cause the
# glyph to point inwards, thus indicating compressive strains.
#
# Strain glyphs are shown at 3 x 3 x 3 locations throughout the cube and each
# vector is coloured by the principal_strain value with the strain spectrum.
# The base "size" of 0*0.02*0.02 indicates that the cones are 0.02 units in
# diameter but have no base length. Scale factors 0.5*0*0 make the glyph half
# as long as the magnitude of the variable_scale (since no scaling is coming
# from the normalised vectors). Their length is thus proportional to the
# principal strain in that direction.
#
gfx modify g_element "Region 1" element_points coordinate deformed_coordinates glyph mirror_cone size "0*0.02*0.02" centre 0,0,0 orientation norm_def_principal_strain_vector1 variable_scale principal_strain1 scale_factors "0.5*0*0" use_elements cell_centres discretization "3*3*3" native_discretization NONE select_on material bluey data principal_strain1 spectrum strain selected_material default_selected
gfx modify g_element "Region 1" element_points coordinate deformed_coordinates glyph mirror_cone size "0*0.02*0.02" centre 0,0,0 orientation norm_def_principal_strain_vector2 variable_scale principal_strain2 scale_factors "0.5*0*0" use_elements cell_centres discretization "3*3*3" native_discretization NONE select_on material bluey data principal_strain2 spectrum strain selected_material default_selected
gfx modify g_element "Region 1" element_points coordinate deformed_coordinates glyph mirror_cone size "0*0.02*0.02" centre 0,0,0 orientation norm_def_principal_strain_vector3 variable_scale principal_strain3 scale_factors "0.5*0*0" use_elements cell_centres discretization "3*3*3" native_discretization NONE select_on material bluey data principal_strain3 spectrum strain selected_material default_selected


# Write the principal strains at the centre of the undeformed cube.
#
gfx modify g_element "Region 1" element_points coordinate coordinates glyph point general size "0*0.02*0.02" centre 0,0,0 label principal_strains use_elements cell_centres discretization "1*1*1" native_discretization NONE select_on material default selected_material default_selected


















