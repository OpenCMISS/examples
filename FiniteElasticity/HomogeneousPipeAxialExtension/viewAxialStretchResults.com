gfx read elem AxialStretch.part0.exelem;
gfx read node AxialStretch.part0.exnode;

gfx define faces;
#Get the deformed coordinate field
gfx define field deform composite Dependent.1 Dependent.2 Dependent.3

gfx edit scene

#Calculate the strains
gfx define field F gradient coordinate coordinates field deform
gfx define field F_transpose transpose source_number_of_rows 3 field F
gfx define field identity3 composite 1 0 0 0 1 0 0 0 1

gfx define field C matrix_multiply number_of_rows 3 fields F_transpose F
gfx define field E2 add fields C identity3 scale_factors 1 -1
gfx define field E coordinate_system rectangular_cartesian scale field E2 scale_factors 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5

gfx define field principal_strains eigenvalues field E
gfx define field principal_strain_vectors eigenvectors eigenvalues principal_strains
gfx define field deformed_principal_strain_vectors matrix_multiply number_of_rows 3 fields principal_strain_vectors F_transpose

gfx define field deformed_principal_strain_vector1 composite deformed_principal_strain_vectors.1 deformed_principal_strain_vectors.2 deformed_principal_strain_vectors.3
gfx define field deformed_principal_strain_vector2 composite deformed_principal_strain_vectors.4 deformed_principal_strain_vectors.5 deformed_principal_strain_vectors.6
gfx define field deformed_principal_strain_vector3 composite deformed_principal_strain_vectors.7 deformed_principal_strain_vectors.8 deformed_principal_strain_vectors.9
# since above vectors have the stretch as their magnitude, normalize them:
gfx define field norm_def_principal_strain_vector1 normalise field deformed_principal_strain_vector1
gfx define field norm_def_principal_strain_vector2 normalise field deformed_principal_strain_vector2
gfx define field norm_def_principal_strain_vector3 normalise field deformed_principal_strain_vector3

gfx define field principal_strain1 composite principal_strains.1
gfx define field principal_strain2 composite principal_strains.2
gfx define field principal_strain3 composite principal_strains.3


#Create some materials
gfx create material bluey ambient 0 0.25 0.5 diffuse 0 0.4 1 emission 0 0 0 specular 0.5 0.5 0.5 alpha 1 shininess 0.3
gfx create material copper ambient 1 0.2 0 diffuse 0.6 0.3 0 emission 0 0 0 specular 0.7 0.7 0.5 alpha 1 shininess 0.3
gfx create material gold ambient 1 0.4 0 diffuse 1 0.7 0 emission 0 0 0 specular 0.5 0.5 0.5 alpha 1 shininess 0.8
gfx create material silver ambient 0.4 0.4 0.4 diffuse 0.7 0.7 0.7 emission 0 0 0 specular 0.7 0.7 0.7 alpha 1 shininess 0.6

#Create a strain spectrum
gfx create spectrum strain
gfx modify spectrum strain clear overwrite_colour
gfx modify spectrum strain linear reverse range -1 0 extend_below red colour_range 1 1 ambient diffuse component 1
gfx modify spectrum strain linear reverse range 0 1 extend_above blue colour_range 1 1 ambient diffuse component 1
gfx modify spectrum strain linear reverse range 0 1 extend_above green colour_range 0.5 0.5 ambient diffuse component 1
gfx modify spectrum default linear reverse range 0.104382 0.737732 extend_above extend_below rainbow colour_range 0 1 component 1;


#Show model and strain vectors
gfx modify g_element "/" general clear;
gfx modify g_element "/" lines domain_mesh1d coordinate coordinates exterior tessellation default LOCAL line line_base_size 0 select_on material gold selected_material default_selected render_shaded;
gfx modify g_element "/" points domain_mesh_highest_dimension coordinate coordinates tessellation default_points LOCAL glyph cone REPEAT_MODE_MIRROR size "0*0.1*0.1" offset 0,0,0 font default orientation norm_def_principal_strain_vector1 variable_scale principal_strain1 scale_factors "1*0*0" cell_centres select_on material bluey data principal_strain1 spectrum strain selected_material default_selected render_shaded;
gfx modify g_element "/" points domain_mesh_highest_dimension coordinate coordinates tessellation default_points LOCAL glyph cone REPEAT_MODE_MIRROR size "0*0.1*0.1" offset 0,0,0 font default orientation norm_def_principal_strain_vector2 variable_scale principal_strain2 scale_factors "1*0*0" cell_centres select_on material bluey data principal_strain2 spectrum strain selected_material default_selected render_shaded;
gfx modify g_element "/" points domain_mesh_highest_dimension coordinate coordinates tessellation default_points LOCAL glyph cone REPEAT_MODE_MIRROR size "0*0.1*0.1" offset 0,0,0 font default orientation norm_def_principal_strain_vector3 variable_scale principal_strain3 scale_factors "1*0*0" cell_centres select_on material bluey data principal_strain3 spectrum strain selected_material default_selected render_shaded;
gfx modify g_element "/" lines domain_mesh1d coordinate deform exterior tessellation default LOCAL line line_base_size 0 select_on material default data principal_strain1 spectrum default selected_material default_selected render_shaded;




gfx cre win
