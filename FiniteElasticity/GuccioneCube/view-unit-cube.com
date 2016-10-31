if (defined $fibres) {
	my $fn = "unit-cube.part0.exnode";
	my $i = 0;
	foreach my $f ("524", "525", "526", "527") {
		my $time = $i++;
		print "time = $time; f = $f/$fn\n";
		gfx read node "$f/$fn" time $time;
	}
	gfx read elem "524/unit-cube.part0.exelem";
} else {
	gfx read node "unit-cube.part0.exnode";
	gfx read elem "unit-cube.part0.exelem";
}

gfx define faces;

$f = "del U_del n";
gfx define field reactions composite "$f.1" "$f.2" "$f.3";

gfx define tessellation eps minimum_divisions "1*2" refinement_factors "1" circle_divisions 12;

gfx modify g_element "/" general clear;
gfx modify g_element "/" points domain_nodes coordinate Geometry tessellation default_points LOCAL glyph sphere size "0.01*0.01*0.01" offset 0,0,0 font default select_on material default selected_material default_selected render_shaded;
gfx modify g_element "/" points domain_nodes coordinate DeformedGeometry tessellation default_points LOCAL glyph sphere size "0.01*0.01*0.01" offset 0,0,0 font default select_on material gold selected_material default_selected render_shaded;
gfx modify g_element "/" lines domain_mesh1d coordinate Geometry tessellation default LOCAL circle_extrusion line_base_size 0.01 select_on material default selected_material default_selected render_shaded;
gfx modify g_element "/" lines domain_mesh1d coordinate DeformedGeometry tessellation default LOCAL circle_extrusion line_base_size 0.01 select_on material gold selected_material default_selected render_shaded;
gfx modify g_element "/" points domain_nodes coordinate DeformedGeometry tessellation default_points LOCAL glyph axes_xyz size "0.2*0.2*0.2" offset 0,0,0 font default select_on material gold selected_material default_selected render_shaded;
gfx modify g_element "/" points domain_nodes coordinate DeformedGeometry tessellation default_points LOCAL glyph arrow_solid size "0.2*0.2*0.2" offset 0,0,0 font default orientation reactions scale_factors "0*0*0" select_on material red selected_material default_selected render_shaded;
gfx modify g_element "/" points domain_mesh_highest_dimension coordinate DeformedGeometry tessellation eps LOCAL glyph cylinder_solid size "0.5*0.02*0.02" offset -0.3,0,0 font default orientation Fibre scale_factors "0*0*0" cell_centres select_on material default selected_material default_selected render_shaded;
gfx modify g_element "/" points domain_mesh_highest_dimension coordinate DeformedGeometry tessellation eps LOCAL glyph sheet size "0.5*0.2*0.2" offset -0.3,0,0 font default orientation Fibre scale_factors "0*0*0" cell_centres select_on material blue selected_material default_selected render_shaded;


gfx create window 1 double_buffer;
gfx modify window 1 image scene "/" filter default infinite_viewer_lighting two_sided_lighting;
gfx modify window 1 image add_light default;
gfx modify window 1 image add_light default_ambient;
gfx modify window 1 layout simple ortho_axes z -y eye_spacing 0.25 width 1408 height 898;
gfx modify window 1 set current_pane 1;
gfx modify window 1 background colour 0 0 0 texture none;
gfx modify window 1 view perspective eye_point 1.60026 -2.92777 0.628807 interest_point 0.5 0.5 0.5 up_vector -0.00543241 0.0358092 0.999344 view_angle 40 near_clipping_plane 0.245779 far_clipping_plane 7.23963 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;
gfx modify window 1 set transform_tool current_pane 1 std_view_angle 40 normal_lines no_antialias depth_of_field 0.0 fast_transparency blend_normal;
