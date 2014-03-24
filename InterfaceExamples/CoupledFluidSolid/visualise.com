#Read in the sequence of nodal positions.
for($i=0; $i < 26; $i=$i+1)
  {
	 $filename = sprintf("./Solid%05d.part0.exnode", $i);
	 print "Reading $filename time $i\n";
	 gfx read node "$filename" time $i;
	 $filename = sprintf("./Fluid%05d.part0.exnode", $i);
	 print "Reading $filename time $i\n";
	 gfx read node "$filename" time $i;
  }
for($i=50; $i < 1500; $i=$i+25)
  {
	 $filename = sprintf("./Solid%05d.part0.exnode", $i);
	 print "Reading $filename time $i\n";
	 gfx read node "$filename" time $i;
	 $filename = sprintf("./Fluid%05d.part0.exnode", $i);
	 print "Reading $filename time $i\n";
	 gfx read node "$filename" time $i;
  }
for($i=1500; $i < 2001; $i=$i+1)
  {
	 $filename = sprintf("./Solid%05d.part0.exnode", $i);
	 print "Reading $filename time $i\n";
	 gfx read node "$filename" time $i;
	 $filename = sprintf("./Fluid%05d.part0.exnode", $i);
	 print "Reading $filename time $i\n";
	 gfx read node "$filename" time $i;
  }
gfx read elem ./Solid00000.part0.exelem;

gfx define field "displacements" component SolidDF.1 SolidDF.2
gfx define field "hydro_pressure" component SolidDF.3
gfx define field "velocities" component FluidDF.1 FluidDF.2
gfx define field "pressure" component FluidDF.3

gfx define field "magnitude_velocity" magnitude field velocities

gfx create spectrum flow
gfx modify spectrum flow clear overwrite_colour
gfx modify spectrum flow linear reverse range 0.0 0.08 extend_above extend_below rainbow colour_range 0 1

gfx define faces egroup SolidRegion
gfx modify g_element SolidRegion lines coordinate SolidGF select_on material green selected_material default_selected
gfx modify g_element SolidRegion lines coordinate displacements select_on material red selected_material default_selected

gfx modify g_element FluidRegion node_points coordinate FluidGF glyph arrow_solid general size "0.05*0.05*0.05" centre 0,0,0 font default orientation velocities select_on material default data magnitude_velocity scale_factors "0.3*0.3*0.3" selected_material default_selected spectrum flow


gfx create window 1 double_buffer;
gfx modify window 1 image scene default light_model default;
gfx modify window 1 image add_light default;
gfx modify window 1 layout 2d ortho_axes z -y eye_spacing 0.25 width 620 height 572;
gfx modify window 1 set current_pane 1;
gfx modify window 1 background colour 0 0 0 texture none;
gfx modify window 1 view parallel eye_point 0.9 1.2 6.11969 interest_point 0.9 1.2 0 up_vector 3.03565e-17 1 7.79185e-09 view_angle 35.1398 near_clipping_plane 0.0611969 far_clipping_plane 21.8697 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;
gfx modify window 1 set transform_tool current_pane 1 std_view_angle 40 normal_lines antialias 8 depth_of_field 0.0 fast_transparency blend_normal;

gfx edit scene
gfx edit spectrum

gfx timekeeper default play speed 100 skip;
gfx create time_editor
