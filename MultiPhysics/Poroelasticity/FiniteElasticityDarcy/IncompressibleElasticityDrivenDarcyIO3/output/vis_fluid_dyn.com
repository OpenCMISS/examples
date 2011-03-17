#Read in the sequence of nodal positions.
for $i (1..2)
  {
	 $filename = sprintf("./output/TIME_STEP_%04d.exnode", $i);
#	 $filename = sprintf("./output/T_001_SUB_%04d.exnode", $i);
	 print "Reading $filename time $i\n";
	 gfx read node "$filename" time $i;
  }

#Read in the element description
gfx read elem ./output/TIME_STEP_0001.exelem;
gfx create window 1

gfx cre spectrum flow
gfx modify spectrum flow clear overwrite_colour;
gfx modify spectrum flow linear reverse range 0 1 extend_above extend_below rainbow colour_range 0 1 component 1;

gfx def field vector_field coord rectangular_cartesian component general.1 general.2 general.3

gfx def field vec_mag mag field vector_field

gfx cre spectrum mass_increase
gfx modify spectrum mass_increase linear reverse range -0.0002 0.0002 extend_above extend_below rainbow colour_range 0 1 component 4;

#gfx cre spectrum material_spec
#gfx modify spectrum material_spec linear reverse range 0.25 0.35 extend_above extend_below rainbow colour_range 0 1 component material 1;

gfx cre spectrum source_spectrum
gfx modify spectrum source_spectrum linear reverse range 4.199 4.201 extend_above extend_below rainbow colour_range 0 1 component source 4;

gfx modify g_element OpenCMISS node_points glyph arrow_solid general size "1.5*1.5*1.5" centre 0,0,0 select_on material default selected_material default_selected data vec_mag orientation vector_field scale_factors "1.0*1.0*1.0" spectrum flow

gfx modify window 1 background colour 1 1 1

gfx define faces egroup OpenCMISS
gfx modify g_element OpenCMISS surfaces select_on material default selected_material default_selected data general spectrum mass_increase

#gfx define faces egroup OpenCMISS
#gfx modify g_element OpenCMISS surfaces select_on material default selected_material default_selected data material spectrum material

#Set the timekeeper playing
gfx timekeeper default play speed 5 skip;
gfx create time_editor

gfx edit scene
gfx edit spectrum

#gfx mod win 1 view near_clipping 6.5

