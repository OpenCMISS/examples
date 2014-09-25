#Read in the sequence of nodal positions.
for $i (0..1000)
  {
	 $filename = sprintf("./output/MainTime_%01d.part0.exnode", $i);
	 print "Reading $filename time $i\n";
	 gfx read node "$filename" time $i;
  }

#Read in the element description
gfx read elem ./output/MainTime_0.part0.exelem;

gfx define field Coordinates.x component Coordinates.x 
gfx define field Coordinates.y component Coordinates.y

gfx define field General.version_1 node_value fe_field General value version 1 
gfx define field General.version_2 node_value fe_field General value version 2
gfx define field General.version_3 node_value fe_field General value version 3

gfx def field flow component General.1
gfx def field area component General.2

gfx modify g_element OpenCMISS general circle_discretization 12 

gfx define field vector_field coord rectangular_cartesian component General.1 General.2

gfx cre spectrum Flow
gfx modify spectrum Flow linear reverse range 2 3 extend_above extend_below rainbow colour_range 0 1 component 1;
gfx cre spectrum Area
gfx modify spectrum Area linear reverse range 0.2 1.2 extend_above extend_below rainbow colour_range 0 1 component 1;

#gfx modify g_element OpenCMISS cylinders radius area data flow spectrum Flow
gfx modify g_element OpenCMISS cylinders constant_radius 0.01

gfx modify g_element OpenCMISS node_points label cmiss_number
gfx edit scene
gfx edit spectrum
gfx cre win

#Set the timekeeper playing
gfx timekeeper default set 1.0;
gfx timekeeper default speed 5;
gfx timekeeper default play;
gfx create time_editor








