#Read in the sequence of nodal positions.
for $i (0..10000)
  {
     $filename = sprintf("./output/MainTime_%01d.part0.exnode", $i);
     print "Reading $filename time $i\n";
     gfx read node "$filename" time $i;
  }

#Read in the element description
gfx read element ./output/MainTime_0.part0.exelem;

gfx define field Coordinates.x component Coordinates.x 
gfx define field Coordinates.y component Coordinates.y
gfx define field Coordinates.z component Coordinates.z

gfx define field General.version_1 node_value fe_field General value version 1 
gfx define field General.version_2 node_value fe_field General value version 2
gfx define field General.version_3 node_value fe_field General value version 3

gfx def field flow component General.1
gfx def field area component General.2

gfx modify g_element OpenCMISS general circle_discretization 12 


gfx define field vector_field coord rectangular_cartesian component General.1 General.2

gfx cre spectrum Flow
gfx modify spectrum Flow linear reverse range -0.1 1.5 extend_above extend_below rainbow colour_range 0 1 component 1;
gfx cre spectrum Pressure
gfx modify spectrum Pressure linear reverse range 0.0 30.0 extend_above extend_below rainbow colour_range 0 1 component 1;

gfx modify g_element OpenCMISS cylinders constant_radius 1.0 data flow spectrum Flow radius_scalar area  scale_factor 1
gfx modify g_element OpenCMISS node_points label cmiss_number

gfx edit scene
#gfx edit spectrum
gfx cre win

#Set the timekeeper playing
gfx timekeeper default set 1.0;
gfx timekeeper default speed 1;
gfx timekeeper default play;
gfx create time_editor
#jpeg2yuv -f 25 -j %d.jpg -I p | mpeg2enc -o mpegfile.m1v
#gfx print window 1 file coronary_rotate$i.sgi

