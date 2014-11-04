for ($i = 22000; $i <= 33000; $i+=10)
#Read in the sequence of nodal positions.
  {
     $filename = sprintf("./output/MainTime_%01d.part0.exnode", $i);
     $time = $i/10
     print "Reading $filename time $time\n";
     gfx read node "$filename" time $i;
  }

#Read in the element description
gfx read element ./output/MainTime_0.part0.exelem;
#gfx read element ./output/MainTime_0.part1.exelem;

gfx define field Coordinates.x component Coordinates.x 
gfx define field Coordinates.y component Coordinates.y
gfx define field Coordinates.z component Coordinates.z

gfx define field General.version_1 node_value fe_field General value version 1 
gfx define field General.version_2 node_value fe_field General value version 2
gfx define field General.version_3 node_value fe_field General value version 3

gfx def field flow component General.1
gfx def field area component General.2
gfx def field velocity divide_components fields flow area

gfx def field tommHg constant 7500
gfx def field pressureMmHg multiply_components fields tommHg Pressure

gfx modify g_element OpenCMISS general circle_discretization 12 

gfx define field vector_field coord rectangular_cartesian component General.1 General.2

gfx cre spectrum Flow
gfx modify spectrum Flow linear reverse range -50.0 350.0 extend_above extend_below rainbow colour_range 0 1 component 1;
gfx cre spectrum Pressure
gfx modify spectrum Pressure linear reverse range 0.0 60.0 extend_above extend_below rainbow colour_range 0 1 component 1;
gfx cre spectrum Velocity
gfx modify spectrum Velocity linear reverse range -0.25 1.25 extend_above extend_below rainbow colour_range 0 1 component 1;
gfx cre spectrum Conc
gfx modify spectrum Conc linear reverse range 0.0 1.0 extend_above extend_below rainbow colour_range 0 1 component 1;

#gfx modify g_element OpenCMISS cylinders constant_radius 0.0001 data flow spectrum Flow radius_scalar area  scale_factor 100000
gfx modify g_element OpenCMISS cylinders constant_radius 1.0 data flow spectrum Flow radius_scalar area  scale_factor 0.05
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

