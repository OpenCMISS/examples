$w=490;                        # width of the graphical window
$h=945;                        # height of the graphical window

# --------------- Reading the arteries of the upper body --------------

gfx create material a_colour ambient 1 0.1 0.1 diffuse 1 0.1 0.1;
gfx create material v_colour ambient 0.1 0.1 1 diffuse 0.1 0.1 1;

#Read in the sequence of nodal positions.
for $i (0..10000)
  {
     $filename = sprintf("./output/MainTime_%01d.part0.exnode", $i);
     print "Reading $filename time $i\n";
     gfx read node "$filename" time $i;
#     $filename = sprintf("./output/MainTime_%01d.part1.exnode", $i);
#     print "Reading $filename time $i\n";
#     gfx read node "$filename" time $i;
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

gfx modify g_element OpenCMISS general circle_discretization 12 

gfx define field vector_field coord rectangular_cartesian component General.1 General.2

gfx cre spectrum Flow
gfx modify spectrum Flow linear reverse range 0.0 1.5 extend_above extend_below rainbow colour_range 0 1 component 1;
gfx cre spectrum Pressure
gfx modify spectrum Pressure linear reverse range 0.0 30.0 extend_above extend_below rainbow colour_range 0 1 component 1;
gfx cre spectrum Conc
gfx modify spectrum Conc linear reverse range 0.0 1.0 extend_above extend_below rainbow colour_range 0 1 component 1;

gfx modify g_element OpenCMISS lines data flow spectrum Flow circle_extrusion line_base_size 5;
gfx modify g_element OpenCMISS node_points;

gfx edit scene
#gfx edit spectrum
gfx create window 1;

# ---------- Creation of timekeeper window  ----------------------------

#Set the timekeeper playing
gfx timekeeper default set 1.0;
gfx timekeeper default speed 1;
gfx timekeeper default play;
gfx create time_editor
#jpeg2yuv -f 25 -j %d.jpg -I p | mpeg2enc -o mpegfile.m1v
#gfx print window 1 file coronary_rotate$i.sgi

