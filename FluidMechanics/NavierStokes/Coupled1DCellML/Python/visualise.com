$w=500;                        # width of the graphical window
$h=1000;                        # height of the graphical window

# --------------- Reading the arteries of the upper body --------------

gfx create material a_colour ambient 1 0.1 0.1 diffuse 1 0.1 0.1;

for ($i=0;$i<24000;$i=$i+10) 
#Read in the sequence of nodal positions.
  {
     $filename = sprintf("./output/MainTime_%01d.part0.exnode", $i);
     $time = $i/10
     print "Reading $filename time $time\n";
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
gfx def field velocity divide_components fields flow area

gfx modify g_element ArterialSystem general circle_discretization 12 

gfx define field vector_field coord rectangular_cartesian component General.1 General.2

gfx cre spectrum Flow
gfx modify spectrum Flow linear reverse range 0.0 200.0 extend_above extend_below rainbow colour_range 0 1 component 1;
gfx cre spectrum Pressure
gfx modify spectrum Pressure linear reverse range 0.0 30.0 extend_above extend_below rainbow colour_range 0 1 component 1;
gfx cre spectrum Velocity
gfx modify spectrum Velocity linear reverse range -0.25 1.25 extend_above extend_below rainbow colour_range 0 1 component 1;
gfx cre spectrum Conc
gfx modify spectrum Conc linear reverse range 0.0 1.0 extend_above extend_below rainbow colour_range 0 1 component 1;

#gfx modify g_element ArterialSystem node_points label cmiss_number
gfx modify g_element ArterialSystem cylinders constant_radius 1.0 data flow spectrum Flow radius_scalar area scale_factor 0.05 tessellation default

# ---------------------- Reading the brain data -----------------------

gfx read nodes example Input/stree/brain;
gfx read elements example Input/stree/brain;

#gfx modify g_element brain surfaces select_on material yellow selected_material default_selected;    
gfx modify g_element brain lines select_on material yellow selected_material default_selected;

# ---------------------- Reading the kidney data ----------------------

gfx read nodes example Input/stree/lkidney;
gfx read elements example Input/stree/lkidney;
gfx read nodes example Input/stree/rkidney;
gfx read elements example Input/stree/rkidney;

#gfx modify g_element lkidney surfaces select_on material gold selected_material default_selected;    
gfx modify g_element lkidney lines select_on material gold selected_material default_selected;
#gfx modify g_element rkidney surfaces select_on material gold selected_material default_selected;    
gfx modify g_element rkidney lines select_on material gold selected_material default_selected;

# ---------------------- Reading the liver data -----------------------

gfx read nodes example Input/stree/liver;
gfx read elements example Input/stree/liver;

#gfx modify g_element liver surfaces select_on material red selected_material default_selected;    
gfx modify g_element liver lines select_on material red selected_material default_selected;

# ---------------------- Reading the spleen data -----------------------

gfx read nodes example Input/stree/spleen;
gfx read elements example Input/stree/spleen;

#gfx modify g_element spleen surfaces select_on material green selected_material default_selected;    
gfx modify g_element spleen lines select_on material green selected_material default_selected;

# ---------------------- Reading the stomach data -----------------------

gfx read nodes example Input/stree/stomach;
gfx read elements example Input/stree/stomach;

#gfx modify g_element stomach surfaces select_on material silver selected_material default_selected;    
gfx modify g_element stomach lines select_on material silver selected_material default_selected;

# ---------------------- Reading the hand data -----------------------

gfx read nodes example Input/stree/rhand;
gfx read elements example Input/stree/rhand;
gfx read nodes example Input/stree/lhand;
gfx read elements example Input/stree/lhand;

#gfx modify g_element lhand surfaces select_on material tissue selected_material default_selected;    
gfx modify g_element lhand lines select_on material tissue selected_material default_selected;
#gfx modify g_element rhand surfaces select_on material tissue selected_material default_selected;    
gfx modify g_element rhand lines select_on material tissue selected_material default_selected;

# ---------------------- Reading the foot data -----------------------

gfx read nodes example Input/stree/rfoot;
gfx read elements example Input/stree/rfoot;
gfx read nodes example Input/stree/lfoot;
gfx read elements example Input/stree/lfoot;

#gfx modify g_element lfoot surfaces select_on material tissue selected_material default_selected;    
gfx modify g_element lfoot lines select_on material tissue selected_material default_selected;
#gfx modify g_element rfoot surfaces select_on material tissue selected_material default_selected;    
gfx modify g_element rfoot lines select_on material tissue selected_material default_selected;

# ---------- Creation of graphical window of appropraite size ---------

gfx edit scene
#gfx edit spectrum
gfx create window 1;
gfx modify window 1 layout width $w height $h;
#gfx modify window 1 view interest_point 0 -1000 0;
gfx modify window 1 view eye_point -1100 1500 300;
#gfx modify window 1 view up_vector -0.207985 0.264624 0.941656;
#gfx modify window 1 view view_angle 28.0;


#Set the timekeeper playing
gfx timekeeper default set 1.0;
gfx timekeeper default speed 1;
gfx timekeeper default play;
gfx create time_editor
#jpeg2yuv -f 25 -j %d.jpg -I p | mpeg2enc -o mpegfile.m1v
#gfx print window 1 file coronary_rotate$i.sgi

