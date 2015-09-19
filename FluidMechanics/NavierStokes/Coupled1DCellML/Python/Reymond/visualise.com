$w=500;                        # width of the graphical window
$h=1000;                        # height of the graphical window

# --------------- Reading the arteries of the upper body --------------

gfx create material a_colour ambient 1 0.1 0.1 diffuse 1 0.1 0.1;

#for ($i=47400;$i<63200;$i=$i+10)
for ($i=0;$i<4000;$i=$i+10)
#for ($i=0;$i<771;$i=$i+1) 
#Read in the sequence of nodal positions.
  {
     $filename = sprintf("./output/MainTime_%01d.part0.exnode", $i);
     $time = $i*0.2
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
gfx modify spectrum Flow log exaggeration 1000 reverse range 0.0 200.0 extend_above extend_below rainbow colour_range 0 1 component 1;
gfx cre spectrum Pressure
gfx modify spectrum Pressure log exaggeration 1000 reverse range 0.0 30.0 extend_above extend_below rainbow colour_range 0 1 component 1;
gfx cre spectrum Velocity
gfx modify spectrum Velocity log exaggeration 1000 reverse range -0.25 1.25 extend_above extend_below rainbow colour_range 0 1 component 1;
gfx cre spectrum Conc
gfx modify spectrum Conc linear reverse range 0.0 1.0 extend_above extend_below rainbow colour_range 0 1 component 1;

gfx modify g_element ArterialSystem node_points label cmiss_number
gfx modify g_element ArterialSystem cylinders constant_radius 1.0 data flow spectrum Flow radius_scalar area scale_factor 0.05 tessellation default

# ---------------------- Reading the brain data -----------------------

gfx read nodes example Input/stree/brain;
gfx read elements example Input/stree/brain;

#gfx modify g_element brain_tree surfaces select_on material yellow selected_material default_selected;    
gfx modify g_element brain_tree lines select_on material yellow selected_material default_selected;



gfx create material brain ambient 1 0.29 0.11 diffuse 0.63 0.58 0.34 emission 0.16 0.11 0.13 specular 0 0 0 alpha 1 shininess 0.0;
gfx create material brain_line ambient 0 0 0 diffuse 0 0 0;

gfx read nodes example Input/Body/brain_organ;
gfx read elements example Input/Body/brain_organ;
gfx modify g_element brain_organ surfaces select_on material brain selected_material default_selected;    
gfx modify g_element brain_organ lines select_on material brain_line selected_material default_selected;

# ---------------------- Reading the heart data -----------------------

gfx create material heart ambient 0 0 0 diffuse 0.79 0.21 0.15 emission 0.39 0 0.0 specular 0.03 0 0;
gfx create material heart_line ambient 0 0 0 diffuse 0 0 0 emission 0.39 0 0.0 specular 0.03 0 0;

gfx read nodes example Input/Body/heart_organ;
gfx read elements example Input/Body/heart_organ;
gfx modify g_element heart_organ surfaces select_on material heart selected_material default_selected;    
gfx modify g_element heart_organ lines select_on material heart_line selected_material default_selected;

# ---------------------- Reading the kidney data ----------------------

gfx read nodes example Input/stree/lkidney;
gfx read elements example Input/stree/lkidney;
gfx read nodes example Input/stree/rkidney;
gfx read elements example Input/stree/rkidney;

#gfx modify g_element lkidney_tree surfaces select_on material gold selected_material default_selected;    
gfx modify g_element lkidney_tree lines select_on material gold selected_material default_selected;
#gfx modify g_element rkidney_tree surfaces select_on material gold selected_material default_selected;    
gfx modify g_element rkidney_tree lines select_on material gold selected_material default_selected;



gfx create material kidney ambient 0.79 0 0 diffuse 0.68 0.21 0.17 emission 0 0 0.05 specular 0.43 0.08 0.16 alpha 1 shininess 0.69;
gfx create material kidney_line ambient 0 0 0 diffuse 0 0 0;

gfx read nodes example Input/Body/kidney_organ;
gfx read elements example Input/Body/kidney_organ;
gfx modify g_element kidney_organ surfaces select_on material kidney selected_material default_selected;    
gfx modify g_element kidney_organ lines select_on material kidney_line selected_material default_selected;

# ---------------------- Reading the liver data -----------------------

gfx read nodes example Input/stree/liver;
gfx read elements example Input/stree/liver;

#gfx modify g_element liver_tree surfaces select_on material red selected_material default_selected;    
gfx modify g_element liver_tree lines select_on material red selected_material default_selected;



gfx create material liver ambient 0.34 0 0 diffuse 0.42 0 0 emission 0.16 0 0 specular 0 0 0 alpha 1 shininess 0.06;
gfx create material liver_line ambient 0 0 0 diffuse 0 0 0;

gfx read nodes example Input/Body/liver_organ;
gfx read elements example Input/Body/liver_organ;
gfx modify g_element liver_organ surfaces select_on material liver selected_material default_selected;    
gfx modify g_element liver_organ lines select_on material liver_line selected_material default_selected;

# ---------------------- Reading the spleen data -----------------------

gfx read nodes example Input/stree/spleen;
gfx read elements example Input/stree/spleen;

#gfx modify g_element spleen_tree surfaces select_on material green selected_material default_selected;    
gfx modify g_element spleen_tree lines select_on material green selected_material default_selected;

# ---------------------- Reading the stomach data -----------------------

gfx read nodes example Input/stree/stomach;
gfx read elements example Input/stree/stomach;

#gfx modify g_element stomach_tree surfaces select_on material silver selected_material default_selected;    
gfx modify g_element stomach_tree lines select_on material silver selected_material default_selected;

# ---------------------- Reading the hand data -----------------------

gfx read nodes example Input/stree/rhand;
gfx read elements example Input/stree/rhand;
gfx read nodes example Input/stree/lhand;
gfx read elements example Input/stree/lhand;

#gfx modify g_element lhand_tree surfaces select_on material tissue selected_material default_selected;    
gfx modify g_element lhand_tree lines select_on material tissue selected_material default_selected;
#gfx modify g_element rhand_tree surfaces select_on material tissue selected_material default_selected;    
gfx modify g_element rhand_tree lines select_on material tissue selected_material default_selected;

# ---------------------- Reading the foot data -----------------------

gfx read nodes example Input/stree/rfoot;
gfx read elements example Input/stree/rfoot;
gfx read nodes example Input/stree/lfoot;
gfx read elements example Input/stree/lfoot;

#gfx modify g_element lfoot_tree surfaces select_on material tissue selected_material default_selected;    
gfx modify g_element lfoot_tree lines select_on material tissue selected_material default_selected;
#gfx modify g_element rfoot_tree surfaces select_on material tissue selected_material default_selected;    
gfx modify g_element rfoot_tree lines select_on material tissue selected_material default_selected;

# ------------ Reading the surface mesh of the entire body ------------

gfx create material skin ambient 1 1 0 diffuse 1 1 1 emission 0.18 0 0 alpha 0.3;

gfx read nodes example Input/Body/surface_organ;
gfx read elements example Input/Body/surface_organ;
gfx modify g_element surface_organ lines select_on invisible;
gfx modify g_element surface_organ surfaces select_on material skin render_shaded;

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

