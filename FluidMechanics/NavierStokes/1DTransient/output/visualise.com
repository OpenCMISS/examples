#Read in the sequence of nodal positions.
gfx read node ./output/MainTime_0.part0.exnode time 1;
gfx read node ./output/MainTime_1.part0.exnode time 2;
gfx read node ./output/MainTime_2.part0.exnode time 3;
gfx read node ./output/MainTime_3.part0.exnode time 4;
gfx read node ./output/MainTime_4.part0.exnode time 5;
gfx read node ./output/MainTime_5.part0.exnode time 6;
gfx read node ./output/MainTime_6.part0.exnode time 7;
gfx read node ./output/MainTime_7.part0.exnode time 8;
gfx read node ./output/MainTime_8.part0.exnode time 9;
gfx read node ./output/MainTime_9.part0.exnode time 10;
gfx read node ./output/MainTime_10.part0.exnode time 11;
gfx read node ./output/MainTime_11.part0.exnode time 12;
gfx read node ./output/MainTime_12.part0.exnode time 13;
gfx read node ./output/MainTime_13.part0.exnode time 14;
gfx read node ./output/MainTime_14.part0.exnode time 15;
gfx read node ./output/MainTime_15.part0.exnode time 16;
gfx read node ./output/MainTime_16.part0.exnode time 17;
gfx read node ./output/MainTime_17.part0.exnode time 18;
gfx read node ./output/MainTime_18.part0.exnode time 19;
gfx read node ./output/MainTime_19.part0.exnode time 20;
gfx read node ./output/MainTime_20.part0.exnode time 21;
gfx read node ./output/MainTime_21.part0.exnode time 22;
gfx read node ./output/MainTime_22.part0.exnode time 23;
gfx read node ./output/MainTime_23.part0.exnode time 24;
gfx read node ./output/MainTime_24.part0.exnode time 25;
gfx read node ./output/MainTime_25.part0.exnode time 26;
gfx read node ./output/MainTime_26.part0.exnode time 27;
gfx read node ./output/MainTime_27.part0.exnode time 28;
gfx read node ./output/MainTime_28.part0.exnode time 29;
gfx read node ./output/MainTime_29.part0.exnode time 30;
#Read in the element description
gfx read elem ./output/MainTime_0.part0.exelem;



gfx define field Coordinates.x component Coordinates.x 
gfx define field Coordinates.y component Coordinates.y

gfx define field General.version_1 node_value fe_field General value version 1 
gfx define field General.version_2 node_value fe_field General value version 2
gfx define field General.version_3 node_value fe_field General value version 3

gfx modify g_element OpenCMISS general circle_discretization 12 
gfx modify g_element OpenCMISS cylinders constant_radius 0.01

gfx define field vector_field coord rectangular_cartesian component General.1 General.2

gfx cre spectrum Flow
gfx modify spectrum Flow linear reverse range -1e-05 1e-06 extend_above extend_below rainbow colour_range 0 1 component 1;
gfx cre spectrum Area
gfx modify spectrum Area linear reverse range 1.0e-06 15.0e-06 extend_above extend_below rainbow colour_range 0 1 component 1;

gfx modify g_element OpenCMISS node_points label General 

gfx cre win

#Set the timekeeper playing
gfx timekeeper default set 1.0;
gfx timekeeper default play;
gfx create time_editor



gfx edit scene
#gfx edit spectrum


