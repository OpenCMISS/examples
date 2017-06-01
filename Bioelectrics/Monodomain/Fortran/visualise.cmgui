gfx read node MonodomainExample.part0.exnode;
gfx read element MonodomainExample.part0.exelem; 
gfx read node Time_2_1.part0.exnode
#gfx read node Time_2_2.part0.exnode time 2
#gfx read node Time_2_3.part0.exnode time 3
#gfx read node Time_2_4.part0.exnode time 4
#gfx read node Time_2_5.part0.exnode time 5
#gfx read node Time_2_6.part0.exnode time 6
#gfx read node Time_2_7.part0.exnode time 7
#gfx read node Time_2_8.part0.exnode time 8
#gfx read node Time_2_9.part0.exnode time 9
#gfx read node Time_2_10.part0.exnode time 10
#gfx read node Time_2_11.part0.exnode time 11
#gfx read node Time_2_12.part0.exnode time 12
#gfx read node Time_2_13.part0.exnode time 13
#gfx read node Time_2_14.part0.exnode time 14
#gfx read node Time_2_15.part0.exnode time 15
#gfx read node Time_2_16.part0.exnode time 16
#gfx read node Time_2_17.part0.exnode time 17
#gfx read node Time_2_18.part0.exnode time 18
#gfx read node Time_2_19.part0.exnode time 19
#gfx read node Time_2_20.part0.exnode time 20
#gfx read node Time_2_21.part0.exnode time 21
#gfx read node Time_2_22.part0.exnode time 22
#gfx read node Time_2_23.part0.exnode time 23
#gfx read node Time_2_24.part0.exnode time 24
#gfx read node Time_2_25.part0.exnode time 25
#gfx read node Time_2_26.part0.exnode time 26
#gfx read node Time_2_27.part0.exnode time 27
#gfx read node Time_2_28.part0.exnode time 28
#gfx read node Time_2_29.part0.exnode time 29
#gfx read node Time_2_30.part0.exnode time 30

#@exnodes=<./Time_*.part*.exnode>;
#$time_index=1;
#foreach $filename (@exnodes) {
#    print "Reading $filename\n";
#    gfx read node "$filename" time $time_index;
#    $time_index++;
#}
gfx define faces egroup Region;
gfx create window 1;
gfx modify window 1 background colour 1.0 1.0 1.0;
gfx modify window 1 view interest_point 0.5,0.5,0.0 eye_point 0.5,0.5,3.0 up_vector 0.0,1.0,0.0;
gfx modify spectrum default clear overwrite_colour;
gfx modify spectrum default linear reverse range -95.0 50.0 extend_above extend_below rainbow colour_range 0 1 component 1;
gfx modify spectrum default linear reverse range -95.0 50.0 banded number_of_bands 10 band_ratio 0.05 component 1;
gfx modify g_element Region lines material black;
gfx modify g_element Region surfaces select_on coordinate Coordinate material default data Vm spectrum default selected_material default_selected render_shaded;
gfx print jpg window 1 file gNa.jpg;
gfx print postscript window 1 file gNa.ps;
