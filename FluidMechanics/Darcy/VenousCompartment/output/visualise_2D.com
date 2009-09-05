gfx read node './cmgui.exnode'
gfx read elem './cmgui.exelem'
gfx create window 1
gfx mod win 1 layout width 1024 height 1024

gfx edit scene
gfx edit spectrum

gfx cre spectrum flow
gfx modify spectrum flow clear overwrite_colour;
gfx modify spectrum flow linear reverse range 0 1.0 extend_above extend_below rainbow colour_range 0 1 component 1;

gfx def field vector_field coord rectangular_cartesian component general.1 general.2 

gfx def field mag_1 magnitude field vector_field

gfx cre spectrum pressure
gfx modify spectrum pressure linear reverse range -13 13 extend_above extend_below rainbow colour_range 0 1 component 3;

gfx modify g_element OpenCMISS node_points glyph arrow_solid general size "0.17*0.17*0.17" centre 0,0,0 select_on material default selected_material default_selected data mag_1 orientation vector_field scale_factors "0.01*0.01*0.01" spectrum flow

gfx modify window 1 background colour 1 1 1

gfx define faces egroup OpenCMISS
gfx modify g_element OpenCMISS surfaces select_on material default selected_material default_selected data general spectrum pressure
