#Read in solution
gfx read node CantileverGrowth.part0.exnode
gfx read element CantileverGrowth.part0.exelem

#Update lines and faces
gfx define faces egroup "Region"

# define deformed geometry and pressure
gfx define field "deformed" component Dependent.1 Dependent.2 Dependent.3
gfx define field "pressure" component Dependent.4

#Create windows
gfx create window 1
gfx modify window 1 set antialias 2
gfx modify window 1 view parallel eye_point 20 -200 20 interest_point 20 20 20 up_vector 0 0 1 view_angle 40 near_clipping_plane 1.5 far_clipping_plane 700 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1
gfx modify window 1 view interest_point 1.0,0.5,0.0 eye_point 1.0,0.5,5.0 up_vector 0.0,1.0,0.0

#Display deformed geometry
gfx modify g_element "Region" lines coordinate deformed select_on material red 
#gfx modify g_element "Region" surfaces coordinate deformed select_on material tissue 
#gfx modify g_element "Region" node_points coordinate deformed glyph sphere General size "2*2*2" centre 0,0,0 font default select_on material default

#Display undeformed lines
gfx modify g_element "Region" lines select_on material green 
gfx modify g_element "/" point  glyph axes general size "2*2*2" centre 0,0,0 font default select_on material default

gfx 
