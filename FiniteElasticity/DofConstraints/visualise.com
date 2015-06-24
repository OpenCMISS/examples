gfx read node "Cantilever.part0.exnode"
gfx read elem "Cantilever.part0.exelem"

gfx create window 1

gfx define faces egroup "Region"
gfx define field deformed coordinate_system rectangular_cartesian composite Dependent.1 Dependent.2 Dependent.3

# View undeformed geometry
gfx modify g_element "Region" lines coordinate Geometry select_on material white

# View deformed geometry
gfx modify g_element "Region" lines coordinate deformed select_on material green
gfx modify g_element "Region" node_points coordinate deformed glyph sphere General size "0.05*0.05*0.05" centre 0,0,0 font default select_on material white

# Add axis
gfx modify g_element "/" point glyph axes_solid general size "10*10*10" centre 0,0,0 font default select_on material default

gfx modify window 1 view up_vector 0 0 1
gfx modify window 1 image view_all
