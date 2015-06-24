gfx read node "prolate_spheroid.part0.exnode"
gfx read elem "prolate_spheroid.part0.exelem"

gfx create window 1

gfx define faces egroup "ProlateSpheroid"

# View undeformed geometry
gfx modify g_element "ProlateSpheroid" lines coordinate Geometry select_on material white

# View deformed geometry
gfx modify g_element "ProlateSpheroid" lines coordinate DeformedGeometry select_on material green
gfx modify g_element "ProlateSpheroid" node_points coordinate DeformedGeometry glyph sphere General size "0.05*0.05*0.05" centre 0,0,0 font default select_on material white

# Add axis
gfx modify g_element "/" point  glyph axes_solid general size "10*10*10" centre 0,0,0 font default select_on material default

# Show internal surface
gfx modify g_element "ProlateSpheroid" surfaces coordinate DeformedGeometry exterior face xi3_0 material white

# View fibre orientations
gfx modify g_element "ProlateSpheroid" element_points coordinate DeformedGeometry discretization "3*3*3" use_elements glyph cylinder_solid size "2.0*0.5*0.5" scale_factors "0*0*0" centre "0.5,0.0,0.0" orientation Fibre material gold
gfx modify g_element "ProlateSpheroid" element_points coordinate DeformedGeometry discretization "3*3*3" use_elements glyph sheet size "2.0*1.2*1.2" scale_factors "0*0*0" centre "0.0,0.0,0.0" orientation Fibre material blue

gfx modify window 1 view up_vector 0 0 -1
gfx modify window 1 image view_all

gfx define tessellation default minimum_divisions "1" refinement_factors "16"
