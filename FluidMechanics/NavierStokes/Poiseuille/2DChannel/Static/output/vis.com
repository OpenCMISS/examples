gfx read node STATICSOLUTION.exnode;
gfx read elem STATICSOLUTION.exelem;
gfx def faces egroup OpenCMISS;

#gfx read node AnalyticNavierStokes.part0.exnode;
#gfx read elem AnalyticNavierStokes.part0.exelem;
#gfx def faces egroup Region_2;

gfx create window 1;

gfx def field ana_velocity component exact.1 exact.2;
gfx def field ana_velmag magnitude field ana_velocity;

gfx def field velocity component general.1 general.2;
gfx def field pressure component general.3;
gfx def field velmag magnitude field velocity;

gfx def field neg1 constant -1;
gfx def field neg_ana multiply_components fields neg1 ana_velocity;
gfx def field error1 add fields neg_ana.1 velocity.1;
gfx def field error2 add fields neg_ana.2 velocity.2;
gfx def field numerical_error component error1 error2;
gfx def field err_mag magnitude field numerical_error;

gfx modify g_element OpenCMISS node_points as node_spheres glyph arrow_solid general size "0.01*0.01*0.01" centre 0,0,0 font default orientation velocity scale_factors "0.00025*0.00025*0.00025" select_on material default; 

gfx modify spectrum default clear overwrite_colour;
gfx modify spectrum default linear reverse range 0.0 1.0 extend_above extend_below rainbow colour_range 0 1 component 1;
#gfx modify spectrum default linear reverse range 0.0 1.0 extend_above extend_below banded number_of_bands 10 band_ratio 0.06 component 1;
gfx modify g_element "OpenCMISS" lines select_on material default selected_material default_selected;
gfx modify g_element "OpenCMISS" surfaces select_on material default data velmag spectrum default selected_material default_selected render_shaded;
#gfx modify spectrum default autorange;

gfx modify window 1 layout 2D;
gfx modify g_element OpenCMISS lines delete;
gfx mod g_e OpenCMISS lines exterior;
