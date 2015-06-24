for ($t=0; $t<=$n; $t+=1){
  $num =  sprintf("%04d", $t);
  $node = "TIME_STEP_".$num.".exnode";
  gfx read node $node time $t;
}

#gfx read node TIME_STEP_0000.exnode;
gfx read elem TIME_STEP_0000.exelem;
gfx def faces egroup OpenCMISS;

gfx create window 1;

gfx def field analytic_velocity component exact.1 exact.2;
gfx def field analytic_pressure component exact.3;
gfx def field analytic_velmag magnitude field analytic_velocity;

gfx def field velocity component general.1 general.2;
gfx def field pressure component general.3;
gfx def field velmag magnitude field velocity;

gfx def field velocity_error component error.1 error.2;
gfx def field pressure_error component error.3;

gfx modify g_element OpenCMISS node_points as node_spheres glyph arrow_solid general size "0.01*0.01*0.01" centre 0,0,0 font default orientation velocity scale_factors "0.05*0.005*0.005" select_on material default; 

gfx modify spectrum default clear overwrite_colour;
gfx modify spectrum default linear reverse range 0.0 1.0 extend_above extend_below rainbow colour_range 0 1 component 1;
#gfx modify spectrum default linear reverse range 0.0 1.0 extend_above extend_below banded number_of_bands 10 band_ratio 0.06 component 1;
gfx modify g_element "OpenCMISS" lines select_on material default selected_material default_selected;
gfx modify g_element "OpenCMISS" surfaces select_on material default data velmag spectrum default selected_material default_selected render_shaded;
#gfx modify spectrum default autorange;

gfx modify window 1 layout 2D;
gfx modify g_element OpenCMISS lines delete;
gfx mod g_e OpenCMISS lines exterior;
