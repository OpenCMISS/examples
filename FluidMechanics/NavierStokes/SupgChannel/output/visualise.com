$numProcs=0;
$tmax=100;
$tinc=10;

for ($t=0; $t<=$tmax; $t+=$tinc){
  for ($proc=$numProcs; $proc>=0; $proc=$proc-1){
    $time =  sprintf("%04d", $t);
    $node = "TIME_STEP_".$time.".exnode";
    gfx read node $node time $t;
  }
}

gfx read elem TIME_STEP_0000.exelem;
gfx def faces egroup "OpenCMISS";

gfx def field velocity component U.1 U.2;
gfx def field pressure component U.3;
gfx def field velmag magnitude field velocity;

gfx modify g_element "OpenCMISS" node_points as node_spheres glyph arrow_solid general size "0.01*0.01*0.01" centre 0,0,0 font default orientation velocity scale_factors "0.05*0.05*0.05" select_on material default; 
#gfx modify g_element "OpenCMISS" node_points as node_spheres glyph arrow_solid general size "0.1*0.1*0.1" centre 0,0,0 font default orientation velocity scale_factors "0.5*0.1*0.1" select_on material default data velmag spectrum default selected_material default_selected;

gfx modify spectrum default clear overwrite_colour;
gfx modify spectrum default linear reverse range 0.0 1.5 extend_above extend_below rainbow colour_range 0 1 component 1;
#gfx modify spectrum default linear reverse range 0.0 1.0 extend_above extend_below banded number_of_bands 10 band_ratio 0.06 component 1;
gfx modify g_element "OpenCMISS" lines select_on material default selected_material default_selected;
gfx modify g_element "OpenCMISS" surfaces select_on material default data velmag spectrum default selected_material default_selected render_shaded;
gfx modify spectrum default autorange;

gfx modify g_element "OpenCMISS" lines delete;
gfx mod g_e "OpenCMISS" lines exterior;

gfx create window 1;


