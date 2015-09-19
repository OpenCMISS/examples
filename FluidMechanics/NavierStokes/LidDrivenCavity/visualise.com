$numProcs=1;

for ($proc=0; $proc<$numProcs; $proc=$proc+1){
    $node = "./output/Re2500Elem20x20_GFEM/LidDrivenCavity.part".$proc.".exnode";
    gfx read node $node time 0;
}

for ($proc=0; $proc<$numProcs; $proc=$proc+1){
    $node = "./output/Re2500Elem20x20_RBS/LidDrivenCavity.part".$proc.".exnode";
    gfx read node $node time 1;
}

for ($proc=0; $proc<$numProcs; $proc=$proc+1){
    $elem = "./output/Re2500Elem20x20_RBS/LidDrivenCavity.part".$proc.".exelem";
    gfx read elem $elem;
}

gfx def faces egroup "Cavity";

gfx def field velocity component U.1 U.2;
gfx def field pressure component U.3;
gfx def field velmag magnitude field velocity;

gfx modify g_element "Cavity" node_points as node_spheres glyph arrow_solid general size "0.01*0.01*0.01" centre 0,0,0 font default orientation velocity scale_factors "0.04*0.03*0.03" select_on material default; 

gfx modify spectrum default clear overwrite_colour;
gfx modify spectrum default linear reverse range 0.0 1.0 extend_above extend_below rainbow colour_range 0 1 component 1;
#gfx modify spectrum default linear reverse range 0.0 1.0 extend_above extend_below banded number_of_bands 10 band_ratio 0.06 component 1;
gfx modify g_element "Cavity" lines select_on material default selected_material default_selected;
gfx modify g_element "Cavity" surfaces select_on material default data velmag spectrum default selected_material default_selected render_shaded;
#gfx modify spectrum default autorange;

gfx create window 1;

gfx modify window 1 layout 2D;
gfx modify g_element Cavity lines delete;
gfx mod g_e "Cavity" lines exterior;



