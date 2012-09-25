gfx create region laplace;
gfx read region LaplaceExample.xml region laplace;

# The input files do not describe faces, so create them here so we can visualise lines.
gfx define face egroup laplace;

gfx modify g_element laplace general clear;
# Visualise element lines
gfx modify g_element laplace lines coordinate laplace.geometric select_on material default selected_material default_selected;
# View surfaces with solution plotted
gfx modify g_element laplace surfaces coordinate laplace.geometric select_on material default data laplace.phi spectrum default selected_material default_selected render_shaded;
gfx modify spectrum default autorange;
# Create axes
gfx modify g_element "/" point  glyph axes_xyz general size "1*1*1" centre 0,0,0 font default select_on material default selected_material default_selected;

# Open window to view solution
gfx create window 1

gfx modify window 1 set antialias 4
