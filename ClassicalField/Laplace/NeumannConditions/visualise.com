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

# Open window to view solution
gfx create window 1
