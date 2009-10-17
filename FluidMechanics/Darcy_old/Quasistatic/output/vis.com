
#Read in the sequence of nodal positions.
for $i (1..5)
  {
	 $filename = sprintf("TIME_STEP_%04d.exnode", $i);
	 print "Reading $filename time $i\n";
	 gfx read node "$filename" time $i;
  }

#Read in the element description
gfx read elem TIME_STEP_0001.exelem;

gfx create window 1
#gfx mod win 1 layout width 1024 height 1024

#Set the timekeeper playing
gfx timekeeper default play speed 1 skip;
gfx create time_editor

gfx edit scene
gfx edit spectrum

