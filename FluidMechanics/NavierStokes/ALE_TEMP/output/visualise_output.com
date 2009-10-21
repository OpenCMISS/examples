#Read in the sequence of nodal positions.
for $i (1..100)
  {
	 $filename = sprintf("./output/TIME_STEP_%04d.exnode", $i);
         gfx read node "$filename"
	 print "Reading $filename time $i\n";
         gfx print window 1 file "$filename.tiff" format RGB
  }