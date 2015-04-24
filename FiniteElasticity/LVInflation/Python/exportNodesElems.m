function [] = exportNodesElems(h, Celems, Cnodes, output_tag)

nab = h(1);
ncirc = h(2);
nr = h(3); % should always be 2, this code cannot handle non-bilayer formations. 


% ---------- MESH OUTPUT

fh =fopen(sprintf('%s.exelem',output_tag),'w');
fprintf(fh,'Group name: Region\nShape. Dimension=3\n #Scale factor sets=1\n l.Lagrange*l.Lagrange*l.Lagrange, #Scale factors=8\n #Nodes=8\n #Fields=1\n 1) coordinates, coordinate, rectangular cartesian, #Components=3\n');
for x='xyz'
  fprintf(fh,['  ' x '. l.Lagrange*l.Lagrange*l.Lagrange, no modify, standard node based.\n  #Nodes=8\n']);
  for i=1:8
    fprintf(fh,['  ' int2str(i) '. #Values=1\n     Value indices: 1\n     Scale factor indices: 0\n']);
  end
end

for e=1:size(Celems,1)
  fprintf(fh,'Element: %d 0 0\nNodes: \n%d %d %d %d %d %d %d %d\nScale factors:\n 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0\n',e,Celems(e,:));
end
fclose(fh);

fh=fopen(sprintf('%s.exnode',output_tag),'w');
fprintf(fh,' Group name: Region\n #Fields=2\n 1) coordinates, coordinate, rectangular cartesian, #Components=3\n');
for x='xyz'
  fprintf(fh,['   ' x '.  Value index=%d, #Derivatives=0\n'],x-'x'+1);
end
fprintf(fh,'2) fibers, coordinate, rectangular cartesian, #Components=3\n');
for x='xyz'
  fprintf(fh,'   %s.  Value index=%d, #Derivatives=0\n',x,x-'x'+4);
end

for ni=1:size(Cnodes,1)
  fprintf(fh,'Node:            %d\n  %.8f\n  %.8f\n  %.8f\n  %.8f\n  0.0\n  0.0\n',ni,Cnodes(ni,:),angle(ni));
end
fclose(fh);
fprintf('Wrote exnode, exelem files for linear hexes to ./out/%s.*\n',output_tag);

fh = fopen(sprintf('%s.ipnode',output_tag),'w');
fprintf(fh,' CMISS Version 1.21 ipnode File Version 2\n');
fprintf(fh,' Heading: Region\n\n');
fprintf(fh, ' The number of nodes is [   1]:   %d\n',size(Cnodes,1));
fprintf(fh, ' Number of coordinates [3]: 3\n');
fprintf(fh, ' Do you want prompting for different versions of nj=1 [N]? N\n');
fprintf(fh, ' Do you want prompting for different versions of nj=2 [N]? N\n');
fprintf(fh, ' Do you want prompting for different versions of nj=3 [N]? N\n');
fprintf(fh, ' The number of derivatives for coordinate 1 is [0]: 0\n');
fprintf(fh, ' The number of derivatives for coordinate 2 is [0]: 0\n');
fprintf(fh, ' The number of derivatives for coordinate 3 is [0]: 0\n\n');
for ni = 1:size(Cnodes,1)
    fprintf(fh,' Node number [    %d]:     %d\n', ni, ni);
    if (ni == 0) || (ni == 0)
        fprintf(fh, ' The number of versions for nj=1 is [1]:   %d\n',ncirc*nr);
        for i = 1:ncirc*nr
            fprintf(fh, ' For version number %d:\n',i);
            fprintf(fh, ' The Xj(1) coordinate is [0.00000E+00]:\t%.8f\n',Cnodes(ni,1));
        end
        fprintf(fh, ' The number of versions for nj=2 is [1]:   %d\n',ncirc*nr);
        for i = 1:ncirc*nr
            fprintf(fh, ' For version number %d:\n',i);
            fprintf(fh, ' The Xj(2) coordinate is [0.00000E+00]:\t%.8f\n',Cnodes(ni,2));
        end
        fprintf(fh, ' The number of versions for nj=3 is [1]:   %d\n',ncirc*nr);
        for i = 1:ncirc*nr
            fprintf(fh, ' For version number %d:\n',i);
            fprintf(fh, ' The Xj(3) coordinate is [0.00000E+00]:\t%.8f\n',Cnodes(ni,3));
        end
        fprintf(fh, '\n');
    else
        fprintf(fh,' The Xj(1) coordinate is [0.00000E+00]:\t%.8f\n',Cnodes(ni,1));
        fprintf(fh,' The Xj(2) coordinate is [0.00000E+00]:\t%.8f\n',Cnodes(ni,2));
        fprintf(fh,' The Xj(3) coordinate is [0.00000E+00]:\t%.8f\n\n',Cnodes(ni,3));
    end
end
fclose(fh);

fh = fopen(sprintf('%s.ipelem',output_tag),'w');
fprintf(fh, ' CMISS Version 2.0 ipelem File Version 2\n');
fprintf(fh, ' Heading:\n\n');
fprintf(fh, ' The number of elements is [1]:\t%d\n',size(Celems,1));

for e = 1:size(Celems,1)
    fprintf(fh, '\n Element number [\t%d]:\t%d\n', e, e);
    fprintf(fh, ' The number of geometric Xj-coordinates is [3]:  3\n');
    fprintf(fh, ' The basis function type for geometric variable 1 is [1]:  1\n');
    fprintf(fh, ' The basis function type for geometric variable 2 is [1]:  1\n');
    fprintf(fh, ' The basis function type for geometric variable 3 is [1]:  1\n');
    fprintf(fh, ' Enter the 8 global numbers for basis 1:  %d  %d  %d  %d  %d  %d  %d  %d\n',Celems(e,:));
    if e <= 0
        fprintf(fh, ' The version number for occurrence  1 of node    18, njj=1 is [1]:\t%d\n',e);
        fprintf(fh, ' The version number for occurrence  1 of node    18, njj=2 is [1]:\t%d\n',e);
        fprintf(fh, ' The version number for occurrence  1 of node    18, njj=3 is [1]:\t%d\n',e);
        if e == ncirc*nr
            fprintf(fh, ' The version number for occurrence  2 of node    18, njj=1 is [1]:\t1\n');
            fprintf(fh, ' The version number for occurrence  2 of node    18, njj=2 is [1]:\t1\n');
            fprintf(fh, ' The version number for occurrence  2 of node    18, njj=3 is [1]:\t1\n');
        else
            fprintf(fh, ' The version number for occurrence  2 of node    18, njj=1 is [1]:\t%d\n',e+1);
            fprintf(fh, ' The version number for occurrence  2 of node    18, njj=2 is [1]:\t%d\n',e+1);
            fprintf(fh, ' The version number for occurrence  2 of node    18, njj=3 is [1]:\t%d\n',e+1);
        end
        fprintf(fh, ' The version number for occurrence  1 of node     1, njj=1 is [1]:\t%d\n',e);
        fprintf(fh, ' The version number for occurrence  1 of node     1, njj=2 is [1]:\t%d\n',e);
        fprintf(fh, ' The version number for occurrence  1 of node     1, njj=3 is [1]:\t%d\n',e);
        if e == ncirc*nr
            fprintf(fh, ' The version number for occurrence  2 of node     1, njj=1 is [1]:\t1\n');
            fprintf(fh, ' The version number for occurrence  2 of node     1, njj=2 is [1]:\t1\n');
            fprintf(fh, ' The version number for occurrence  2 of node     1, njj=3 is [1]:\t1\n');
        else
            fprintf(fh, ' The version number for occurrence  2 of node     1, njj=1 is [1]:\t%d\n',e+1);
            fprintf(fh, ' The version number for occurrence  2 of node     1, njj=2 is [1]:\t%d\n',e+1);
            fprintf(fh, ' The version number for occurrence  2 of node     1, njj=3 is [1]:\t%d\n',e+1);
        end
    end
end
fclose(fh);







