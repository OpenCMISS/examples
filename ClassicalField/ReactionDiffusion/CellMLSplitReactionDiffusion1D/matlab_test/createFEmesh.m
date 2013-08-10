function [nodes elems] = createFEmesh(numEl,numNod,xL)

%just having a simple function to create a simple mesh of a 1-D rod of lenght xL

dx = xL/numEl;

%store node num in first column, and coord in second
x = 0;
for node = 1:numNod
	nodes(node,:) = [node x+dx*(node-1)];
end

%elements stored as a two-d array with each row being an element, and the columns stating the node numbers
%I have just come up with the following element entries based on simple observation of node and element numbering patterns on a 1D rod.
elems = [[1:1:numEl]' [2:1:numNod]'];


