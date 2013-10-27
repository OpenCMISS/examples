% generates a 2D mesh for FSI problem: plate in a cross flow
% author: Andreas Hessenthaler (C)
%

clear all
close all

format longe

% decide which mesh will be plotted:
% 0 - none // 1 - linear mesh // 2 - quadratic mesh
showMeshes = 2;

D=0.2;

SolidOriginX = 3;
SolidOriginY = 4*D;
SolidOriginZ = 0;

SolidX = 2;
SolidY = D;1.25E-3;
SolidZ = 10*D;

FluidX = 8;
FluidY = 8*D+SolidY;
FluidZ = 12*D;

SolidElementNumberX = 10;
SolidElementNumberY = 1;
SolidElementNumberZ = 2;

SolidElementSizeX = SolidX / SolidElementNumberX;
SolidElementSizeY = SolidY / SolidElementNumberY;
SolidElementSizeZ = SolidZ / SolidElementNumberZ;

FluidElementSizeX = SolidElementSizeX;
FluidElementSizeY = SolidElementSizeX;
FluidElementSizeZ = SolidElementSizeZ;

% in x = 4 plane!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


SolidGeometryY = kron(ones(1, SolidElementNumberZ*2+1), SolidOriginY:SolidElementSizeY/2:SolidOriginY+SolidY);
SolidGeometryZ = zeros(1, size(SolidGeometryY, 2));
SolidGeometryPressureY = kron(ones(1, SolidElementNumberZ+1), SolidOriginY:SolidElementSizeY:SolidOriginY+SolidY);
SolidGeometryPressureZ = zeros(1, size(SolidGeometryPressureY, 2));
% converts index of SolidGeometryY/Z into global node number
SolidContinuousToGlobal = zeros(1, size(SolidGeometryY, 2));

for i = 3:2:SolidElementNumberZ*2+1
    SolidGeometryPressureZ(i:i+1) = SolidGeometryPressureZ(1:2) + (i-1)/2*SolidElementSizeZ;
    SolidContinuousToGlobal([1, 3, 7, 9]+6*(i-1)/2-6) = [1, 2, 3, 4]+i-3;
end
for i = 4:3:((SolidElementNumberZ+1)*3-1)*2-3
    SolidGeometryZ(i:i+2) = SolidGeometryZ(1:3) + (i-1)/3*SolidElementSizeZ/2;
end
number = SolidContinuousToGlobal(end);
for i = 1:size(SolidGeometryZ,2)
    if(SolidContinuousToGlobal(i)==0)
        number = number + 1;
        SolidContinuousToGlobal(i) = number;
    end
end
SolidElementNodes = zeros(SolidElementNumberY*SolidElementNumberZ, 9);
SolidElementNodesXi1 = SolidElementNodes;
SolidElementNodesXi2 = SolidElementNodesXi1;
SolidInterfaceNodesContinuousToGlobal = zeros(2, size(SolidContinuousToGlobal, 2)-(SolidElementNumberZ*2));
SolidInterfaceNodesGlobalToContinuous = SolidInterfaceNodesContinuousToGlobal(1,:);
xi1 = [0.0, 0.5, 1.0, 0.0, 0.5, 1.0, 0.0, 0.5, 1.0];
xi2 = [0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0];
for i = 1:SolidElementNumberY*SolidElementNumberZ
    SolidElementNodes(i,:) = SolidContinuousToGlobal([1:9]+(i-1)*6);
    SolidElementNodesXi1(i,:) = xi1;
    SolidElementNodesXi2(i,:) = xi2;
end
xi1s = zeros(size(SolidInterfaceNodesContinuousToGlobal, 2));
for i = 1:2*SolidElementNumberZ
    SolidInterfaceNodesContinuousToGlobal(1, 2*i-1:2*i) = SolidContinuousToGlobal([1, 3]+3*(i-1));
    SolidInterfaceNodesGlobalToContinuous(1, 2*i-1:2*i) = [1, 3]+3*(i-1);
end
SolidInterfaceNodesContinuousToGlobal(1, end-2:end) = SolidContinuousToGlobal(end-2:end);
SolidInterfaceNodesGlobalToContinuous(1, end-2:end) = size(SolidContinuousToGlobal,2)-2:size(SolidContinuousToGlobal,2);

% converts the global node number into the index of SolidGeometryY/Z
SolidGlobalToContinuous = zeros(1, size(SolidContinuousToGlobal, 2));
for i = 1:size(SolidContinuousToGlobal, 2)
    SolidGlobalToContinuous(SolidContinuousToGlobal(i)) = i;
end
InterfaceNodes = 1:size(SolidGlobalToContinuous, 2);
InterfaceNodes = [InterfaceNodes(1:2:end-2), InterfaceNodes(end-1), InterfaceNodes(end), InterfaceNodes(end-3:-2:2)];
InterfaceGeometryY = SolidGeometryY(SolidInterfaceNodesGlobalToContinuous);
InterfaceGeometryZ = SolidGeometryZ(SolidInterfaceNodesGlobalToContinuous);

InterfaceElementNodes = 1:3;
InterfaceElementNodesXi1 = [0.0, 0.5, 1.0];
for i = 1:(2*SolidElementNumberZ+SolidElementNumberY)-1
    InterfaceElementNodes(i+1, 1:3) = InterfaceElementNodes(i, 1:3)+2;
end

NumberOfInterfaceNodes = max(max(InterfaceElementNodes));

InterfaceGlobalNodesForGeometry = zeros(1, NumberOfInterfaceNodes);

InterfaceGlobalNodesForGeometry(1:2:NumberOfInterfaceNodes-1) = 1:floor(NumberOfInterfaceNodes/2);
InterfaceGlobalNodesForGeometry(2:2:NumberOfInterfaceNodes-2) = NumberOfInterfaceNodes:-1:floor(NumberOfInterfaceNodes/2)+3;
InterfaceGlobalNodesForGeometry(end-1) = InterfaceGlobalNodesForGeometry(end-2)+1;
InterfaceGlobalNodesForGeometry(end) = InterfaceGlobalNodesForGeometry(end-1)+1;

% plot(SolidGeometryY, SolidGeometryZ, '*')
%plot(InterfaceGeometryY, InterfaceGeometryZ, '*')

FluidGeometryY = kron(ones(1,SolidElementNumberZ*2),[0:SolidElementSizeX/2:SolidOriginY, SolidOriginY+SolidY:SolidElementSizeX/2:FluidY]);
FluidGeometryZ = kron([0:SolidElementSizeZ/2:SolidZ-SolidElementSizeZ/2],ones(1,size([0:SolidElementSizeX/2:SolidOriginY, SolidOriginY+SolidY:SolidElementSizeX/2:FluidY],2)));
numberOfFluidNodesInLowerRow = size([0:SolidElementSizeX/2:SolidOriginY, SolidOriginY+SolidY:SolidElementSizeX/2:FluidY], 2);
temp = size(FluidGeometryY, 2);
FluidGeometryPressureY = kron(ones(1,SolidElementNumberZ),[0:SolidElementSizeX:SolidOriginY, SolidOriginY+SolidY:SolidElementSizeX:FluidY]);
FluidGeometryPressureZ = kron([0:SolidElementSizeZ:SolidZ-SolidElementSizeZ],ones(1,size([0:SolidElementSizeX:SolidOriginY, SolidOriginY+SolidY:SolidElementSizeX:FluidY],2)));
% plot(FluidGeometryY, FluidGeometryZ, '*')
% plot(FluidGeometryPressureY, FluidGeometryPressureZ, '*')
tempPressure = size(FluidGeometryPressureY, 2);

FluidGeometryZ = [FluidGeometryZ, kron(SolidZ:D/2:FluidZ,ones(1,size([0:SolidElementSizeX/2:SolidOriginY, SolidOriginY+SolidElementSizeY/2, SolidOriginY+SolidY:SolidElementSizeX/2:FluidY],2)))];
FluidGeometryY = [FluidGeometryY, kron(ones(1,size(SolidZ:D/2:FluidZ,2)),[0:SolidElementSizeX/2:SolidOriginY, SolidOriginY+SolidElementSizeY/2, SolidOriginY+SolidY:SolidElementSizeX/2:FluidY])];
temp2 = temp + size([0:SolidElementSizeX/2:SolidOriginY, SolidOriginY+SolidElementSizeY/2, SolidOriginY+SolidY:SolidElementSizeX/2:FluidY],2);
% plot(FluidGeometryY, FluidGeometryZ, '*')
FluidGeometryPressureZ = [FluidGeometryPressureZ, kron(SolidZ:D:FluidZ,ones(1,size([0:SolidElementSizeX:SolidOriginY, SolidOriginY+SolidY:SolidElementSizeX:FluidY],2)))];
FluidGeometryPressureY = [FluidGeometryPressureY, kron(ones(1,size(SolidZ:D:FluidZ,2)),[0:SolidElementSizeX:SolidOriginY, SolidOriginY+SolidY:SolidElementSizeX:FluidY])];
% plot(FluidGeometryPressureY, FluidGeometryPressureZ, '*')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% pressure nodes first, then numbers of other nodes!!!

FirstFluidInterfaceNode = size(0:SolidElementSizeX/2:SolidOriginY,2);
NumberOfFluidNodesLowerY = size([0:SolidElementSizeX/2:SolidOriginY, SolidOriginY+SolidY:SolidElementSizeX/2:FluidY],2);

FluidInterfaceNodeArrayIndices(1:2:NumberOfInterfaceNodes-4) = FirstFluidInterfaceNode:NumberOfFluidNodesLowerY:temp;
FluidInterfaceNodeArrayIndices(2:2:NumberOfInterfaceNodes-3) = FirstFluidInterfaceNode+1:NumberOfFluidNodesLowerY:temp;
FluidInterfaceNodeArrayIndices(end+1) = FluidInterfaceNodeArrayIndices(end-1)*2-FluidInterfaceNodeArrayIndices(end-3);
FluidInterfaceNodeArrayIndices(end+1) = FluidInterfaceNodeArrayIndices(end)+1;
FluidInterfaceNodeArrayIndices(end+1) = FluidInterfaceNodeArrayIndices(end)+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% pressure nodes first, then numbers of other nodes!!!


NumberOfFluidNodes = size(FluidGeometryY, 2);
NumberOfFluidPressureNodes = size(FluidGeometryPressureY, 2);
FluidNodes = zeros(1, NumberOfFluidNodes);
index = 1;
for i = 1:NumberOfFluidNodes
    for j = 1:NumberOfFluidPressureNodes
        if(abs(FluidGeometryY(i)-FluidGeometryPressureY(j))<0.001/2 && abs(FluidGeometryZ(i)-FluidGeometryPressureZ(j))<0.001/2)
            if(FluidNodes(i)~=0)
                index = index - 1;
            end
            FluidNodes(i) = j;
            break
        else
            if(FluidNodes(i)==0)
                FluidNodes(i) = NumberOfFluidPressureNodes + index;
                index = index + 1;
            end
        end
    end
end

SolidInterfaceNodesContinuousToGlobal(2,:) = FluidNodes(FluidInterfaceNodeArrayIndices);
% SolidInterfaceNodesContinuousToGlobal & InterfaceGlobalNodesForGeometry
% have same node ordering. Bottom left to top right.


% bottom element row from left end to solid
FluidIndexArray = [1:3, 1+numberOfFluidNodesInLowerRow:3+numberOfFluidNodesInLowerRow, 1+numberOfFluidNodesInLowerRow*2:3+numberOfFluidNodesInLowerRow*2];
k = 1;
j = 1;
m = 1;
check = 0;
FluidElementNumbersOnInterface = 0;
% possibly no need for the breaks.. leave em in for now just to be sure
while (FluidNodes(FluidIndexArray(end)) <= FluidNodes(temp2))
    while (FluidNodes(FluidIndexArray(end))<=SolidInterfaceNodesContinuousToGlobal(2, 1+4*m))
        if(FluidNodes(FluidIndexArray(end)) > FluidNodes(temp2))
            check = 1;
            disp 'break1'
            break
        else
            FluidElementNodes(k,:) = FluidNodes(FluidIndexArray);
            k = k + 1;
            if (m > 1)
                FluidIndexArray = FluidIndexArray + 2;
            else
                j = j + 2;
                FluidIndexArray = [j:j+2, j+numberOfFluidNodesInLowerRow:j+2+numberOfFluidNodesInLowerRow, j+numberOfFluidNodesInLowerRow*2:j+2+numberOfFluidNodesInLowerRow*2];
            end
        end
    end
    if (FluidElementNumbersOnInterface(1)==0)
        FluidElementNumbersOnInterface = [k-1, k];
    else
        FluidElementNumbersOnInterface = [FluidElementNumbersOnInterface, k-1, k];
    end
    % from solid to right end
    if (check)
        disp 'break check'
        break
    end
    if (m > 1)
        if (m == SolidElementNumberZ)
            FluidIndexArray = FluidIndexArray + [1, 1, 1, 1, 1, 1, 2, 2, 2];
        else
            FluidIndexArray = FluidIndexArray + 1;
        end
    else
        j = j + 1;
        FluidIndexArray = [j:j+2, j+numberOfFluidNodesInLowerRow:j+2+numberOfFluidNodesInLowerRow, j+numberOfFluidNodesInLowerRow*2:j+2+numberOfFluidNodesInLowerRow*2];
    end
    while (FluidNodes(FluidIndexArray(end))<=FluidNodes(numberOfFluidNodesInLowerRow*(2*m+1)))
        if(FluidNodes(FluidIndexArray(end)) > FluidNodes(temp2))
            check = 1;
            disp 'break2'
            break
        else
            FluidElementNodes(k,:) = FluidNodes(FluidIndexArray);
            k = k + 1;
            if (m > 1)
                storeIndex = FluidIndexArray(end);
                FluidIndexArray = FluidIndexArray + 2;
            else
                j = j + 2;
                FluidIndexArray = [j:j+2, j+numberOfFluidNodesInLowerRow:j+2+numberOfFluidNodesInLowerRow, j+numberOfFluidNodesInLowerRow*2:j+2+numberOfFluidNodesInLowerRow*2];
            end
        end
    end
    if (check)
        disp 'break check 2'
        break
    end
    m = m + 1;
    if (m > 1)
        FluidIndexArray = [1+numberOfFluidNodesInLowerRow*2*(m-1):3+numberOfFluidNodesInLowerRow*2*(m-1), ...
            1+numberOfFluidNodesInLowerRow*(2*m-1):3+numberOfFluidNodesInLowerRow*(2*m-1), ...
            1+numberOfFluidNodesInLowerRow*2*m:3+numberOfFluidNodesInLowerRow*2*m];
    end
end
% complete FluidElementNumbersOnInterface
FluidElementNumbersOnInterface = [FluidElementNumbersOnInterface, 2*FluidElementNumbersOnInterface(end-1)-FluidElementNumbersOnInterface(end-3)+1];
% correct fluid index array
FluidIndexArray = FluidIndexArray + [0, 0, 0, 1, 1, 1, 2, 2, 2];
storeFluidIndexArray = FluidIndexArray;
temp3 = FluidElementNodes(end,end);
mm = 1;
check = 0;
% all fluid elements for z = 0 --> z = SolidZ done, now --> z = FluidZ
% NOTE: breaks are necessary here!
while (FluidNodes(FluidIndexArray(end)) <= FluidNodes(end))
    if (storeIndex+(numberOfFluidNodesInLowerRow+1)*2*mm > size(FluidNodes,2))
        disp 'break'
        break
    end
    while (FluidNodes(FluidIndexArray(end))<=FluidNodes(storeIndex+(numberOfFluidNodesInLowerRow+1)*2*mm))
        if(FluidNodes(FluidIndexArray(end)) > FluidNodes(end))
            check = 1;
            disp 'break1'
            break
        else
            FluidElementNodes(k,:) = FluidNodes(FluidIndexArray);
            k = k + 1;
            FluidIndexArray = FluidIndexArray + 2;
            if (FluidIndexArray(end) > size(FluidNodes, 2))
                check = 1;
                break
            end
        end
    end
    % from solid to right end
    if (check)
        disp 'break check'
        break
    end
    mm = mm + 1;
    FluidIndexArray = storeFluidIndexArray + (numberOfFluidNodesInLowerRow+1)*2*(mm-1);
end

SortedInterfaceSolidFluidNodes = zeros(2,size(InterfaceGlobalNodesForGeometry,2));
for i = 1:size(SortedInterfaceSolidFluidNodes,2)
    for j = 1:size(SortedInterfaceSolidFluidNodes,2)
        if (InterfaceGlobalNodesForGeometry(j) == i)
            SortedInterfaceSolidFluidNodes(:,i) = SolidInterfaceNodesContinuousToGlobal(:,j);
            break
        end
    end
end
SortedFluidElementNumbersOnInterface = [FluidElementNumbersOnInterface(1:2:end), FluidElementNumbersOnInterface(end-1:-2:2)];

SolidXi1 = zeros(SolidElementNumberZ*2+1,3);
SolidXi2 = SolidXi1;
FluidXi1 = SolidXi1;
FluidXi2 = SolidXi1;
for i = 1:SolidElementNumberZ
    SolidXi1(i, :) = [0.0, 0.0, 0.0];
    SolidXi2(i, :) = [0.0, 0.5, 1.0];
    FluidXi1(i, :) = [1.0, 1.0, 1.0];
    FluidXi2(i, :) = [0.0, 0.5, 1.0];
end
SolidXi1(SolidElementNumberZ+1, :) = [0.0, 0.5, 1.0];
SolidXi2(SolidElementNumberZ+1, :) = [1.0, 1.0, 1.0];
FluidXi1(SolidElementNumberZ+1, :) = [0.0, 0.5, 1.0];
FluidXi2(SolidElementNumberZ+1, :) = [0.0, 0.0, 0.0];
for i = SolidElementNumberZ+2:2*SolidElementNumberZ+1
    SolidXi1(i, :) = [1.0, 1.0, 1.0];
    SolidXi2(i, :) = [1.0, 0.5, 0.0];
    FluidXi1(i, :) = [0.0, 0.0, 0.0];
    FluidXi2(i, :) = [1.0, 0.5, 0.0];
end

SolidInterfaceElements = [1:SolidElementNumberZ, SolidElementNumberZ, ...
    SolidElementNumberZ:-1:1];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% BOUNDARY CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define no displacement BC nodes
NoDisplacement = SolidContinuousToGlobal(1:SolidElementNumberY*2+1);
% define slip BC nodes
SlipFluid = unique((FluidGeometryZ==FluidZ).*FluidNodes);
if (SlipFluid(1)==0)
    SlipFluid = SlipFluid(2:end);
end
%SlipFluid = SlipFluid(2:end);
% define no-slip BC nodes
NoSlipFluid = unique((FluidGeometryZ==0).*FluidNodes);
if (NoSlipFluid(1)==0)
    NoSlipFluid = NoSlipFluid(2:end);
end
% define inlet BC nodes
InletFluid = unique((FluidGeometryY==0).*FluidNodes);
if (InletFluid(1)==0)
    InletFluid = InletFluid(2:end);
end
if (InletFluid(1)==1)
    InletFluid = InletFluid(2:end);
end
% define outlet nodes
ZeroPressure = unique((FluidGeometryY==FluidGeometryY(end)).*FluidNodes);
if (ZeroPressure(1)==0)
    ZeroPressure = ZeroPressure(2:end);
end
OutletFixed = ZeroPressure;
% also remove first node because it is a no-slip BC, last node which is a
% slip BC
% and all nodes that are not part of the linear pressure elements
ZeroPressure = unique((ZeroPressure<=NumberOfFluidPressureNodes).*ZeroPressure);
if (ZeroPressure(1)==0)
    ZeroPressure = ZeroPressure(2:end);
end
%ZeroPressure = ZeroPressure(2:end-1);

% define fixed wall BC (laplace smoothing) nodes
FixedNodes = unique([NoSlipFluid, OutletFixed, InletFluid, SlipFluid]);
if (FixedNodes(1)==0)
    FixedNodes = FixedNodes(2:end);
end
FixedNodesZ = unique((FluidGeometryZ==FluidGeometryZ(end)).*FluidNodes);
if (FixedNodesZ(1)==0)
    FixedNodesZ = FixedNodesZ(2:end);
end
% define moved wall BC (laplace smoothing) nodes
MovedNodes = sort(SortedInterfaceSolidFluidNodes(2,:));
MovedNodesY = FixedNodesZ;

% usage:
% MovedNodesY, FixedNodesZ, MovedNodes, FixedNodes
% ZeroPressure, NoSlipFluid, SlipFluid, InletFluid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% DUPLICATE CHECK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if the geometry contains two nodes with the same location in space
% (time consuming!)
DuplicateCheck = false;
if(DuplicateCheck)
    for i = 1:size(SolidContinuousToGlobal,2)
        for j = 1:size(SolidContinuousToGlobal,2)
            if (i~=j)
                if (SolidGeometryY(i)==SolidGeometryY(j) ...
                        && SolidGeometryZ(i)==SolidGeometryZ(j))
                    print('Duplicate nodes in solid geometry.')
                    print(i,j)
                end
            end
        end
    end
    for i = 1:size(FluidNodes,2)
        for j = 1:size(FluidNodes,2)
            if (i~=j)
                if (FluidGeometryY(i)==FluidGeometryY(j) ...
                        && FluidGeometryZ(i)==FluidGeometryZ(j))
                    print('Duplicate nodes in fluid geometry.')
                    print(i,j)
                end
            end
        end
    end
    for i = 1:size(InterfaceGlobalNodesForGeometry,2)
        for j = 1:size(InterfaceGlobalNodesForGeometry,2)
            if (i~=j)
                if (InterfaceGeometryY(i)==InterfaceGeometryY(j) ...
                        && InterfaceGeometryZ(i)==InterfaceGeometryZ(j))
                    print('Duplicate nodes in interface geometry.')
                    print(i,j)
                end
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

intSpec = '%d ';
realSpec = '%1.16E\n';

solidNodes = fopen('Plate2DSolidNodes.txt','w');
fluidNodes = fopen('Plate2DFluidNodes.txt','w');
interfaceNodes = fopen('Plate2DInterfaceNodes.txt','w');
solidElements = fopen('Plate2DSolidElementNodeNumbers.txt','w');
fluidElements = fopen('Plate2DFluidElementNodeNumbers.txt','w');
interfaceElements = fopen('Plate2DInterfaceElementNodeNumbers.txt','w');
interfaceSolidElements = fopen('Plate2DInterfaceSolidElementNodeNumbers.txt','w');
interfaceFluidElements = fopen('Plate2DInterfaceFluidElementNodeNumbers.txt','w');
interfaceSolidElement = fopen('Plate2DInterfaceSolidElements.txt','w');
interfaceFluidElement = fopen('Plate2DInterfaceFluidElements.txt','w');
solidY = fopen('Plate2DSolidY.txt','w');
solidZ = fopen('Plate2DSolidZ.txt','w');
fluidY = fopen('Plate2DFluidY.txt','w');
fluidZ = fopen('Plate2DFluidZ.txt','w');
interfaceY = fopen('Plate2DInterfaceY.txt','w');
interfaceZ = fopen('Plate2DInterfaceZ.txt','w');

meshConnectivity = fopen('Plate2DSortedInterfaceNodes.txt','w');
solidXi1 = fopen('Plate2DSolidXi1.txt','w');
solidXi2 = fopen('Plate2DSolidXi2.txt','w');
fluidXi1 = fopen('Plate2DFluidXi1.txt','w');
fluidXi2 = fopen('Plate2DFluidXi2.txt','w');

displacement = fopen('Plate2DdisplacementBC.txt','w');
movedY = fopen('Plate2DmovedYNodesBC.txt','w');
fixedZ = fopen('Plate2DfixedZNodesBC.txt','w');
moved = fopen('Plate2DmovedNodesBC.txt','w');
fixed = fopen('Plate2DfixedNodesBC.txt','w');
pressure = fopen('Plate2DpressureBC.txt','w');
noslip = fopen('Plate2DnoSlipBC.txt','w');
inlet = fopen('Plate2DinletBC.txt','w');
slip = fopen('Plate2DslipBC.txt','w');


% % solidNodes = fopen('input/Cantilever2D/SolidNodes.txt','w');
% % fluidNodes = fopen('input/Cantilever2D/FluidNodes.txt','w');
% % interfaceNodes = fopen('input/Cantilever2D/InterfaceNodes.txt','w');
% % solidElements = fopen('input/Cantilever2D/SolidElementNodeNumbers.txt','w');
% % fluidElements = fopen('input/Cantilever2D/FluidElementNodeNumbers.txt','w');
% % interfaceElements = fopen('input/Cantilever2D/InterfaceElementNodeNumbers.txt','w');
% % interfaceSolidElements = fopen('input/Cantilever2D/InterfaceSolidElementNodeNumbers.txt','w');
% % interfaceFluidElements = fopen('input/Cantilever2D/InterfaceFluidElementNodeNumbers.txt','w');
% % interfaceSolidElement = fopen('input/Cantilever2D/InterfaceSolidElements.txt','w');
% % interfaceFluidElement = fopen('input/Cantilever2D/InterfaceFluidElements.txt','w');
% % solidY = fopen('input/Cantilever2D/SolidY.txt','w');
% % solidZ = fopen('input/Cantilever2D/SolidZ.txt','w');
% % fluidY = fopen('input/Cantilever2D/FluidY.txt','w');
% % fluidZ = fopen('input/Cantilever2D/FluidZ.txt','w');
% % interfaceY = fopen('input/Cantilever2D/InterfaceY.txt','w');
% % interfaceZ = fopen('input/Cantilever2D/InterfaceZ.txt','w');
% % 
% % meshConnectivity = fopen('input/Cantilever2D/SortedInterfaceNodes.txt','w');
% % solidXi1 = fopen('input/Cantilever2D/SolidXi1.txt','w');
% % solidXi2 = fopen('input/Cantilever2D/SolidXi2.txt','w');
% % fluidXi1 = fopen('input/Cantilever2D/FluidXi1.txt','w');
% % fluidXi2 = fopen('input/Cantilever2D/FluidXi2.txt','w');
% % 
% % displacement = fopen('input/Cantilever2D/displacementBC.txt','w');
% % movedY = fopen('input/Cantilever2D/movedYNodesBC.txt','w');
% % fixedZ = fopen('input/Cantilever2D/fixedZNodesBC.txt','w');
% % moved = fopen('input/Cantilever2D/movedNodesBC.txt','w');
% % fixed = fopen('input/Cantilever2D/fixedNodesBC.txt','w');
% % pressure = fopen('input/Cantilever2D/pressureBC.txt','w');
% % noslip = fopen('input/Cantilever2D/noSlipBC.txt','w');
% % inlet = fopen('input/Cantilever2D/inletBC.txt','w');
% % slip = fopen('input/Cantilever2D/slipBC.txt','w');




fprintf(solidNodes,'%d\n',size(SolidContinuousToGlobal,2));
fprintf(fluidNodes,'%d\n',size(FluidNodes,2));
fprintf(interfaceNodes,'%d\n',size(InterfaceGlobalNodesForGeometry,2));
fprintf(solidElements,'%d\n',size(SolidElementNodes,1));
fprintf(fluidElements,'%d\n',size(FluidElementNodes,1));
fprintf(interfaceElements,'%d\n',size(InterfaceElementNodes,1));
fprintf(interfaceSolidElements,'%d\n',size(InterfaceElementNodes,1));
fprintf(interfaceFluidElements,'%d\n',size(InterfaceElementNodes,1));
fprintf(interfaceSolidElement,'%d\n',size(SolidInterfaceElements,2));
fprintf(interfaceFluidElement,'%d\n',size(SortedFluidElementNumbersOnInterface,2));
fprintf(solidY,'%d\n',size(SolidGeometryY,2));
fprintf(solidZ,'%d\n',size(SolidGeometryZ,2));
fprintf(fluidY,'%d\n',size(FluidGeometryY,2));
fprintf(fluidZ,'%d\n',size(FluidGeometryZ,2));
fprintf(interfaceY,'%d\n',size(InterfaceGeometryY,2));
fprintf(interfaceZ,'%d\n',size(InterfaceGeometryZ,2));
fprintf(meshConnectivity,'%d\n',size(SortedInterfaceSolidFluidNodes,2));
fprintf(solidXi1,'%d\n',size(InterfaceElementNodes,1));
fprintf(solidXi2,'%d\n',size(InterfaceElementNodes,1));
fprintf(fluidXi1,'%d\n',size(InterfaceElementNodes,1));
fprintf(fluidXi2,'%d\n',size(InterfaceElementNodes,1));
fprintf(displacement,'%d\n',size(NoDisplacement,2));
fprintf(movedY,'%d\n',size(MovedNodesY,2));
fprintf(fixedZ,'%d\n',size(FixedNodesZ,2));
fprintf(moved,'%d\n',size(MovedNodes,2));
fprintf(fixed,'%d\n',size(FixedNodes,2));
fprintf(pressure,'%d\n',size(ZeroPressure,2));
fprintf(noslip,'%d\n',size(NoSlipFluid,2));
fprintf(inlet,'%d\n',size(InletFluid,2));
fprintf(slip,'%d\n',size(SlipFluid,2));


fprintf(interfaceSolidElement,intSpec,SolidInterfaceElements);
fprintf(interfaceFluidElement,intSpec,SortedFluidElementNumbersOnInterface);
fprintf(solidNodes,'%d\n',SolidContinuousToGlobal);
fprintf(fluidNodes,'%d\n',FluidNodes);
fprintf(interfaceNodes,'%d\n',InterfaceGlobalNodesForGeometry);
for i = 1:size(SolidElementNodes,1)
    fprintf(solidElements,intSpec,SolidElementNodes(i,:));
    fprintf(solidElements,'\n');
end
for i = 1:size(FluidElementNodes,1)
    fprintf(fluidElements,intSpec,FluidElementNodes(i,:));
    fprintf(fluidElements,'\n');
end
for i = 1:size(InterfaceElementNodes,1)
    fprintf(interfaceElements,intSpec,InterfaceElementNodes(i,:));
    fprintf(interfaceElements,'\n');
    fprintf(interfaceSolidElements,intSpec,SortedInterfaceSolidFluidNodes(1,1+2*(i-1):3+2*(i-1)));
    fprintf(interfaceSolidElements,'\n');
    fprintf(interfaceFluidElements,intSpec,SortedInterfaceSolidFluidNodes(2,1+2*(i-1):3+2*(i-1)));
    fprintf(interfaceFluidElements,'\n');
    fprintf(solidXi1,'%1.16E ',SolidXi1(i,:));
    fprintf(solidXi1,'\n');
    fprintf(solidXi2,'%1.16E ',SolidXi2(i,:));
    fprintf(solidXi2,'\n');
    fprintf(fluidXi1,'%1.16E ',FluidXi1(i,:));
    fprintf(fluidXi1,'\n');
    fprintf(fluidXi2,'%1.16E ',FluidXi2(i,:));
    fprintf(fluidXi2,'\n');
end
fprintf(solidY,realSpec,SolidGeometryY);
fprintf(solidZ,realSpec,SolidGeometryZ);
fprintf(fluidY,realSpec,FluidGeometryY);
fprintf(fluidZ,realSpec,FluidGeometryZ);
fprintf(interfaceY,realSpec,InterfaceGeometryY);
fprintf(interfaceZ,realSpec,InterfaceGeometryZ);
for i = 1:size(SortedInterfaceSolidFluidNodes,2)
    fprintf(meshConnectivity,'%d ',SortedInterfaceSolidFluidNodes(:,i));
    fprintf(meshConnectivity,'\n');
end
fprintf(displacement,intSpec,NoDisplacement);
fprintf(movedY,intSpec,MovedNodesY);
fprintf(fixedZ,intSpec,FixedNodesZ);
fprintf(moved,intSpec,MovedNodes);
fprintf(fixed,intSpec,FixedNodes);
fprintf(pressure,intSpec,ZeroPressure);
fprintf(noslip,intSpec,NoSlipFluid);
fprintf(inlet,intSpec,InletFluid);
fprintf(slip,intSpec,SlipFluid);


fclose(solidNodes);
fclose(fluidNodes);
fclose(interfaceNodes);
fclose(solidElements);
fclose(fluidElements);
fclose(interfaceElements);
fclose(interfaceSolidElements);
fclose(interfaceFluidElements);
fclose(interfaceSolidElement);
fclose(interfaceFluidElement);
fclose(solidY);
fclose(solidZ);
fclose(fluidY);
fclose(fluidZ);
fclose(interfaceY);
fclose(interfaceZ);
fclose(meshConnectivity);
fclose(solidXi1);
fclose(solidXi2);
fclose(fluidXi1);
fclose(fluidXi2);
fclose(displacement);
fclose(movedY);
fclose(fixedZ);
fclose(moved);
fclose(fixed);
fclose(pressure);
fclose(noslip);
fclose(inlet);
fclose(slip);


if (showMeshes == 1)
    plot(SolidGeometryPressureY, SolidGeometryPressureZ, 'r*', FluidGeometryPressureY, FluidGeometryPressureZ, 'g+')
elseif (showMeshes == 2)
    plot(SolidGeometryY, SolidGeometryZ, 'r*', FluidGeometryY, FluidGeometryZ, 'g+')
end
