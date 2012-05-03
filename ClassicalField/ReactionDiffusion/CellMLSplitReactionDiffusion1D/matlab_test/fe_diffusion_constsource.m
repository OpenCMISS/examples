%Author: Vijay Rajagopal

%fe_diffusion.m: Finite element solution to the 1D diffusion equation u_t-d/dx(Du_x)=0, where D is a constant.

close all
clear all

%set up the 1d problem defined in the finite difference code: x going from 0 to 1; t goes from 0 to 0.25
%initial condition, u(0,x) = sin(pi*x), boundary condition, u(t,0) = u(t,1) = 0

%set up diffusivity params
D = 0.5;

%set up geometry, element parameters
xL = 100.0;
numEl = 10;
numNod = numEl+1;
[nodes elems] = createFEmesh(numEl,numNod,xL);
xi = [0.2113248654 0.7886751346];
w = [0.5 0.5];


%set up time parameters
T = 0.49;
dt = 0.0001;
t = [0:dt:T];

u = zeros(numNod,size(t,2));
%set up initial condition
for node = 1:numNod
	u(node,1) = 0; %sin(pi*nodes(node,2));
end
%set up BC

u(1,1:size(t,2)) = 1.5;
u(numNod,1:size(t,2)) = 1.5;

%set up source - constant through space and time
source = 0.5;
source_vec = ones(numNod,1)*source;

%source with linear function in x - source = 0.01x
%source_vec = [0:0.1:1]'; %for x in range 0 to 100
%source_gradinx = 0.01;
%loop through each time step
for tstep = 1:(size(t,2)-1)
	%we need to solve for the set of uj's at this tstep
	%using the equation:
	%Adudx+Bu = f (see note book)
	%Using the backward difference approximation to the time derivative,
	%Lapidus & Pinder set up the system of ODEs for solving the diffusion eqn
	A = zeros(numNod,numNod);
	B = zeros(numNod,numNod);
	C = zeros(numNod,numNod);
	f = zeros(numNod,1);

	%loop through each element
	for elem = 1:numEl
		Ae = zeros(2,2);
		%calculate element A coeffs
		Je = nodes(elems(elem,2),2)-nodes(elems(elem,1),2);
		Ae(1,1) = Je*(((1-xi(1))*(1-xi(1))*w(1))+((1-xi(2))*(1-xi(2))*w(2)));
		Ae(1,2) = Je*(((1-xi(1))*(xi(1))*w(1))+((1-xi(2))*(xi(2))*w(2)));
		Ae(2,1) = Je*(((1-xi(1))*(xi(1))*w(1))+((1-xi(2))*(xi(2))*w(2)));
		Ae(2,2) = Je*(w(1)*(xi(1)^2)+w(2)*(xi(2)^2));

		%assemble these into global A
		globNod1 = nodes(elems(elem,1));
		globNod2 = nodes(elems(elem,2));
		A(globNod1,globNod1) = A(globNod1,globNod1)+Ae(1,1);
		A(globNod1,globNod2) = A(globNod1,globNod2)+Ae(1,2);
		A(globNod2,globNod1) = A(globNod2,globNod1)+Ae(2,1);
		A(globNod2,globNod2) = A(globNod2,globNod2)+Ae(2,2);
		
		Be = zeros(2,2);
		Be(1,1) = D/Je;
		Be(1,2) = -1*D/Je;
		Be(2,1) = -1*D/Je;
		Be(2,2) = D/Je;
		B(globNod1,globNod1) = B(globNod1,globNod1)+Be(1,1);
		B(globNod1,globNod2) = B(globNod1,globNod2)+Be(1,2);
		B(globNod2,globNod1) = B(globNod2,globNod1)+Be(2,1);
		B(globNod2,globNod2) = B(globNod2,globNod2)+Be(2,2);

		fe = zeros(2,1);
		ur = [u(nodes(elems(elem,1)),tstep);u(nodes(elems(elem,2)),tstep)];
                source_e = zeros(2,1);
                %constant source in space, time, and independent of current conc.
	        source_e(1,1) = source_vec(globNod1,1)*Je*(((1-xi(1))*w(1)) + ((1-xi(2))*w(2)));
	        source_e(2,1) = source_vec(globNod2,1)*Je*(((xi(1))*w(1)) + ((xi(2))*w(2)));
		%source described as a linear function in space; source = 0.01*x (x in range 0 to 100)
                % 
	        %source_e(1,1) = source_vec(globNod1,1)*Je*(((1-xi(1))*w(1)) + ((1-xi(2))*w(2)))
	        %source_e(2,1) = source_vec(globNod2,1)*Je*(((xi(1))*w(1)) + ((xi(2))*w(2)))
		fe(1,1) = fe(1,1)+source_e(1,1)+(1/dt)*(Ae(1,:)*ur);
		fe(2,1) = fe(2,1)+source_e(2,1)+(1/dt)*(Ae(2,:)*ur);
		f(globNod1,1) = f(globNod1,1)+fe(1,1);
		f(globNod2,1) = f(globNod2,1)+fe(2,1);
	end
	C = ((1/dt)*A)+B;

	%solve for current time step
    %since Dirchlet's have been prescribed at BCs, first and last row of
    %the C matrix must be removed or replaced with the appropriate
    %coefficients
    f(1,1) = u(1,2);
    f(numNod,1) = u(numNod,2);
    C(1,1) = 1;
    C(1,2:numNod) = 0;
    C(numNod,numNod) = 1;
    C(numNod,1:numNod-1) = 0;
    u(:,tstep+1) = inv(C)*(f);
    
end		
