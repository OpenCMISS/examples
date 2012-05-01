%Author: Vijay Rajagopal

%fe_diffusion_strang_constantsource.m: Finite element solution to the 1D diffusion equation u_t-d/dx(Du_x)=0, where D is a constant.
%solving using first order splitting techniques to test code against opencmiss code.

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
dt = 0.001;
odestep = 0.00001;
t = [0:dt:T];

u = zeros(numNod,size(t,2));
%set up initial condition
for node = 1:numNod
	u(node,1) = 0; %sin(pi*nodes(node,2));
end
%set up BC

u(1,1:size(t,2)) = 1.5;
u(numNod,1:size(t,2)) = 1.5;

%source with linear function in x - source = 0.01x
%source_vec = [0:0.1:1]'; %for x in range 0 to 100
%source_gradinx = 0.01;
%loop through each time step

for tstep = 1:(size(t,2)-1)
%_______________________________________________________________________________________________________________________________
%ODE SOLVE 1
        %First solve the ode part for tstep<timestep<tstep+dt
        initT = t(tstep);
        halfT = t(tstep)+(dt*0.5);
        finalT = t(tstep)+dt;

	time_interval = finalT-initT;
	%evaluate the ode at each of the spatially discretized nodes.
	tintvl = [initT:odestep:finalT];
	uintmd_one = zeros(numNod,size(tintvl,2));
	
	for node = 1:numNod
		uintmd_one(node,:) = lsode("constfunc",u(node,tstep),tintvl);
	end
%________________________________________________________________________________________________________________________________
%PDE SOLVE
	%we now need to solve for the set of uj's over the interval tstep to tstep+dt (initT to finalT) now
	%on the 2nd order diffusion part; i.e. delU/delT = del2U/delXsquared
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
		%set the initial values of u for the pde solver to be the u values at dt/2 from the previous ode solve.
		ur = [uintmd_one(nodes(elems(elem,1)),size(uintmd_one,2));uintmd_one(nodes(elems(elem,2)),size(uintmd_one,2))];
		%note source term not evaluated anymore because it is in the ode component of the strang split
		fe(1,1) = fe(1,1)+(1/dt)*(Ae(1,:)*ur);
		fe(2,1) = fe(2,1)+(1/dt)*(Ae(2,:)*ur);
		f(globNod1,1) = f(globNod1,1)+fe(1,1);
		f(globNod2,1) = f(globNod2,1)+fe(2,1);
	end
	C = ((1/dt)*A)+B;

	%solve pde for current full time step
        %since Dirchlet's have been prescribed at BCs, first and last row of
        %the C matrix must be removed or replaced with the appropriate
        %coefficients
        f(1,1) = u(1,2);
        f(numNod,1) = u(numNod,2);
        C(1,1) = 1;
        C(1,2:numNod) = 0;
        C(numNod,numNod) = 1;
        C(numNod,1:numNod-1) = 0;
	uintmd_two = inv(C)*(f);
	%the final u value for the overall pde is the value from this last ode solver at dt/2
	u(:,tstep+1) = uintmd_two(1:numNod,size(uintmd_two,2));
end		
