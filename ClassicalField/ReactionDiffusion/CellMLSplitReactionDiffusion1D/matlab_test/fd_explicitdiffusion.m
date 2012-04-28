%Author: Vijay Rajagopal

%fd_diffusion.m: Finite difference solution to the 1D diffusion equation u_t-d/dx(Du_x)=0, where D is spatially varying like in cells with different diffusivity 
%coefficients in different spaces

%Define a 2D grid of points of length and time xdim and tdim and set up empty matrix to fill concentration field
close all;
clear all;
xdim=1;
tdim=0.25;

dx=0.1;
dt=0.01;

x=[0:dx:xdim];
t=[0:dt:tdim];
u=zeros(length(t),length(x));
%medium properties: diffusion coefficeint, D as spatially varying.

% D=1
%D=ones(1,length(x));

% D=ax
D=0.4*x;
%set up initial conditions u(0,x)=sin(pi*x)
for xi=1:length(x)
  u(1,xi)=sin(pi*x(xi));
end
%set up BC conditions u(t,0)=0, u(t,1)=0
u(1:length(t),1)=0;
u(1:length(t),length(x))=0;

%fd approx of u_t=[u(t+dt,x)-u(t,x)]/dt
%central diff fd approx of u_xx=[u(t,x+dx)-2u(t,x)+u(t,x-dx)]/dx^2
%d/dx(Du_x) = D_x*u_x+Du_xx
%therefore diffusion equation with spatial variation of diffusion coeff:
%u_t=D_x*u_x+D*u_xx
%therefore value of u at next time step aftersome finite difference approximation workings in my notebook:
% u(t+dx,x)=u(t,x)+(dt/dx^2)*[(D(t,x+dx/2)-D(t,x-dx/2))(u(t,x+dx/2)-u(t,x-dx/2))+D*(u(t,x+dx)-2u(t,x)+u(t,x-dx))]
 
stepcoeff=(dt/(dx^2));
%solve for the remainder of time after t=0
for ti=1:(length(t)-1)
  for xi=2:(length(x)-1)
    utxplushalf=(u(ti,xi)+u(ti,xi+1))/2;
    utxminushalf=(u(ti,xi)+u(ti,xi-1))/2;
    Dtxplushalf=(D(xi)+D(xi+1))/2;
    Dtxminushalf=(D(xi)+D(xi-1))/2;
    u(ti+1,xi)=u(ti,xi)+(stepcoeff*D(xi)*(u(ti,xi+1)-2*u(ti,xi)+u(ti,xi-1)))+(stepcoeff*(Dtxplushalf-Dtxminushalf)*(utxplushalf-utxminushalf));
  end
end
surf(x,t,u);




