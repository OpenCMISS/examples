%Author: Vijay Rajagopal

%fd_diffusion.m: Finite difference solution to the 1D diffusion equation u_t-d/dx(Du_x)=0, where D is a constant.
%In this code we use implicit finite difference method for solving the problem
% Setting up spatially varying D is pretty easy. I have calculated the formula
% for implicit method with spatial varying D in my notebook.

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
%medium properties: diffusion coefficeint, D is constant

 D=1;
%D=ones(1,length(x));

% D=ax
%D=0.4*x;
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
%therefore value of u at next time step after some finite difference approximation workings in my notebook:
% u(t,x)=(1+2r)*u(t+dt,x)-ru(t+dt,x+dx)-ru(t+dt,x-dx)
%where r=dt/dx^2
%u(t+dt,xi's) are all unknowns. Sets up system of linear equations to be solved.
 
r=(dt/(dx^2));
b=u(1,2:length(x)-1);
b=b';
%set of matrix of eqn. 2 rows and columns missing because of bc's having been set up.
A=zeros(length(x)-2,length(x)-2);
%solve for the remainder of time after t=0
for ti=2:(length(t))
  for xi=1:(length(x)-2)
   %if statements to set up equations only for unknowns. 2 of the u variables at the ends of the 1D line are known cos they are set as BCs. 
   %this means we can delete the row and the column in A that correspond to each of these known variables. 
    if xi~=1
      A(xi,xi-1)=-1*r;
    else
      b(xi)=b(xi)+r*u(ti,1);
    end
    A(xi,xi)=(1+2*r);
    if xi~=(length(x)-2)
      A(xi,xi+1)=-1*r;
    else
      b(xi)=b(xi)+r*u(ti,length(x));
    end
  end
%solve for X - u(t+dt,xs) using system of equations AX=b
  X=inv(A)*b;
  X=X';
  u(ti,2:(length(x)-1))=X;
  b=X';
end
surf(x,t,u);




