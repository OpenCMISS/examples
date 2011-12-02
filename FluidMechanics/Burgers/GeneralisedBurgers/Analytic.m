global a;
global b;
global c;
global D;
global E;
global L;
global ne ;
global np;
global theta;
global deltat;
global startt;
global stopt;
global t;

a = 1.0;
b = -1.0;
c = 1.0;
D = 1.0;
E = 1.0;
L = 3.0;
ne = 6;
np = ne+1;
theta = 0.5;
deltat = 0.005;
startt = 0.0;
stopt = 1.0;

t = startt;

global x = zeros(np,1);

global K = zeros(np,np);
global C = zeros(np,np);

K1 = b*ne/L;
K2 = -2*b*ne/L;
C1 = a*L/(6*ne);
C2 = (2*a*L)/(3*ne);

x(1) = 0.0;
alpha(1) = 0.0;

K(1,1) = K2/2;
K(1,2) = K1;

C(1,1) = C2/2;
C(1,2) = C1;

for i = 2:np-1

  x(i) = (i-1)*L/ne;

  K(i,i-1) = K1;
  K(i,i) = K2;
  K(i,i+1) = K1;

  C(i,i-1) = C1;
  C(i,i) = C2;
  C(i,i+1) = C1;
    
endfor

x(np) = L;

K(np,np-1) = K1;
K(np,np) = K2/2;

C(np,np-1) = C1;
C(np,np) = C2/2;

#function [ analu, analv ] = analytic( time )
#
#  global a;
#  global b;
#  global c;
#  global D;
#  global E;
#  global x;
#
#  analu = (D+a*x)/(E+c*time);
#  analv = -c*(D+a*x)/((E+c*time)*(E+c*time));
#
#endfunction

function [ analu, analv ] = analytic( t )

  global a;
  global b;
  global c;
  global D;
  global E;
  global x;
  global np;

  analu = zeros(np,1);
  analv = zeros(np,1);

  for i=1:np
    analu(i) = a*D+2*b/(c*(x(i)-c*D*t+E));
    analv(i) = 2*b*c*D/(c*(x(i)-c*D*t+E)*(x(i)-c*D*t+E));
  endfor

endfunction

function [ err, normerr ] = error( z, analz )

  [ numrow, numcol ] = size( z );
  err = zeros(numrow,1);
  sum = 0.0;
  for i = 1:numrow
    if( abs(analz(i)) > 0.00000001 )
      err(i) = 100.0*(z(i)-analz(i))/analz(i);
    else
      err(i) = 0.0;
    endif
    sum = sum + err(i)*err(i);
  endfor
  normerr = sqrt(sum);

endfunction

function y = g ( z )

   global c;
   global np;

   y = zeros(np,1);
  
   y(1)=c/6.0*(z(2)*z(2)+z(1)*z(2)-2*z(1)*z(1));

   for i = 2:np-1

     y(i)= c/6.0*(2*z(i+1)*z(i+1)+z(i-1)*z(i)-z(i)*z(i+1)-2*z(i-1)*z(i-1));

   endfor
   
   y(np) = c/6.0*(2*z(np)*z(np)-z(np-1)*z(np)-z(np-1)*z(np-1));

endfunction 

function J = delgdelu ( z )

   global c;
   global np;

   J = zeros(np,np);
   
   J(1,1) = c/6.0*(z(2)-4*z(1));
   J(1,2) = c/6.0*(2*z(2)+z(1));

   for i = 2:np-1

     J(i,i-1) = c/6.0*(-z(i)-2*z(i-1));
     J(i,i) = c/6.0*(z(i+1)-z(i-1));
     J(i,i+1) = c/6.0*(2*z(i+1)+z(i));

   endfor

   J(np,np-1) = c/6.0*(-z(np)-2*z(np-1));
   J(np,np) = c/6.0*(4*z(np)-z(np-1));

endfunction

function y = newg ( myalpha )
  
   global predictedu;
   global deltat;

   newu = predictedu + deltat*myalpha;

   y = g( newu );

endfunction 

function r = residual ( myalpha )

   global A;
   global prevg;
   global beta;
   global theta;

   r = A*myalpha + theta*newg( myalpha ) + (1-theta)*prevg + beta;
	 
endfunction

function jacnewg = delnewgdelu ( alpha )

   global predictedu;
   global deltat;

   newu = predictedu + deltat*alpha;

   jacnewg = delgdelu ( newu );

endfunction

function J = Jacobian ( myalpha )

   global A;
   global theta;
   global deltat;

   J = A + theta*deltat*delnewgdelu ( myalpha );

endfunction

function [ reducedr, reducedJacobian ] = reducedfunction ( reducedalpha )

   global alpha;
   global np;

   myalpha = zeros(np,1); 
  
   myalpha(1) = alpha(1);
   myalpha(2:np-1) = reducedalpha(1:np-2);
   myalpha(np) = alpha(np);

   r = residual( myalpha );

   reducedr = r(2:np-1);

   if(nargout > 1)

      J = Jacobian( myalpha );

      reducedJ = J(2:np-1,2:np-1);

   endif

endfunction

global A = zeros(np,np);
global u = zeros(np,1);
global v = zeros(np,1);
global analyticu = zeros(np,1);
global analyticv = zeros(np,1);
global meanpredictedu = zeros(np,1);
global predictedu = zeros(np,1);
global prevu = zeros(np,1);
global currentg = zeros(np,1);
global analyticg = zeros(np,1);
global prevg = zeros(np,1);
global alpha = zeros(np,1);
global reducedalpha = zeros(np-2,1);
global beta = zeros(np,1);
global resid = zeros(np,1);
global psi = zeros(np,1);
global analyticresid = zeros(np,1);
global uerr = zeros(np,1);
global uerrnorm = 0.0;

A = C + theta*deltat*K;

#cond(A)

[ u, analyticv ] = analytic( t );

currentg = g( u );

while t < stopt

  t = t + deltat

  [ analyticu, analyticv ] = analytic( t )

  prevu = u

  prevg = currentg

  meanpredictedu = u

  predictedu = u

  alpha = ( analyticu - predictedu )/deltat

  beta = K*meanpredictedu;

  startalpha = zeros(np-2,1);
  
  #startalpha = analyticv(2:np-1);

  [ reducedalpha, fval, info ] = fsolve( @reducedfunction, startalpha );

  alpha(2:np-1) = reducedalpha

  u = predictedu + deltat*alpha

  v = (u - prevu)/deltat;

  currentg = g( u );

  analyticg = g( analyticu );

  #A*alpha

  #beta

  psi = A*alpha + theta*currentg + (1-theta)*prevg + beta;

  resid = C*v + K*u + currentg;

  analyticresid =  C*analyticv + K*analyticu + analyticg;

  [ uerr, unormerr ] = error( u, analyticu )

  unormerr
  
  plot(x, u, x, analyticu,x,resid)
  
  pause(0.05)

  #quit

endwhile
