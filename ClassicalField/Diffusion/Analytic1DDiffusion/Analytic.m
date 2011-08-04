global a;
global b;
global c;
global k;
global L;
global ne ;
global np;
global theta;
global deltat;
global startt;
global stopt;
global t;

a = 1.0;
b = pi/2;
c = 0;
k = -1;
L = 3;
ne = 20;
np = ne+1;
theta = 0.5;
deltat = 0.01;
startt = 0.0;
stopt = 0.5;

t = startt;

global x = zeros(np,1);
global alpha = zeros(np,1);

global K = zeros(np,np);
global C = zeros(np,np);

K1 = -k*ne/L;
K2 = 2*k*ne/L;
C1 = a*L/(6*ne);
C2 = (2*a*L)/(3*ne);

x(1) = 0.0;

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

function [ analu, analv ] = analytic( time )

  global a;
  global b;
  global c;
  global k;
  global L;
  global x;

  analu = a*exp(4*pi*pi*k*time/(L*L))*cos(2*pi*x/L+b)+c;
  analv = 4*pi*pi*k*a/(L*L)*exp(4*pi*pi*k*time/(L*L))*cos(2*pi*x/L+b);

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
  normerr(1) = sqrt(sum);

endfunction

global A = zeros(np,np);
global prevu = zeros(np,1);
global meanpredictedu = zeros(np,1);
global predictedu = zeros(np,1);
global beta = zeros(np,1);
global r1 = zeros(np,1);
global r2 = zeros(np,1);
global rhs = zeros(np,1);
global reducedA = zeros(np-2,np-2);
global reducedrhs = zeros(np-2,1);
global reducedalpha = zeros(np-2);
global u = zeros(np,1);
global v = zeros(np,1);
global analyticu = zeros(np,1);
global analyticv = zeros(np,1);
global resid = zeros(np,1);
global analyticresid = zeros(np,1);
global uerr = zeros(np,1);
global uerrnorm = 0.0;

A = C + theta*deltat*K;

unit = ones(np,1)

A*unit

[ u, analyticv ] = analytic( t );

while t < stopt

  t = t + deltat

  prevu = u;

  [ analyticu, analyticv ] = analytic( t );

  meanpredictedu = u;

  predictedu = u;

  alpha = ( analyticu - predictedu )/deltat;

  beta = K*meanpredictedu;

  r1 = A(:,1)*alpha(1);

  r2 = A(:,np)*alpha(np);

  rhs = -beta - r1 - r2;

  reducedA = A(2:np-1,2:np-1);

  reducedrhs = rhs(2:np-1)

  reducedalpha = reducedA\reducedrhs

  alpha(2:np-1) = reducedalpha(1:np-2);

  u = predictedu + deltat*alpha;

  v = (u - prevu)/deltat;

  resid = C*v + K*u;

  analyticresid = C*analyticv + K*analyticu;

  [ uerr, unormerr ] = error( u, analyticu );

  unormerr
  
  max(u)
  max(analyticu)

  plot(x, u, x, analyticu)
  
  pause(1)
  
  #quit

endwhile
