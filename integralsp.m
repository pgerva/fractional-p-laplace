function [Isp]=integralsp(u,s,p,x)
% INTEGRALSP: evaluate P.V. \int_R  g_x(y) dy 
% with g_x(y) = abs(u(x)-u(y))^(p-2)*(u(x)-u(y))/abs(x-y)^(1+sp)
%        s in (0,1), p>1, x in a interval
% [Isp]=integralsp(u,s,p,x)
%
% Input: u = @(x).... function handle of the function appearing in g_x
%        s =  fractional order of the differential operator
%        p =  index of the p-Laplace operator
%        x = point in (-1,1) at which evaluate the integral
%
% Output: Isp = P.V. \int_R  g_x(y) dy 
%        with g_x(y) = abs(u(x)-u(y))^(p-2)*(u(x)-u(y))/abs(x-y)^(1+sp)
%        It holds: (-\Delta_p)^s u(x) =  C_{s,p} Isp
%        with
%        C_{s,p}=s*p/2*(1-s)*2^(2*s-1)*Gamma((1+s*p)/2)/
%                                           (Gamma((p+1)/2)*Gamma(2-s))
% 
%   References: 
%   [1] F. Colasuonno, F. Ferrari, P. Gervasio, and A. Quarteroni.
%       "Some evaluations of the fractional p-Laplace operator on radial
%       functions" (2021). 
%       Preprint: https://arxiv.org/abs/2112.08239
%
%   [2] C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%       "Spectral Methods. Fundamentals in Single Domains"
%       Springer Verlag, Berlin Heidelberg New York, 2006.
%
%   Copyright (C) 2007 - 2022  by Paola Gervasio
% 
%   This work is licensed under a Creative Commons 
%   Attribution-NonCommercial-ShareAlike 4.0 International License
%   http://creativecommons.org/licenses/by-nc-sa/4.0/
%
%   This file is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  

delta=1/50; 
np=256; np1=np+1;

gx_16=@(y)(abs(u(x))).^(p-2).*u(x)./((abs(x-y)).^(1+p*s));

gx=@(y)(abs(u(x)-u(y))).^(p-2).*(u(x)-u(y))./((abs(x-y)).^(1+p*s));

% integral: MATLAB function
I(1)=integral(gx_16,-Inf,-1,'RelTol',1e-12,'AbsTol',1e-15);
I(2)=gl_quadrature(gx,-1,x-delta,np1);
I(3)=gl_quadrature(gx,x-delta,x,np);
I(4)=gl_quadrature(gx,x,x+delta,np);
I(5)=gl_quadrature(gx, x+delta,1,np1);
I(6)=integral(gx_16,1,+Inf,'RelTol',1e-12,'AbsTol',1e-15);

Isp= sum(I);

return

function [gl_integral]=gl_quadrature(f,a,b,np)
% GL_QUADRATURE: Gauss Legendre quadrature formula to integrate f on (a,b)
% [gl_integral]=gl_quadrature(f,a,b,np)
%
% Input: f = @(y)... function handle of the function to integrate
%        a, b = end points of the integration interval (a < b)
%        np = number of points for the quadrature formula
%
% Output: gl_integral = value of the computed integral
%
%
% Reference: F. Colasuonno, F. Ferrari, P. Gervasio, and A. Quarteroni.
% Some evaluations of the fractional p-Laplace operator on radial
% functions. (2022) Mathematics in Engineering
%
%   Copyright (C) 2007 - 2022  by Paola Gervasio


if a > b
  warning("The end-points are not sorted correctly")
  gl_integral= [] ;
  return
end

% compute nodes and weights of the GL quadrature formula
[x,w]=xwlg(np,a,b);
x=x(:);
gl_integral=f(x)'*w;

return

function [x,w] = xwlg(np,a,b)
%XWLG  Computes nodes and weights of the Legendre-Gauss  quadrature formula.
%
%    [x,w]=xwlg(np) returns the np weigths and nodes 
%    of the corresponding Legendre Gauss quadrature 
%    formula in the reference interval (-1,1).
%
%    [x,w]=xwlg(np,a,b) returns the np weigths and the nodes 
%    of the corresponding Legendre Gauss quadrature 
%    formula in the  interval (a,b).
%
% Input: np = number of nodes
%        a, b = extrema of the interval
%
% Output: x(np,1) = LG nodes  (CHQZ2, (2.3.10), pag. 76)
%         w(np,1) = LG weigths (CHQZ2, (2.3.10), pag. 76)
%
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.
%
%   Copyright (C) 2007  by Paola Gervasio


if np<=1
  x=0;w=2;
  return
end
x=jacobi_roots(np,0,0);
w=2./(plegendre_der(x,np).^2.*(1-x.^2));

%
% map on (a,b)
%
if nargin == 3
  bma=(b-a)*.5;
  bpa=(b+a)*.5;
  x=bma*x+bpa;
  w=w*bma;
end
return

function [x,flag] = jacobi_roots(n,alpha,beta) 
% JACOBI_ROOTS: Computes the n zeros of the Jacoby polynomial P_n^{(\alpha,\beta)}(x)
%   by Newton method and deflation process.
%
%    [x] = jacobi_roots(n,alpha,beta) 
%
% Input: n = polynomial degree
%        alpha,beta = parameters of Jacobi polynomial
%
% Output: x = zeros of Jacobi polynomial (column array of size n)
%         flag  = -1: The polynomial degree should be greater than 0
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.
%
%   Copyright (C) 2007  by Paola Gervasio


flag=0;
if n < 1 
  es='The polynomial degree should be greater than 0';
  disp(es); flag = -1; x = []; 
return 
 
else
x=zeros(n,1);

x0=cos(pi/(2*n));
tol=1e-14; kmax=15;
for j=1:n
 diff=tol+1;kiter=0;
    while kiter <=kmax && diff>=tol
        [p,pd]=jacobi_eval(x0,n,alpha,beta);
% deflation process q(x)=p(x)/((x-x_1)*... (x-x_{j-1}))
% q(x)/q'(x)=p(x)/[p'(x)-p(x)*\sum_{i<j} 1/(x-x_i)]
        ss=sum(1./(x0-x(1:j-1)));
        x1=x0-p/(pd-ss*p);
        diff=abs(x1-x0);
        kiter=kiter+1;
        x0=x1;
    end
    x0=0.5*(x1+cos((2*(j+1)-1)*pi/(2*n)));
    x(j)=x1;
end 
x=sort(x);
end


return

function [polder,pol] = plegendre_der (x, n) 
% PLEGENDRE_DER    Evaluates the first derivative of Legendre polynomial
%                  of degree n 
%
%    [polder,pol]=plegendre_der(x,n)  evaluates (L_n)'(x), L_n(x)
%    at the node(s) x, using the three term relation  (2.3.3), pag. 75 CHQZ2.
%
% Input: x = scalar or column array
%        n = degree of Legendre polynomial
%
% Output: polder = scalar or column array (same dimension as x)
%             with the evaluation of (L_n)'(x)
%         pol  = scalar or column array (same dimension as x)
%             with the evaluation of L_n(x)
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.
%
%   Copyright (C) 2007  by Paola Gervasio


nn=size(x);

if nn==1
  % x is a scalar
  if n==0
    pol=1; polder=0;
  elseif n==1
    pol=x; polder=1;
  else
    polkp1=0; polkp1der=0;
    polkm1=1;  polk=x;
    polkm1der=0; polkder=1;
    for k =1:n-1
      duekp1=2*k+1; kp1=1/(k+1);
      polkp1=(duekp1*x*polk-k*polkm1)*kp1;
      polkp1der=(duekp1*(x*polkder+polk)-k*polkm1der)*kp1;
      polkm1der=polkder; polkder=polkp1der;
      polkm1=polk; polk=polkp1;
    end
    polder=polkp1der; pol=polkp1;
  end

else
x=x(:);

if n==0
pol=ones(nn);polder=zeros(nn);
elseif n==1
pol=x; polder=ones(nn);
else

polk=x; polkder=ones(nn);
polkm1=ones(nn); polkm1der=zeros(nn);
polkp1=zeros(nn); polkp1der=zeros(nn);
for k =1:n-1
   duekp1=2*k+1;  kp1=1/(k+1);
   polkp1=(duekp1*x.*polk-k*polkm1)*kp1; 
   polkp1der=(duekp1*(x.*polkder+polk)-k*polkm1der)*kp1; 
   polkm1der=polkder; polkder=polkp1der; 
   polkm1=polk; polk=polkp1; 
end 
polder=polkp1der; pol=polkp1;
end
end
return 


function [p,pd] = jacobi_eval(x,n,alpha,beta) 
% JACOBI_EVAL: Evaluates Jacobi polynomial P_n^{(\alpha,\beta)} at x\in(-1,1)
%
%              (formula (2.5.3) pag. 92, CHQZ2)
%  [p,pd] = jacobi_eval(x,n,alpha,beta) 
%
% Input: x = scalar or one-dimensional array of length (m) 
%        n =  polynomial degree
%        alpha, beta= parameters of Jacoby polynomial
%
% Output: p(m,3) = [P_n^{(\alpha,\beta)}(x),
%                   P_(n-1)^{(\alpha,\beta)}(x),
%                   P_(n-2)^{(\alpha,\beta)}(x)];
%
%         pd(m,3) = [(P_n^{(\alpha,\beta)})'(x),
%                    (P_(n-1)^{(\alpha,\beta)})'(x),
%                    (P_(n-2)^{(\alpha,\beta)})'(x)];
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.
%
%   Copyright (C) 2007  by Paola Gervasio

apb=alpha+beta; ab2=alpha^2-beta^2; 

nn=length(x);
if nn==1
% x is a scalar
p=1;   pd=0; 
if n == 0    
   return 
elseif n==1
p = (alpha-beta+(apb+2)*x)*0.5; 
 
pd = 0.5*(apb+2); 
else
p1 = p;
p = (alpha-beta+(apb+2)*x)*0.5; 

pd1 = pd; 
pd = 0.5*(apb+2); 
for k = 1:n-1 
   k1=k+1; k2=k*2; k2ab=k2+alpha+beta;
   k2ab1=k2ab+1; k2ab2=k2ab1+1;
   p2=p1; p1=p; 
   pd2=pd1; pd1=pd; 
   a1 = 2*k1*(k1+apb)*k2ab; 
% Gamma(x+n)/Gamma(x) = (x+n-1) * ... * (x+1) * x 
   a21 = k2ab1*ab2;
   a22 = k2ab2*k2ab1*k2ab; 
   a3=2*(k+alpha)*(k+beta)*k2ab2;
   p = ((a21+a22*x)*p1-a3*p2)/a1; 
   pd= (a22*p1+(a21+a22*x)*pd1-a3*pd2)/a1; 
end 
end

else
% x is an array

[m1,m2]=size(x);
if m1<m2
    x=x';
end
m=max(m1,m2);
p=[ones(m,1),zeros(m,1),zeros(m,1)];   pd=zeros(m,3); 
if n == 0    
   return 
elseif n==1
p(:,2) = p(:,1); p(:,3)=p(:,2);
p(:,1) = (alpha-beta+(apb+2)*x)*0.5; 

pd(:,2) = pd(:,1); pd(:,3)=pd(:,2);
pd(:,1) = 0.5*(apb+2); 
else
p(:,2) = p(:,1); p(:,3)=p(:,2);
p(:,1) = (alpha-beta+(apb+2)*x)*0.5; 

pd(:,2) = pd(:,1); pd(:,3)=pd(:,2);
pd(:,1) = 0.5*(apb+2);     
for k = 1:n-1 
    k2=k*2; k2ab=k2+alpha+beta;
   p(:,3)=p(:,2); p(:,2)=p(:,1); 
   pd(:,3)=pd(:,2); pd(:,2)=pd(:,1); 
   a1 = 2*(k+1)*(k+apb+1)*k2ab; 
% Gamma(x+n)/Gamma(x) = (x+n-1) * ... * (x+1) * x 
   a21 = (k2ab+1)*ab2;
   a22 = (k2ab+2)*(k2ab+1)*k2ab; 
   a3=2*(k+alpha)*(k+beta)*(k2ab+2);
   p(:,1) = ((a21+a22*x).*p(:,2)-a3*p(:,3))/a1; 
   pd(:,1)= (a22*p(:,2)+(a21+a22*x).*pd(:,2)-a3*pd(:,3))/a1; 
end 
end

end
return 
 
 

