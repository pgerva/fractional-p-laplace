function Iex=exact_integral_value(s,p)
% EXACT_INTEGRAL_VALUE: evaluation of the exact value of 
%   P.V. \int_R  g_x(y) dy 
%   with g_x(y) = abs(u(x)-u(y))^(p-2)*(u(x)-u(y))/abs(x-y)^(1+sp)
%   and with  u(x)=(1-|x|^m)_+^s, m=p/(p-1)
%   It holds: (-\Delta_p)^s u(x) =  C_{s,p} * P.V. \int_R  g_x(y) dy  
%   with
%   C_{s,p}=s*p/2*(1-s)*2^(2*s-1)*Gamma((1+s*p)/2)/
%                                           (Gamma((p+1)/2)*Gamma(2-s))
%   Iex=exact_integral_value(s,p)
%   Input: s = fractional order 
%          p = index of the p-Laplace operator
%   Output: Iex = P.V. \int_R  g_x(y) dy 
%
%   Reference: F. Colasuonno, F. Ferrari, P. Gervasio, and A. Quarteroni.
%   Some evaluations of the fractional p-Laplace operator on radial
%   functions. (2021)
%   Preprint: https://arxiv.org/abs/2112.08239
%
%   Copyright (C) 2007 - 2022  by Paola Gervasio


if p==2
Iex=pi/sin(pi*s);

elseif p==3

  if abs(s-2/15)<eps
    % s=2/15
    inte01=-5/2+5*gamma(11/30)*gamma(17/15)/(2^(4/15)*sqrt(pi))-2/3*pi*sec(7*pi/30);
  elseif abs(s-1/5)<eps
    % s=1/5
    inte01=1/3*(-5-2*sqrt(2-2/sqrt(5))*pi-4*gamma(-2/5)*gamma(6/5)/gamma(4/5));
  elseif abs(s-0.4)<eps
    % s=0.4
    inte01=1/6*(-5+4*gamma(-4/5)*(-2*gamma(7/5)/gamma(3/5)+gamma(9/5)));
  elseif abs(s-1/2)<eps
    % s=1/2
    inte01=2/3*(-1+log(4));
  elseif abs(s-3/4)<eps
    % s=3/4
    inte01=2/9*(-2+3*pi-8*sqrt(pi)*gamma(7/4)/gamma(1/4));
  elseif abs(s-7/12)<eps
    % s=7/12
    inte01=4/21*(-3+7*pi-7*gamma(-7/6)*gamma(19/12)/gamma(5/12));
  else
    warning('the exact value of the integral is not available')
    Iex=[]; return
  end
Iex=2/(s*p)+2*inte01;  

elseif p==4
  % s=2/5
  if abs(s-2/15)<eps
    inte01=0;
  elseif abs(s-1/5)<eps
    % s=1/5
    inte01= -5/4+3*pi/sqrt(2*(5+sqrt(5)))-27*(gamma(-3/5))^2/(50*gamma(4/5))-...
      9*gamma(-3/5)*gamma(11/10)/(2^(9/5)*sqrt(pi));
  elseif abs(s-0.4)<eps
    % s=0.4
    inte01=-5/8-3*pi/sqrt(10-2*sqrt(5))+54*(gamma(-6/5))^2/(125*gamma(3/5))-...
      9*gamma(-6/5)*gamma(7/5)/(4*gamma(1/5));
  elseif abs(s-1/2)<eps
    % s=1/2
    inte01=5/2-3/4*pi;
  elseif abs(s-3/4)<eps
    % s=3/4
    inte01=1/12*(-4+9*sqrt(2)*pi+6*(2+sqrt(2))*sqrt(pi)*gamma(-5/4)/gamma(-3/4));
  elseif abs(s-7/12)<eps
    % s=7/12
    inte01=-3/28*(4+7*sqrt(2)*pi*(1-6*gamma(-7/4)/((1+sqrt(3))*gamma(-7/12)*gamma(-1/6)))+...
      42*gamma(-7/4)*gamma(7/6)/gamma(-7/12));
  else
    warning('the exact value of the intrgral is not available')
    Iex=[]; return
  end
Iex=2/(s*p)+2*inte01;
end

