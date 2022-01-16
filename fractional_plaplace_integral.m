% 
% FRACTIONAL_PLAPLACE_INTEGRAL: MATLAB script for the numerical valuation of 
%   P.V. \int_R  g_x(y) dy 
%   with g_x(y) = abs(u(x)-u(y))^(p-2)*(u(x)-u(y))/abs(x-y)^(1+sp)
%   and with  u(x)=(1-|x|^m)_+^s, m=p/(p-1)
%   It holds: (-\Delta_p)^s u(x) =  C_{s,p} * P.V. \int_R  g_x(y) dy  
%   with
%   C_{s,p}=s*p/2*(1-s)*2^(2*s-1)*Gamma((1+s*p)/2)/
%                                           (Gamma((p+1)/2)*Gamma(2-s))
% 
%   x in a interval, p>1, s in (0,1)
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

% ----------------------------------------------
% Set the parameter p
p=3;
% Set an array of values of the parameter s (in (0,1)) 
S=[2/15, 0.2 , 0.4, 0.5, 7/12] ; 
% ----------------------------------------------
 
np=101;
% set of the points at which the integral I^(s,p) is  computed
xx=linspace(-1,1,np);

% The array I will contain the values of I^(s,p)(x)
% the rows span the points x, the columns span the index s
I=zeros(np,length(S));

% The array fractional_plaplace will contain the values of (-Delta_p)^s u(x)
% the rows span the points x, the columns span the index s
fractional_plaplace=zeros(np,length(S));

figure(1); clf


for ks=1:length(S)
  s=S(ks);
  m=p/(p-1);
  
  % costant C_(N,s,p) defined in [1], page 1  
  C_sp=@(N,s,p)s*p/2*(1-s)*2^(2*s-1)/pi^((N-1)/2)*...
    gamma((N+s*p)/2)/(gamma((p+1)/2)*gamma(2-s));
  % C_(1,s,p)
  Csp=C_sp(1,s,p);
  
  % definition of the function u(x)
  u=@(x)(1-abs(x).^m).^s.*(x<1).*(x>-1);
  
  % loop on the point x in (-1,1)
  for i=2:np-1
    % evaluation of the integral I^(s,p)(x) (formula (5.1) of [1])
    % by the algorithm described in Sect 5 of [1]
    I(i,ks)=integralsp(u,s,p,xx(i));
    
    % print the computed value of I^(s,p)(x)
    fprintf('s=%f, p=%f, x= %f, I^{(s,p)}(x)=%e\n',s,p,xx(i),I(i,ks))
    
    % evaluate (-Delta_p)^s u(x) =C_{1,s,p} *I^(s,p)(x)
    fractional_plaplace(i,ks)=Csp*I(i,ks);
  end
  
  % plot the computed values of the integrals
  plot(xx(2:np-1),I(2:np-1,ks),'.','Markersize',10,...
    'Displayname',['s=',num2str(s)])
  legend('-dynamiclegend')
  hold on
  
end

grid on
legend('AutoUpdate','off')
color=colormap('lines');

for ks=1:length(S)
  s=S(ks);
  % evaluate the exact integral I^(s,p)(0)
  % when p=2, I^(s,2)(0)=pi/sin(pi*s)
  % when p\neq 2, I^(s,p) is preevaluated by Wolfram Mathematica, here we
  % only report the results for some values of s
  Iex=exact_integral_value(s,p);
  
  % plot it 
  if ~isempty(Iex)
    plot(0,Iex,'s','Color',color(ks,:),'Markersize',12)

  % evaluate the error |I_appx^(s,p)(0)-I^(s,p)(0)|
    ind0=find(xx==0);
    I_num=I(ind0,ks);
    err=abs(I_num-Iex);
    fprintf('p=%d, s=%e, error at x=0: %e\n',p,s,err);
  end

end

xlabel('$x$','Interpreter','Latex')
ylabel(['$\tilde I^{(s,',num2str(p),')}(x)$'],'Interpreter','Latex')
set(gca, 'Fontsize',12)

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if p==2
% Evaluate and plot the errors  |I_appx^(s,p)(x)-I^(s,p)(x)| at each x
% available only for p=2
  figure(2);clf
  
  for ks=1:length(S)
    s=S(ks);
    Iex=exact_integral_value(s,p);
    err=abs(I(2:np-1,ks)-Iex);
    semilogy(xx(2:np-1),err,'.','Markersize',10,'Displayname',['s=',num2str(s)])
    legend('-dynamiclegend')
    hold on
  end
  
  grid on
  xlabel('$x$','Interpreter','Latex')
  ylabel(['$|I^{(s,',num2str(p),')}(x) -\tilde I^{(s,',...
    num2str(p),')}(x)|$'],'Interpreter','Latex')
  set(gca, 'Fontsize',12)
  
end
