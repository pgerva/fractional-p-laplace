# Fractional p-Laplace calculator

## A MATLAB tool for evaluating the fractional p-laplace operator in 1d domains.

This repository contains a MATLAB script, a MATALB app and the related
functions to numerically compute 


<img src="https://latex.codecogs.com/svg.image?I^{(s,p)}(x)=\lim_{\varepsilon&space;\to&space;0}\int_{(B_\varepsilon(x))^c}&space;\frac{|(u(x)-u(y)|^{p-2}(u(x)-u(y))}{|x-y|^{1&plus;sp}}dy" title="I^{(s,p)}(x)=\lim_{\varepsilon \to 0}\int_{(B_\varepsilon(x))^c} \frac{|(u(x)-u(y)|^{p-2}(u(x)-u(y))}{|x-y|^{1+sp}}dy" />
and <img src="https://latex.codecogs.com/svg.image?(-\Delta_p)^su(x)=C_{s,p}&space;I^{(s,p)}(x)" title="(-\Delta_p)^su(x)=C_{s,p} I^{(s,p)}(x)" />
with:
<img src="https://latex.codecogs.com/svg.image?C_{s,p}=\frac{\frac{sp}{2}(1-s)2^{2s-1}\Gamma\left(\frac{N&plus;sp}{2}\right)}{\Gamma\left(\frac{p&plus;1}{2}\right)\Gamma(2-s)}," title="C_{s,p}=\frac{\frac{sp}{2}(1-s)2^{2s-1}\Gamma\left(\frac{N+sp}{2}\right)}{\Gamma\left(\frac{p+1}{2}\right)\Gamma(2-s)}," />
<img src="https://latex.codecogs.com/svg.image?u(x)=(1-|x|^m)_&plus;^s,&space;\qquad&space;m=\frac{p}{p-1}" title="u(x)=(1-|x|^m)_+^s, \qquad m=\frac{p}{p-1}" />
<img src="https://latex.codecogs.com/svg.image?\mbox{for&space;any&space;}&space;p>1,&space;\&space;s\in(0,1),\&space;x\in(-1,1).&space;" title="\mbox{for any } p>1, \ s\in(0,1),\ x\in(-1,1). " />


The integrals are approximated by Gauss-Legendre quadrature formulas (see Reference  [2]) and the MATLAB function integral.m following the strategy given
in Reference [1].

[1] F. Colasuonno, F. Ferrari, P. Gervasio, and A. Quarteroni. "Some evaluations of the fractional p-Laplace operator on radial functions" (2021).  Preprint: https://arxiv.org/abs/2112.08239

[2] C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,  "Spectral Methods. Fundamentals in Single Domains"  Springer Verlag, Berlin Heidelberg New York, 2006.



## List of files

<b>fractional_plaplace_integral.m</b> is a script that computes and plot the values of `I^{(s,p)}(x)` for `p=3` and `s`  in `{2/15, 0.2 , 0.4, 0.5,
7/12}` at 99 points in `(-1,1)`.

The values of `s` and `p` can be modified by editing the script.

The values of `I^{(s,p}(x)` are computed and stored in the variable I.

The values of `(-\Delta_p)^s u(x)` are computed and stored in the variable fractional_plaplace.

The exact value of `I^{(s,p}(x)` is furnished for `p=2` at each `x` and for any `s`, and only for some values of `s` and at `x=0` when `p=1` or
`p=4`.

<b>integralsp.m</b> is a MATLAB function that computes `I^{(s,p)}(x)` for given `p`, `s`, `u` and `x`. The function `u` can be modified by the user.

<b>app.mlapp</b> is a MATLAB app that computes `I^{(s,p)}(x)` for given `p`, `s`, and `x` when `u(x)=(1-|x|^m)_+^s`,  and `m=p/(p-1)`




## Examples of usage 

You can digit the following commands in a MATLAB Command Window

```bash
>> fractional_plaplace_integral  
```

```bash
>> s=0.5; p=3; x=0.2; m=p/(p-1); u=@(x)(1-abs(x).^m).^s.*(x<1).*(x>-1);
>> Isp=integralsp(u,s,p,x)
```

```bash
>> app
```

