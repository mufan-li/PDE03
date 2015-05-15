% function [rhs, coefu, coefux, coefuxx, coefuxxx, coefuxxxx, rhsd] = pde1(x)
%
% returns the values of the PDE coefficient functions and right side

function [rhs, coefu, coefux, coefuxx, coefuxxx, coefuxxxx, rhsd, coefut] = pde1(x, t)

global PDEno PDEname Uno Uname Dim;
global T SigmaC Rf K Smin Smax;
global etaA etaB etaC nu mu R eee;

% note only take/return 1d vectors
o = ones(max(size(x),size(t))); z = zeros(max(size(x),size(t)));
% watch out!!! (default)
coefuxxx = 0; coefuxxxx = 0; rhsd = 0;

switch PDEno
% SP Problems
case {1003}
%from Ascher+M+R pg 390 ex 10.3 -- boundary layer at 0, when 0 <~ eee
%-1, 2+cos(pi*x), eee, [0, 1]
coefu = -1; coefux = 2+cos(pi*x); coefuxx = eee;
case {1012}
%from Ascher+M+R pg 442 ex 10.12 -- boundary layers at -1 and 1, when 0 <~ eee
%-0.5, -x, eee, [-1, 1]
coefu = -1/2; coefux = -x; coefuxx = eee;
case {1092}
%from Ascher+M+R pg 371 ex 9.2 -- interior schock layer at 0, when 0 <~ eee
%0, x, eee, [-1, 1]
coefu = 0; coefux = x; coefuxx = eee;
case {1004}
%from Ascher+M+R pg 394 ex 10.4 -- interior schock layer at 0, when 0 <~ eee
% 0, 2*x, eee, [-1, 1]
coefu = 0; coefux = 2*x; coefuxx = eee;

% layer problems
case {905}
%PB 5 in [320] from Carey and Dinh PB 2 -- interior layer at mu, when nu large
coefu = 0; coefux = -2*nu*(x-mu); coefuxx = -1/nu - nu*(x-mu)^2;
case {904}
%PB 4 in [320]
%adapted from Celia+Gray -- boundary layers at 0 and 1, when etaB large
coefu =-1; coefux = 1; coefuxx = 1;
case {903}
%PB 3 in [320] from Celia+Gray pg 180 -- boundary layer at 0, when etaB large
coefu = 0; coefux = etaB; coefuxx = etaA+etaB*x;
% nonlinear problem
case {906}
coefu = -1; coefux = 0; coefuxx = 1;

case{'PS2-4(iii)'}
PDEname = '0.01*Uxx - Ux = Ut';
coefu     = 0;
coefux    = -1;
coefuxx   = 0.01;
coefuxxx  = 0;
coefuxxxx = 0;
coefut    = 1;

case{'PS2-4(iv)'}
PDEname = '0.01*Uxx + Ux = Ut';
coefu     = 0;
coefux    = 1;
coefuxx   = 0.01;
coefuxxx  = 0;
coefuxxxx = 0;
coefut    = 1;

case{15}
PDEname = 'x^2*Uxx + x*Ux - U = Ut';
coefu     = -1; % x^(3/2);
coefux    = x; % sin(x) + 1;
coefuxx   = x^2; %exp(x/2);
coefuxxx  = 0; % 1/(1+x^2);
coefuxxxx = 0;
coefut    = 1;

case{14}
PDEname = 'x^2*Uxx - U = Ut';
coefu     = -1; % x^(3/2);
coefux    = 0; % sin(x) + 1;
coefuxx   = x^2; %exp(x/2);
coefuxxx  = 0; % 1/(1+x^2);
coefuxxxx = 0;
coefut    = 1;

case{13}
PDEname = 'x^2*Uxx = Ut';
coefu     = 0; % x^(3/2);
coefux    = 0; % sin(x) + 1;
coefuxx   = x^2; %exp(x/2);
coefuxxx  = 0; % 1/(1+x^2);
coefuxxxx = 0;
coefut    = 1;

case{12}
PDEname = 'x*Uxx - x*Ux + x*U = Ut';
coefu     = x; % x^(3/2);
coefux    = -x; % sin(x) + 1;
coefuxx   = x; %exp(x/2);
coefuxxx  = 0; % 1/(1+x^2);
coefuxxxx = 0;
coefut    = 1;

case{11}
PDEname = 'x*Uxx + x*U = Ut';
coefu     = x; % x^(3/2);
coefux    = 0; % sin(x) + 1;
coefuxx   = x; %exp(x/2);
coefuxxx  = 0; % 1/(1+x^2);
coefuxxxx = 0;
coefut    = 1;

case{10}
PDEname = 'x*Uxx = Ut';
coefu     = 0; % x^(3/2);
coefux    = 0; % sin(x) + 1;
coefuxx   = x; %exp(x/2);
coefuxxx  = 0; % 1/(1+x^2);
coefuxxxx = 0;
coefut    = 1;

case{2}
PDEname = 'Uxx - Ux + U = Ut';
coefu     = 1; % x^(3/2);
coefux    = -1; % sin(x) + 1;
coefuxx   = 1; %exp(x/2);
coefuxxx  = 0; % 1/(1+x^2);
coefuxxxx = 0;
coefut    = 1;

case{1}
PDEname = 'Uxx + U = Ut';
coefu     = 1; % x^(3/2);
coefux    = 0; % sin(x) + 1;
coefuxx   = 1; %exp(x/2);
coefuxxx  = 0; % 1/(1+x^2);
coefuxxxx = 0;
coefut    = 1;

case{-1}
PDEname = 'Ut-0.5*SigmaC^2*x^2*Uxx-Rf*x*Ux+Rf*U = 0';
coefu     = -Rf+z;
coefux    = +Rf.*x;
coefuxx   = +0.5.*SigmaC^2.*x.^2;
coefuxxx  = z;
coefuxxxx = z;
coefut    = o;

otherwise
PDEname = 'Uxx = Ut';    
%coefu   =-1;        % 0;                  % 0;           %-1/(x+2);%-1/(x+1);
%coefux  = 1;        %-2*nu*(x-mu);        % etaB;        % sin(x); % x^2 + 1;
%coefuxx = 1;        %-1/nu - nu*(x-mu)^2; % etaA+etaB*x; % exp(x);
coefu     = z; % x^(3/2);
coefux    = z; % sin(x) + 1;
coefuxx   = o; %exp(x/2);
coefuxxx  = z; % 1/(1+x^2);
coefuxxxx = z;
coefut    = o;
%if x > 0.5, coefu = -2;, end;
end

[true1, true2, true3, true4, true5, t6] = truevd(x, t);

if nargout < 5, coefuxxx = 0; coefuxxxx = 0; end
rhs = coefut.*t6 - (coefu.*true1 + coefux.*true2 + coefuxx.*true3 ...
    + coefuxxx.*true4 + coefuxxxx.*true5);
% derivative of p*u'' + q*u' + r*u for p, q, r constant
rhsd = coefu.*true2 + coefux.*true3 + coefuxx.*true4;
