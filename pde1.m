% function [rhs, coefu, coefux, coefuxx, coefuxxx, coefuxxxx, rhsd] = pde1(x)
%
% returns the values of the PDE coefficient functions and right side

function [rhs, coefu, coefux, coefuxx, coefuy, coefuyy, coefuxy] = pde1(x, y)

global PDEno PDEname Uno Uname;
% global T SigmaC Rf K Smin Smax;
% global etaA etaB etaC nu mu R eee;

% note only take/return 1d vectors
o = ones(max(size(x),size(y))); z = zeros(max(size(x),size(y)));
% watch out!!! (default)
coefuxxx = 0; coefuxxxx = 0; rhsd = 0;

% PDEnoVec = [0:4];

switch PDEno
case{4}
	PDEname = 'Uxy + Uxx + Uyy = 0';  
	coefu     = z; % x^(3/2);
	coefux    = z; % sin(x) + 1;
	coefuxx   = o; %exp(x/2);
	coefuy    = z;
	coefuyy	  = o;
	coefuxy   = o;

case{3}
	PDEname = 'Uy + Uxx + Uyy = 0';  
	coefu     = z; % x^(3/2);
	coefux    = z; % sin(x) + 1;
	coefuxx   = o; %exp(x/2);
	coefuy    = o;
	coefuyy	  = o;
	coefuxy   = z;

case{2}
	PDEname = 'Ux + Uxx + Uyy = 0';  
	coefu     = z; % x^(3/2);
	coefux    = o; % sin(x) + 1;
	coefuxx   = o; %exp(x/2);
	coefuy    = z;
	coefuyy	  = o;
	coefuxy   = z;

case{1}
	PDEname = 'U + Uxx + Uyy = 0';  
	coefu     = o; % x^(3/2);
	coefux    = z; % sin(x) + 1;
	coefuxx   = o; %exp(x/2);
	coefuy    = z;
	coefuyy	  = o;
	coefuxy   = z;

case{0}
	PDEname = 'Uxx + Uyy = 0';  
	coefu     = z; % x^(3/2);
	coefux    = z; % sin(x) + 1;
	coefuxx   = o; %exp(x/2);
	coefuy    = z;
	coefuyy	  = o;
	coefuxy   = z;

otherwise
PDEname = 'Uxx + Uyy = 0';  
coefu     = z; % x^(3/2);
coefux    = z; % sin(x) + 1;
coefuxx   = o; %exp(x/2);
coefuy    = z;
coefuyy	  = o;
coefuxy   = z;
%if x > 0.5, coefu = -2;, end;
end

[true1, true2, true3, true4, true5, t6] = truevd(x, t);

if nargout < 5, coefuxxx = 0; coefuxxxx = 0; end
rhs = coefut.*t6 - (coefu.*true1 + coefux.*true2 + coefuxx.*true3 ...
    + coefuxxx.*true4 + coefuxxxx.*true5);
% derivative of p*u'' + q*u' + r*u for p, q, r constant
rhsd = coefu.*true2 + coefux.*true3 + coefuxx.*true4;
