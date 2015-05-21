% function [rhs, coefu, coefux, coefuxx, coefuxxx, coefuxxxx, rhsd] = pde1(x)
%
% returns the values of the PDE coefficient functions and right side

function [rhs, coefu, coefux, coefuxx, coefuy, ...
			coefuyy, coefuxy, coefut] = pde2(x, y, t)

global PDEno PDEname Uno Uname;
% global T SigmaC Rf K Smin Smax;
% global etaA etaB etaC nu mu R eee;

% note only take/return 1d vectors
s = max(size(x),max(size(y),size(t)));
o = ones(s); 
z = zeros(s);
% watch out!!! (default)
coefuxxx = 0; coefuxxxx = 0; rhsd = 0;

% PDEnoList = [0:5,10:18,20:22];

switch PDEno

case{22}
	PDEname = 't.^2 * Uxx + t.^2 * Uyy = Ut + g';  
	coefu     = z; % x^(3/2);
	coefux    = z; % sin(x) + 1;
	coefuxx   = t.^2; %exp(x/2);
	coefuy    = z;
	coefuyy	  = t.^2;
	coefuxy   = z;
	coefut    = o;

case{21}
	PDEname = 'Uxx + Uyy = t.^2 * Ut + g';  
	coefu     = z; % x^(3/2);
	coefux    = z; % sin(x) + 1;
	coefuxx   = o; %exp(x/2);
	coefuy    = z;
	coefuyy	  = o;
	coefuxy   = z;
	coefut    = t.^2;

case{20}
	PDEname = 'Uxx + Uyy = t * Ut + g';  
	coefu     = z; % x^(3/2);
	coefux    = z; % sin(x) + 1;
	coefuxx   = o; %exp(x/2);
	coefuy    = z;
	coefuyy	  = o;
	coefuxy   = z;
	coefut    = t;


case{18}
	PDEname = 'x * y * U + Uxx + Uyy = Ut + g';  
	coefu     = x.*y; % x^(3/2);
	coefux    = z; % sin(x) + 1;
	coefuxx   = o; %exp(x/2);
	coefuy    = z;
	coefuyy	  = o;
	coefuxy   = z;
	coefut    = o;

case{17}
	PDEname = 'x * y * Uxy + Uxx + Uyy = Ut + g';  
	coefu     = z; % x^(3/2);
	coefux    = z; % sin(x) + 1;
	coefuxx   = o; %exp(x/2);
	coefuy    = z;
	coefuyy	  = o;
	coefuxy   = x.*y;
	coefut    = o;

case{16}
	PDEname = 'y * Ux + x * Uy + Uxx + Uyy = Ut + g';  
	coefu     = z; % x^(3/2);
	coefux    = y; % sin(x) + 1;
	coefuxx   = o; %exp(x/2);
	coefuy    = x;
	coefuyy	  = o;
	coefuxy   = z;
	coefut    = o;

case{15}
	PDEname = 'x * Ux + y * Uy + Uxx + Uyy = Ut + g';  
	coefu     = z; % x^(3/2);
	coefux    = x; % sin(x) + 1;
	coefuxx   = o; %exp(x/2);
	coefuy    = y;
	coefuyy	  = o;
	coefuxy   = z;
	coefut    = o;

case{14}
	PDEname = 'x * y.^2 * Uxx + x.^2 * y * Uyy = Ut + g';  
	coefu     = z; % x^(3/2);
	coefux    = z; % sin(x) + 1;
	coefuxx   = x.*y.^2; %exp(x/2);
	coefuy    = z;
	coefuyy	  = x.^2.*y;
	coefuxy   = z;
	coefut    = o;

case{13}
	PDEname = 'x^2 * y * Uxx + x * y^2 * Uyy = Ut + g';  
	coefu     = z; % x^(3/2);
	coefux    = z; % sin(x) + 1;
	coefuxx   = x.^2.*y; %exp(x/2);
	coefuy    = z;
	coefuyy	  = x.*y.^2;
	coefuxy   = z;
	coefut    = o;

case{12}
	PDEname = 'x * y * Uxx + x * y * Uyy = Ut + g';  
	coefu     = z; % x^(3/2);
	coefux    = z; % sin(x) + 1;
	coefuxx   = x.*y; %exp(x/2);
	coefuy    = z;
	coefuyy	  = x.*y;
	coefuxy   = z;
	coefut    = o;

case{11}
	PDEname = 'y * Uxx + x * Uyy = Ut + g';  
	coefu     = z; % x^(3/2);
	coefux    = z; % sin(x) + 1;
	coefuxx   = y; %exp(x/2);
	coefuy    = z;
	coefuyy	  = x;
	coefuxy   = z;
	coefut    = o;

case{10}
	PDEname = 'x * Uxx + y * Uyy = Ut + g';  
	coefu     = z; % x^(3/2);
	coefux    = z; % sin(x) + 1;
	coefuxx   = x; %exp(x/2);
	coefuy    = z;
	coefuyy	  = y;
	coefuxy   = z;
	coefut    = o;

case{5}
	PDEname = 'U - Ux - Uy - Uxy + Uxx + Uyy = Ut + g';  
	coefu     = o; % x^(3/2);
	coefux    = -o; % sin(x) + 1;
	coefuxx   = o; %exp(x/2);
	coefuy    = -o;
	coefuyy	  = o;
	coefuxy   = -o;
	coefut    = o;

case{4}
	PDEname = 'Uxy + Uxx + Uyy = Ut + g';  
	coefu     = z; % x^(3/2);
	coefux    = z; % sin(x) + 1;
	coefuxx   = o; %exp(x/2);
	coefuy    = z;
	coefuyy	  = o;
	coefuxy   = o;
	coefut    = o;

case{3}
	PDEname = 'Uy + Uxx + Uyy = Ut + g';
	coefu     = z; % x^(3/2);
	coefux    = z; % sin(x) + 1;
	coefuxx   = o; %exp(x/2);
	coefuy    = -o;
	coefuyy	  = o;
	coefuxy   = z;
	coefut    = o;

case{2}
	PDEname = 'Ux + Uxx + Uyy = Ut + g';  
	coefu     = z; % x^(3/2);
	coefux    = -o; % sin(x) + 1;
	coefuxx   = o; %exp(x/2);
	coefuy    = z;
	coefuyy	  = o;
	coefuxy   = z;
	coefut    = o;

case{1}
	PDEname = 'U + Uxx + Uyy = Ut + g';  
	coefu     = -o; % x^(3/2);
	coefux    = z; % sin(x) + 1;
	coefuxx   = o; %exp(x/2);
	coefuy    = z;
	coefuyy	  = o;
	coefuxy   = z;
	coefut    = o;

case{0}
	PDEname = 'Uxx + Uyy = Ut + g';  
	coefu     = z; % x^(3/2);
	coefux    = z; % sin(x) + 1;
	coefuxx   = o; %exp(x/2);
	coefuy    = z;
	coefuyy	  = o;
	coefuxy   = z;
	coefut    = o;

otherwise
PDEname = 'Uxx + Uyy = Ut + g';  
coefu     = z; % x^(3/2);
coefux    = z; % sin(x) + 1;
coefuxx   = o; %exp(x/2);
coefuy    = z;
coefuyy	  = o;
coefuxy   = z;
coefut    = o;
%if x > 0.5, coefu = -2;, end;
end

[t1,t2,t3,t4,t5,t6,t7] = truevd2(x, y, t);

if nargout < 5, coefuxxx = 0; coefuxxxx = 0; end
rhs = (coefu.*t1 + coefux.*t2 + coefuxx.*t3 ...
    + coefuy.*t4 + coefuyy.*t5 + coefuxy.*t6 - coefut.*t7);
% derivative of p*u'' + q*u' + r*u for p, q, r constant
% rhsd = coefu.*true2 + coefux.*true3 + coefuxx.*true4;






