
% classdef grid < handle
%	class defined for handling grid related 
%	parameters and functions
%
classdef grid < handle

properties
	nx % number of nodes
	neqx % number of equations, also matrix size
	hx % unif increment size
	gridxu % unif grid
	gridx % general grid

	ny
	neqy
	hy
	gridyu
	gridy

	nt
	nitr % number of iterations, each solves a linear system
	ngridt % = nt + 1
	stepj % current step number
	ht % unif time step
	gridt % time grid
	nitr_gt % number of iterations at each time step
	htj % current time step size
	tj % current time value
	tj1 % next time value
	dnorm % adapt time param
	theta % Crank-Nicolson coefficient
end

methods

	% function g = grid(ni, opt_grid, x_center, y_center, penaltyno)
	%	initialize the grid class
	%
	% inputs:
	%	ni: int, experiment iteration count
	%	opt_grid: struct, containing grid options
	%	x_center: float, center of x-non-unif grid
	%	y_center: float, center of y-non-unif grid
	%	penaltyno: int, choice of penalty type
	%
	% outputs:
	%	g: grid class
	%
	function g = grid(ni, opt_grid, x_center, y_center, penaltyno)
		g.nx = opt_grid.nodex(ni); % for specific node counts
		g.neqx = g.nx+1; % num of eq = matrix size
		
		% unif grid increment
		g.hx = (opt_grid.xmax - opt_grid.xmin) / g.nx; 
		g.gridxu = opt_grid.xmin + g.hx * [0 : g.nx]; % uniform grid

		g.ny = opt_grid.nodey(ni); % for specific node counts
		g.neqy = g.ny + 1; 
		g.hy = (opt_grid.ymax - opt_grid.ymin) / g.ny;
		g.gridyu = opt_grid.ymin + g.hy * [0 : g.ny];

		% init non-uniform grids
		[g.gridx, g.gridy] = g.nugrid(opt_grid, x_center, y_center);

		g.nt = opt_grid.nodet(ni);
		g.nitr = 0; % total # of iterations
		g.ngridt = g.nt + 1;
		g.stepj = 0;
		g.ht = (opt_grid.tmax - opt_grid.tmin) / g.nt; 
		% initialize larger for non-uniform stepping
		g.gridt = opt_grid.tmin + g.ht * [0 : g.nt * 3]; 
		g.nitr_gt = zeros(size(g.gridt));

		g.htj = g.ht;
		g.tj = g.gridt(1); 
		g.tj1 = g.tj; % for init
		g.dnorm = opt_grid.dnorm0(penaltyno) / 2^(ni-1); % half each time
		g.theta = opt_grid.theta;
	end % end function grid

	% [gridx, gridy] = nugrid(g, opt_grid, x_center, y_center)
	%	creates non-uniform grids for both x and y direction
	%
	function [gridx, gridy] = nugrid(g, opt_grid, x_center, y_center)
		
		gridx = g.each_grid(g.gridxu, opt_grid.xmin, opt_grid.xmax, ...
							opt_grid.xno(opt_grid.NonUnif + 1), ...
							x_center);

		gridy = g.each_grid(g.gridyu, opt_grid.ymin, opt_grid.ymax, ...
							opt_grid.yno(opt_grid.NonUnif + 1), ...
							y_center);

	end % end function nugrid

	% function [gridx] = each_grid(g, gridxu, xmin, xmax, ...
	%						 Gridno, x_center)
	%	creates non-uniform grid for one axis
	%
	function [gridx] = each_grid(g, gridxu, xmin, xmax, Gridno, x_center)
		m = length(gridxu)-1;
		
		switch Gridno
			case {31}
				% centered strike
				a0 = 0.2;
				a = grid21a(gridxu,x_center,xmax,a0);
				gridx = grid21(gridxu,x_center,xmax,a0);

			case {30}
				% volatility grid for the Heston model
				d2 = xmax/50;
				gridx = d2 * sinh(asinh(xmax/d2)/m * (0:m));

			case {21}
				% centered strike
				a0 = 20;
				a = grid21a(gridxu,x_center,xmax,a0);
				gridx = grid21(gridxu,x_center,xmax,a0);
				
			case {20}
				%Black Scholes Grid 1
				a = (ceil(0.38*m)+0.5)/m;
				b = fzero(@(z) sinh(z*(1-a))/sinh(z*a) - ...
							(xmax-x_center)/x_center,10);
				gridx = (1 + sinh(b.*(gridxu/xmax-a)) / ...
							sinh(b*a) )*x_center;

			otherwise
				gridx = gridxu;
		end % end switch

		gridx(1) = 0; % floating point error causing x<0

	end % end function each_grid

	% function nugridt(g, opt_grid, uj1, uj0, Rt)
	%  - modifies the step size of the next step
	%  - change in place, tj, tj1, htj, gridt
	%
	% inputs:
	% 	g.htj - the current step size, to be modified
	%	g.stepj - the current step number
	%	g.tj - the current time value, to be modified
	%	g.tj1 - the next time value, to be modified, to be modified
	%	g.gridt - the time step grid, to be modified
	%	g.dnorm - adaptive time step coefficient
	%	Rt - adaptive time step coefficient - cap
	%	s - the solution object
	%	s.uj1 - the solution at time tj1 (before modifying)
	%	s.uj0 - the solution at time tj (before modifying)
	%
	% changes (no output):
	%	g.htj - the current time step size
	%	g.tj1 - the next time point
	%
	% note:
	%  - grid contains htj, stepj, pt (tj), gridt, dnorm
	%  - opt_grid contains tmax (T), adapt_time, Rt
	%
	function nugridt(g, s, opt_grid)
		% move forward one step
        g.tj = g.tj1;

		switch opt_grid.adapt_time
			case {1}
				% simple non-uniform
				if (g.stepj==1)
					g.htj = g.htj/1e3;
				elseif (g.stepj>2)
					% element wise
					ht1 = g.dnorm * g.htj * ...
						min(max(max(opt_grid.Rt, s.uj0), s.uj1) ./ ...
						abs(s.uj0-s.uj1));
					% dnorm = dnorm / 2^(ht1>ht*2);
					g.htj = min(ht1, g.htj*3); % restricting growth
					% ht = ht1;
				end
				
				% check tmax
				if (g.tj + g.htj > opt_grid.tmax)
					g.htj = opt_grid.tmax - g.tj;
				end
				
			otherwise
				% uniform ht
				g.htj = g.htj; % for explicitness
		end

		% update tj1, gridt
		g.tj1 = g.gridt(g.stepj) + g.htj;
		g.gridt(g.stepj+1) = g.tj1;

	end % end function nugridt

end % end methods
end % end classdef

% function [a] = grid21a(gridxu,K,Smax,a0)
%	determine the a that (almost) centers the strike
%	in between two nodes
%
function [a] = grid21a(gridxu,K,Smax,a0)
	ng = 2e3; % # of guesses
	plk.x = zeros(ng+1,1);
	plk.y = zeros(ng+1,2);
	
	for i = 1:ng+1
		a = a0 * (1+(i-ng/2-1)/(ng/2*10));% restrict to 10%
		gridx = grid21(gridxu,K,Smax,a);
		nk = find(max(gridx-K,0),1);
		yk = abs(gridx(nk-1:nk)-K);
		plk.x(i) = a/a0;
		plk.y(i,:) = yk;
	end
	
	[~,ik] = sort(abs(plk.y(:,1)-plk.y(:,2)));
	[~,ik2] = min(abs(plk.x(ik(1:10))-1));
	a = plk.x(ik(ik2))*a0;
%	 plot(plk.x,plk.y);
%	 disp(a);
end

% function [gridx] = grid21(gridxu,K,Smax,a)
%	Black-Scholes non-uniform grid 2
%
function [gridx] = grid21(gridxu,K,Smax,a)
%	 a = 4; % a=4 sets grid 21 ~ grid 20
	c1 = asinh((Smax-K)/a);
	c2 = asinh((0-K)/a);
	gridx = K + a*sinh(c1*gridxu/Smax+c2*(1-gridxu/Smax));
end








