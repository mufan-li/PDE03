
% classdef solution < handle
%	class defined for handling solution related 
%	parameters and functions
%
classdef solution < handle

properties
	opt_pde % struct for returning change in opt.pde settings

	uj0 % solution at previous step
	uj1 % solution at current (next) step
	uj00 % initial condition
	aux % auxiliary variable (operator splitting)
	P % penalty matrix
	f % exercise value, only consider constant throughout time
	ucomp % stored solutions (can be empty)

	rhs0 % right hand side of PDE, extra terms in BC
	rhs1
	coefs0 % PDE coefficients of derivatives
	coefs1
	b0 % extra terms from Dirichlet BC
	b1

	A0
	A1
	A0d
	A1d
	A0b
	A1b
	A0n
	A1n
	Am
	Im

	theta
	Aim
	Aex
	rhs
end

methods
	
	% 
	% function [s, opt_pde] = solution(pgrid, opt_pde, opt_var, StoreU)
	%	- initialize a grid object
	%
	% inputs:
	%	pgrid - grid object, containing grid required for 
	%			boundary conditions
	%	opt_pde, opt_var - structs, containing input options
	%	StoreU - boolean, decides whether to store solution
	%
	% outputs:
	%	s - a solution object, with initialized values
	%	opt_pde - pde options, note problem names are assigned here
	%
	function [s] = solution(pgrid, opt_pde, opt_var, StoreU)
		% initial condition
		[s.uj0,s.opt_pde] = DirichletBC(pgrid.gridx, pgrid.gridy, ...
						 pgrid.tj, opt_pde, opt_var);

		s.f = s.uj0; % exercise value

		s.uj1 = s.uj0; % for initialization
		s.uj00 = s.uj0;
		s.aux = zeros(size(s.uj0));

		% store data, includes IC
		if StoreU
			s.ucomp = zeros(length(s.uj0),pgrid.nt*3+1);
			s.ucomp(:,1) = s.uj0;
		end
		
		% P = spdiag(aux);
	end % end function solution

	% function setup_mat(s, pgrid, opt_pde, opt_var)
	% 	- set values in matrices before solve the linear system
	%		uj1 = Aim\(Aex*uj0 - rhs);
	%
	% inputs:
	%	s - solution class, containing the matrices to be set
	%	pgrid - grid class
	%	opt_pde - struct, PDE related options
	%	opt_var - struct, variable related options
	%
	% changes: (no output)
	%	s.rhs0 - vector, RHS of the linear system at t=tj0
	%	s.coefs0 - matrix, coefficients in the PDE at t=tj0
	%	s.b0 - vector, extra values due to BC at t=tj0
	%
	%	s.A0d - matrix, derivative related terms
	%			Ad = Im * (Ad0+Adx+Adxx+Ady+Adyy+Adxy); 
	%	s.A0b - matrix, Dirichlet BC terms
	%			Ab = Ic - Im - In;
	%	s.A0n - matrix, Neumann BC terms (Heston)
	%			An = Anu + Anx + Anxx + Any + Anyy;
	%	s.A0 - matrix, combining with BC terms
	%			A = Ab + Ad + An
	%	s.Im - matrix, non-BC coordinates
	%	s.Am - struct of matrices, all extra matrices above
	%
	%	s.theta - float in [0,1], Crank-Nicolson parameter
	%	s.Aim - matrix, the implicit component of 
	%			the linear system to be solved
	%	s.Aex - matrix, the explicit component
	%	s.rhs - vector, the RHS terms of the linear system
	%	
	function setup_mat(s, pgrid, opt_pde, opt_var)
		% assume time dependent coefficients
		[s.rhs0, s.coefs0, s.b0] = ...
			rhscfd2(pgrid, pgrid.tj, opt_pde, opt_var);
		[s.rhs1, s.coefs1, s.b1] = ...
			rhscfd2(pgrid, pgrid.tj1, opt_pde, opt_var);

		[s.A0, s.A0d, s.A0b, s.A0n] = ...
			cfd2(s.coefs0, pgrid, opt_pde, opt_var);
			% cfd2(nx, ny, gridx, gridy, coefs0, bet);
		[s.A1, s.A1d, s.A1b, s.A1n, s.Am] = ...
			cfd2(s.coefs0, pgrid, opt_pde, opt_var);
		s.Im = s.Am.Im;

		% Rannacher Smoothing
		%	choose theta = 1 (implicit) stepping for j <= 3
		s.theta = max(ismember(pgrid.stepj,[1:3]), pgrid.theta); 

		s.Aim = s.Im - s.theta * pgrid.htj *(s.A1d + s.A1n) + s.A1b;
		s.Aex = s.Im + (1 - s.theta) * pgrid.htj * (s.A0d + s.A0n);
		s.rhs = pgrid.htj * (s.theta * s.rhs1 + ...
				(1 - s.theta) * s.rhs0) - s.b1;
	end % end function setup_mat

	% function [opt_pen] = solve_mat(s, pgrid, opt_pde, opt_pen)
	function [opt_pen] = solve_mat(s, pgrid, opt_pde, opt_pen)
		% global Penalty opt_pen.Name OptionType tol;

		m = length(s.uj0); 
		mx = length(pgrid.gridx); 
		my = length(pgrid.gridy);

		aux0 = s.aux;
		aux1 = aux0; % unassigned
		pgrid.nitr = pgrid.nitr+1; % only penalty is different
		nitk = 0; % number of LU iterations 

		s.P = 0*spdiag(s.f); % init

		switch opt_pde.Amer
		case {1} % if the problem is LCP and requires penalty

			% choose penalty type
			switch opt_pen.no;

			% % Perform poorly
			% case {5,'Quadratic Penalty'}
			% 	% Zvan, Forsyth, Vetzal (1998)
			% 	opt_pen.Name = 'Quadratic Penalty';
			% 	% tol = 1e-6;
			% 	ujk0 = uj0;
			% 	P = spdiag((uj0-f)<0)/tol;

			% 	for k = 1:10
			% 		ujk1 = (Aim-P*spdiag(ujk0-2*f,m)) \ ...
			% 			(Aex*uj0 - rhs + P*f.^2);
			% 		P0 = P;
			% 		P = spdiag((ujk1-f)<0)/tol;

			% 		if (max( ... % component-wise division
			% 			abs(ujk1-ujk0) ./ max(ones(size(ujk1)),ujk1) ...
			% 		) < tol ) ...
			% 			|| isequal(P0,P)
			% 			break
			% 		end
			% 		ujk0 = ujk1;
			% 	end
			% 	uj1 = ujk1;
			% 	nit = nit+k-1; % additional iterations
			% 	nitk = k;

			% % Equivalent to case 2
			% case {4,'Splitting HaHo13'}
			% 	opt_pen.Name = 'Splitting HaHo13';
			% 	vj1 = Aim\(Aex*uj0 - rhs + htj*aux0);
			% 	aux1 = max( aux0 - 1/htj*(vj1 - f) ,0);
			% 	uj1 = max( vj1-htj*aux0 ,f);

			case {3,'Discrete Penalty'}
				opt_pen.Name = 'Discrete Penalty';
				
				ujk0 = s.uj0;
				s.P = spdiag((s.uj0 - s.f)<0)/opt_pen.tol;

				for k = 1:10
					ujk1 = (s.Aim + s.P) \ ...
							(s.Aex * s.uj0 - s.rhs + s.P * s.f);
					P0 = s.P;
					s.P = spdiag((ujk1 - s.f) < 0) / opt_pen.tol;

					if (max( ... % component-wise division
						abs(ujk1-ujk0) ./ max(ones(size(ujk1)),ujk1) ...
					) < opt_pen.tol ) ...
						|| isequal(P0,s.P)
						break
					end
					ujk0 = ujk1;
				end
				s.uj1 = ujk1;
				pgrid.nitr = pgrid.nitr+k-1; % additional iterations
				nitk = k;

			case {2,'Operator Splitting'}
				opt_pen.Name = 'Operator Splitting';
				vj1 = s.Aim \ (s.Aex * s.uj0 - s.rhs + pgrid.htj * aux0);
				aux1 = s.Im * max( aux0 - 1/pgrid.htj * (vj1 - s.f) ,0);
				s.uj1 = max( vj1 + pgrid.htj * (aux1 - aux0), s.f);
				s.aux = aux1;

			case {1,'Explicit Penalty'}
				opt_pen.Name = 'Explicit Penalty';
				s.uj1 = s.Aim \ (s.Aex * s.uj0 - s.rhs);
				ind = s.uj1 < s.f;
				s.uj1(ind) = s.f(ind);

			otherwise
				opt_pen.Name = 'None';
				s.uj1 = s.Aim \ (s.Aex * s.uj0 - s.rhs);
			end % end switch opt_pen.no
			
		otherwise
			opt_pen.Name = 'None';
			s.uj1 = s.Aim \ (s.Aex * s.uj0 - s.rhs);
		end % end switch opt_pde.Amer

		pgrid.nitr_gt(pgrid.stepj) = nitk;

	end % end function solve_mat

end % end methods
end % end classdef












