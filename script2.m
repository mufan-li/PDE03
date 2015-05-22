% higher level script to run script.m

% initialize global vars
global Uno Uname BCno PDEno PDEname;

Uno = 6;
BCno = 0;
PDEno = 0;
ntimes = 3;
nodex = 2.^(2:ntimes+1);
nodey = nodex;
nodet = nodex;

UnoList = [0:8]; % basic debugging
% PDEnoList = [0:5]; % basic
% UnoList = [10:18,30:31]; % complex functions
% PDEnoList = [10:18]; % complex PDEs
% UnoList = [40:43,50]; % IVP specific functions
PDEnoList = [20:25]; % time dependent coef PDEs
% UnoList=1;
% PDEnoList=22;

for PDEno = PDEnoList
	for Uno = UnoList
		script;

		disp(strcat([PDEname,', u = ',Uname]));
		if (max(abs(errg))<1e-10)
			disp('Solution Exact.');
			disp(' ');
		else
			disp(errg);
			errgr = errg(1:ntimes-1) ./ errg(2:ntimes);
			disp(errgr);
		end
	end
end

