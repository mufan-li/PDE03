% higher level script to run script.m

% initialize global vars
global Uno Uname BCno PDEno PDEname;

Uno = 6;
BCno = 0;
PDEno = 0;
ntimes = 5;
nodex = 2.^(2:ntimes+1);
nodey = nodex;

% UnoList = [0:8]; % basic debugging
% PDEnoList = [0:5]; % basic
UnoList = [10:18,30:31]; % complex functions
PDEnoList = [10:18]; % complex PDEs
% UnoList=0;
% PDEnoList=2;

for PDEno = PDEnoList
	for Uno = UnoList
		script;

		disp(strcat([PDEname,', u = ',Uname]));
		if (max(abs(errg))<1e-10)
			disp('Solution Exact.');
		end
		
		disp(errg);
		errgr = errg(1:ntimes-1) ./ errg(2:ntimes);
		disp(errgr);
	end
end

