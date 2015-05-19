% higher level script to run script.m

% initialize global vars
global Uno Uname BCno PDEno PDEname;

Uno = 0;
BCno = 0;
PDEno = 0;
ntimes = 5;
nodex = 2.^(2:ntimes+1);
nodey = nodex;

% UnoList = [-1:8]; % basic debugging
PDEnoList = 0:4; % basic
UnoList = [10:13,30:31]; % complex functions

% UnoList=0;
% PDEnoList=4;

for PDEno = PDEnoList
	for Uno = UnoList
		script;

		disp(strcat([PDEname,', u = ',Uname]));
		disp(errg);
		errgr = errg(1:ntimes-1) ./ errg(2:ntimes);
		disp(errgr);
	end
end

