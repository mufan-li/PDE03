% higher level script to run script.m

% initialize global vars
global Uno Uname BCno PDEno PDEname;

Uno = 0;
BCno = 0;
PDEno = 0;
ntimes = 5;
nodex = 2.^(2:ntimes+1);
nodey = nodex;

UnoList = 0:2;
PDEnoList = 0:4;
% UnoList=0;
% PDEnoList=0;

for Uno = UnoList
	for PDEno = PDEnoList
		script;

		disp(strcat([PDEname,', u = ',Uname]));
		disp(errg);
		errgd = errg(2:ntimes)-errg(1:ntimes-1);
		disp(errgd);
		errgr = errgd(1:ntimes-2) ./ errgd(2:ntimes-1);
		disp(errgr);
	end
end

