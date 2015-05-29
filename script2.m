% higher level script to run script.m

% initialize global vars
global Uno Uname BCno PDEno PDEname Rbno;

Uno = 6; BCno = 0; PDEno = 0; ntimes = 5;
nodex = 2.^(2:ntimes+1); nodey = nodex; nodet = nodex;

global T Sx Sy rho Rf K Smin Smax Penalty;
% init variables
T = 1; Sx = 0.3; Sy = 0.4; rho = 0.3; Rf = 0.05; K = 10;
Smin = 0; Smax = 100;

% UnoList = [0:8]; % basic debugging
% PDEnoList = [0:5]; % basic
% UnoList = [10:18,30:31]; % complex functions
% PDEnoList = [10:18]; % complex PDEs
% UnoList = [40:43,50]; % IVP specific functions
% PDEnoList = [20:25]; % time dependent coef PDEs

% European Rainbow
UnoList=-1; % zero 
PDEnoList=100;
Rbno = 10; % margrabe
Penalty = 0;

for PDEno = PDEnoList
	for Uno = UnoList
		script;
	end
end