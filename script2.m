% higher level script to run script.m

% initialize global vars
global Uno Uname BCno PDEno PDEname Rbno RbName;

Uno = 6; BCno = 0; PDEno = 0; ntimes = 4;
nodex = 2.^(2:ntimes+1); nodey = nodex; nodet = nodex;

global T Sx Sy rho Rf K Smin Smax Penalty PenaltyName;
% init variables
T = 2; Sx = 0.05; Sy = 0.2; rho = -0.8; Rf = 0.05; K = 30;
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
RbnoList = [10:12]; % Euro 0:2, Amer 10:12
PenaltyList = 1:3;

for PDEno = PDEnoList
	for Uno = UnoList
		for Rbno = RbnoList
			for Penalty = PenaltyList
				script;
			end
		end
	end
end