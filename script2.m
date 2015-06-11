% higher level script to run script.m

% initialize global vars
global Uno Uname BCno PDEno PDEname Rbno RbName;

Uno = 6; BCno = 0; PDEno = 0;
ntimes = 4;
Gridno = 21;
nodex = 2.^((1:ntimes)+2); nodey = nodex; nodet = nodex;

global T Sx Sy rho Rf K q1 q2 xp yp;
global Smin Smax Penalty PenaltyName OptionType;
% init variables
T = 3; Sx = 0.2; Sy = Sx; rho = 0; Rf = 0.05; K = 100;
q1 = 0.1; q2 = q1; % dividend
xp = K; yp = xp; % points to evaluate
Smin = 0; Smax = 500;

% UnoList = [0:8]; % basic debugging
% PDEnoList = [0:5]; % basic
% UnoList = [10:18,30:31]; % complex functions
% PDEnoList = [10:18]; % complex PDEs
% UnoList = [40:43,50]; % IVP specific functions
% PDEnoList = [20:25]; % time dependent coef PDEs

% European Rainbow
UnoList=-1; % zero 
PDEnoList=100;
RbnoList = [12]; % Euro 0:2, Amer 10:19
% RbnoList = 0;
PenaltyList = 2:3;

for PDEno = PDEnoList
	for Uno = UnoList
		for Rbno = RbnoList
			for Penalty = PenaltyList
				script;
			end
		end
	end
end