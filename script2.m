% higher level script to run script.m

% initialize global vars
global Uno Uname BCno PDEno PDEname Rbno RbName Gridno;

Uno = -1; BCno = 1; PDEno = 101;
ntimes = 3;
Gridno = 0;
nodex = 2.^((1:ntimes)+2); nodey = nodex; nodet = nodex;

global T Sx Sy rho Rf K q1 q2 xp yp;
global alp bet gam;
global Smin Smax Unift Penalty PenaltyName OptionType;

% init variables - two asset
T = 1; Sx = 0.2; Sy = Sx; rho = 0.2; Rf = 0.05; K = 100;
q1 = 0; q2 = q1; % dividend
% heston
alp = 1; bet = Sy^2; gam = Sy;

xp = K; yp = 0.2; % points to evaluate
Smin = 0; Smax = 500;
ymin = 0; ymax = 1;
Unift = 1; dnorm0 = 0.05;

% European Rainbow
UnoList=-1; % zero 
PDEnoList=101;
% RbnoList = [0,10]; % Euro 0:2, Amer 10:19, Heston 30:31
RbnoList = 30;
PenaltyList = 3;

for Rbno = RbnoList
	for Penalty = PenaltyList
		for Gridno = [0]
			for Unift = 0
				script;
			end
		end
	end
end