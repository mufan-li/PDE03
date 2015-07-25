% higher level script to run script.m

% initialize global vars
global Uno Uname BCno PDEno PDEname Rbno RbName Gdno;

Uno = -1; BCno = 1; PDEno = 101;
ntimes = 3;
Gdno = 0;
Gridxno = [0 21]; Gridyno = [0 31];
nodex = 2.^((1:ntimes)+2); nodey = nodex; nodet = nodex;

global T Sx Sy rho Rf K q1 q2 xp yp;
global alp bet gam;
global Smin Smax ymax Unift Penalty PenaltyName OptionType tol;

% init variables - two asset
T = 1; Sx = 0.2; Sy = Sx; rho = 0.5; Rf = 0.05; K = 100;
q1 = 0.0; q2 = q1; % dividend
% heston
alp = 0.2; bet = 0.5; gam = 0.9;

xp = K; yp = bet; % points to evaluate
Smin = 0; Smax = 14*K;
ymin = 0; ymax = 5;
Unift = 1;
dnorm0 = [5 1.3 5 5 5];
tol = 1e-6;
TrackTime = 1; Display = 1;

% European Rainbow
UnoList=-1; % zero 
PDEno = 101;
% RbnoList = [0,10]; % Euro 0:2, Amer 10:19, Heston 30:31
RbnoList = 31;
PenaltyList = 2;

for Rbno = RbnoList
	for Penalty = PenaltyList
		for Gdno = 1
			for Unift = 0
				script;
			end
		end
	end
end