% higher level script to run script.m

opt.pde.Uno = -1;
opt.pde.BCno = 1;
opt.pde.PDEno = 101;
% IC/BC choices, Euro 0:2, Amer 10:19, Heston 30:31
opt.pde.Rbno = 31;
% will be reassigned in DirichletBC.m
opt.pde.Amer = true;
opt.pde.RbName = 'Heston American Put';

opt.var.T = 0.25;
opt.var.sigmax = 0.2; % vol of stock x
opt.var.sigmay = 0.2; % vol of stock y
opt.var.rho = 0.1; % correlation
opt.var.Rf = 0.02; % risk-free
opt.var.K = 10; % strike
opt.var.q1 = 0; % dividend of x
opt.var.q2 = 0; % dividend of y
% Heston Notation following Ikonen & Toivanen (2009)
opt.var.alpha = 5; % mean-revert coef
opt.var.beta = 0.16; % mean-revert level
opt.var.gamma = 0.9; % vol coefficient

opt.grid.ntimes = 1;
opt.grid.NonUnif = false; % toggles non-unif init
opt.grid.xno = [0 21]; % options of both unif and non-unif
opt.grid.yno = [0 31];
opt.grid.nodex = 2.^((1:opt.grid.ntimes)+2); % number of nodes
opt.grid.nodey = opt.grid.nodex;
opt.grid.nodet = opt.grid.nodex;
opt.grid.xmin = 0;
opt.grid.xmax = 15 * opt.var.K;
opt.grid.ymin = 0;
opt.grid.ymax = 5;
opt.grid.tmin = 0;
opt.grid.tmax = opt.var.T;
opt.grid.adapt_time = true;
opt.grid.dnorm0 = [5 1.3 5 5 5]; % adapt time step parameter
opt.grid.Rt = 1; % adapt time step parameter
opt.grid.theta = 1/2; % Crank-Nicolson

% penalty parameters
opt.pen.tol = 1e-6;
opt.pen.no = 2;
opt.pen.Name = 'None';

% summary control
opt.sum.xp = opt.var.K; % point of evaluation
opt.sum.yp = 0.25;
opt.sum.xv = 8:12 * 10 / opt.var.K; % multiple points
opt.sum.yv = [.0625 .25];
opt.sum.TrackTime = true;
opt.sum.Display = true;
opt.sum.StoreU = 1 * (opt.grid.ntimes <= 4);

script;

% % initialize global vars
% global Uno Uname BCno PDEno PDEname Rbno RbName Gdno;
% global T Sx Sy rho Rf K q1 q2 xp yp;
% global alp bet gam;
% global Smin Smax ymax Unift Penalty PenaltyName OptionType tol;

% Uno = -1; BCno = 1; PDEno = 101;
% ntimes = 3;
% Gdno = 0;
% Gridxno = [0 21]; Gridyno = [0 31];
% % nodex = 2.^((1:ntimes)+2); nodey = nodex; nodet = nodex;
% nodex = [1 2 3]*100; nodey = nodex/2; nodet = nodex/4;

% % init variables
% % T = 1; Sx = 0.2; Sy = Sx; rho = 0.5; Rf = 0.05; K = 100;
% T = 0.25; Sx = 0.2; Sy = Sx; rho = 0.1; Rf = 0.1; K = 10;
% q1 = 0.0; q2 = q1; % dividend
% % heston
% alp = 5; bet = 0.16; gam = 0.9;

% xp = K; yp = 0.25; % points to evaluate
% xv = 8:12; yv = [0.0625 0.25];
% Smin = 0; Smax = 14*K; ymin = 0; ymax = 5;
% dnorm0 = [5 1.3 5 5 5]; tol = 1e-6;
% Unift = 1; TrackTime = 1; Display = 1;
% StoreU = 1 * (ntimes<=4);

% % European Rainbow
% UnoList=-1; % zero 
% PDEno = 101;
% % RbnoList = [0,10]; % Euro 0:2, Amer 10:19, Heston 30:31
% RbnoList = 31;
% PenaltyList = 2:3;

% for Rbno = RbnoList
% 	for Penalty = PenaltyList
% 		for Gdno = 1
% 			for Unift = 0:1
% 				script;
% 			end
% 		end
% 	end
% end