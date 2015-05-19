% New Script file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables init outside of script.m %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uno
% BCno
% ntimes
% nodex
% nodey
% nintx

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions need changes %
%%%%%%%%%%%%%%%%%%%%%%%%%%

% rhscfd2
% cfd2
% errorfd
% pde1

BC = 'DirechletBC';
Dim = 2;
Uname = '';
% ax = Smin; bx = Smax;
ax = 0; bx = 1;
ay = 0; by = 1; % y dim
errg = zeros(1, ntimes);

for ni = 1:ntimes
	%%%%%%%%%%%%%%%%%%%
	% init space vars %
	%%%%%%%%%%%%%%%%%%%

	nn = ni; % modify later for flexible node counts

	nx = nodex(ni); % for specific node counts
    ngridx = nx+1; nintx(nn) = nx;
    neqx = nx-1; 
    numeqx = neqx;
    hx = (bx-ax)/nx;
    gridx = ax + hx*[0:nx]; gridxu = gridx;
    % gridx = nugrid(gridxu,ax,bx); %non-uniform spacing

    ny = nodey(ni); % for specific node counts
    ngridy = ny+1; ninty(nn) = ny;
    neqy = ny-1; 
    numeqy = neqy;	% n-1 for D, n for DN and periodic, n+1 for N
    hy = (by-ay)/ny;
    gridy = ay + hy*[0:ny]; gridyu = gridy;
    % gridy = nugrid(gridyu,ay,by); %non-uniform spacing

    %%%%%%%%%%%%
    % time dim %
    %%%%%%%%%%%%

    % to be added

    %%%%%%%%%%%%%%%%%%%%%%%
    % Solve linear system %
    %%%%%%%%%%%%%%%%%%%%%%%

	[rhs, coefs] = rhscfd2(nx, ny, gridx, gridy);
    [A, A2, A1, A0] = cfd2(nx, ny, gridx, gridy, coefs);
	uj1 = A\(rhs);

	% Calculate error
    errg = errorfd(ngridx, ngridy, gridx, gridy, ni, uj1, ...
    				errg, truevd);
end




















