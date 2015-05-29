% New Script file

% consider changing script.m and script2.m to functions
% to including helper functions within the same file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables init outside of script.m %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uno
% BCno
% PDEno
% ntimes
% nodex
% nodey
% nodet

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions need changes %
%%%%%%%%%%%%%%%%%%%%%%%%%%

BC = 'DirechletBC';
% Dim = 2;
Uname = '';
% ax = Smin; bx = Smax;
ax = Smin; bx = Smax;
ay = Smin; by = Smax; % y dim
at = 0; bt = T; % t dim (IVP)
theta = 1/2; % Crank-Nicolson
errg = zeros(1, ntimes);

for ni = 1:ntimes
    tic
	%%%%%%%%%%%%%%%%%%%
	% init space vars %
	%%%%%%%%%%%%%%%%%%%

	nn = ni; % modify later for flexible node counts

	nx = nodex(ni); % for specific node counts
    ngridx = nx+1; %nintx(nn) = nx;
    neqx = nx-1; 
    numeqx = neqx;
    hx = (bx-ax)/nx;
    gridx = ax + hx*[0:nx]; gridxu = gridx;
    % gridx = nugrid(gridxu,ax,bx,1); %non-uniform spacing

    ny = nodey(ni); % for specific node counts
    ngridy = ny+1; %ninty(nn) = ny;
    neqy = ny-1; 
    numeqy = neqy;	% n-1 for D, n for DN and periodic, n+1 for N
    hy = (by-ay)/ny;
    gridy = ay + hy*[0:ny]; gridyu = gridy;
    % gridy = nugrid(gridyu,ay,by,1); %non-uniform spacing

    %%%%%%%%%%%%
    % time dim %
    %%%%%%%%%%%%

    nt = nodet(ni);
    ngridt = nt + 1;
    ht = (bt-at)/nt;
    gridt = at + ht*[0:nt];
    % note - adaptive grid changed inside loop

    uj0 = DirechletBC(gridx(2:nx),gridy(2:ny),at); % IC
    uj00 = uj0;
    %%%%%%%%%%%%%%%%%%%%%%%
    % Solve linear system %
    %%%%%%%%%%%%%%%%%%%%%%%

    for j = 1:nt
        tj = gridt(j);
        tj1 = gridt(j+1);
        htj = tj1 - tj; % adaptive time step

        % assume time dependent coefficients
        [rhs0, coefs0] = rhscfd2(nx, ny, gridx, gridy, tj);
        [rhs1, coefs1] = rhscfd2(nx, ny, gridx, gridy, tj1);

        A0 = cfd2(nx, ny, gridx, gridy, coefs0);
        A1 = cfd2(nx, ny, gridx, gridy, coefs1);
        Im = speye(size(A0));

        Aim = Im - theta*htj*A1;
        Aex = Im + (1-theta)*htj*A0;
        rhs = htj*(theta*rhs1 + (1-theta)*rhs0);

        uj1 = Aim\(Aex*uj0-rhs); % note rhs here is not on the rhs
        uj0 = uj1; % move to next step
    end
	
	% Calculate error - at the final step
    [errg,trueval] = errorfd2(ngridx, ngridy, gridx, gridy, ...
        ni, uj1, errg, bt);
    toc
end

