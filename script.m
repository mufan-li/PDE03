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
ay = ymin; by = ymax; % y dim
at = 0; bt = T; % t dim (IVP)
theta = 1/2; % Crank-Nicolson
errg = zeros(1, ntimes);
m = summary(ntimes); % summary class

for ni = 1:ntimes
    tic
	%%%%%%%%%%%%%%%%%%%
	% init space vars %
	%%%%%%%%%%%%%%%%%%%

	nn = ni; % modify later for flexible node counts

	nx = nodex(ni); % for specific node counts
    ngridx = nx+1; %nintx(nn) = nx;
    neqx = nx+1; 
    numeqx = neqx;
    hx = (bx-ax)/nx;
    gridx = ax + hx*[0:nx]; gridxu = gridx;
    gridx = nugrid(gridxu,ax,bx,Gridno); %non-uniform spacing

    ny = nodey(ni); % for specific node counts
    ngridy = ny+1; %ninty(nn) = ny;
    neqy = ny+1; 
    numeqy = neqy;	% n-1 for D, n for DN and periodic, n+1 for N
    hy = (by-ay)/ny;
    gridy = ay + hy*[0:ny]; gridyu = gridy;
    gridy = nugrid(gridyu,ay,by,Gridno); %non-uniform spacing

    %%%%%%%%%%%%
    % time dim %
    %%%%%%%%%%%%

    nt = nodet(ni);
    nit = 0; % total # of iterations
    ngridt = nt + 1;
    ht = (bt-at)/nt; htj = ht;
    gridt = at + ht*[0:nt*3]; % initialize larger
    tj = gridt(1); tj1 = tj; % for init
    dnorm = dnorm0 / 2^(ni-1); % half each time
    % note - adaptive grid changed inside loop

    uj0 = DirechletBC(gridx,gridy,at); % IC
    uj1 = uj0; % for initialization
    uj00 = uj0;
    aux = zeros(size(uj0));
    % P = spdiag(aux);
    %%%%%%%%%%%%%%%%%%%%%%%
    % Solve linear system %
    %%%%%%%%%%%%%%%%%%%%%%%

    for j = 1:nt
        tj = tj1;
        [htj, gridt] = nugridt(htj, j, tj, T, uj1, uj0,...
                                gridt, dnorm, 1);
        tj1 = tj + htj;

        uj0 = uj1; % after adaptive time step

        % assume time dependent coefficients
        [rhs0, coefs0, b0] = rhscfd2(nx, ny, gridx, gridy, tj);
        [rhs1, coefs1, b1] = rhscfd2(nx, ny, gridx, gridy, tj1);

        [A0,A0d,A0b] = cfd2(nx, ny, gridx, gridy, coefs0);
        [A1,A1d,A1b,Am] = cfd2(nx, ny, gridx, gridy, coefs1);
        Im = Am.Im;

        Aim = Im - theta*htj*A1d + A1b;
        Aex = Im + (1-theta)*htj*A0d;
        rhs = htj*(theta*rhs1 + (1-theta)*rhs0) - b1;

        [uj1,aux,nit,P] = t_step(uj0, rhs, Aim, Aex, ...
                    gridx, gridy, aux, htj, tj1, nit);

        if (tj1>=T)
            break
        end
    end
	nt = j;
    gridt = gridt(1:j+1);

    Nm.nx = nx+1; % grid points
    Nm.ny = ny+1;
    Nm.nt = nt; % steps
    Nm.nit = nit;
    Nm.ni = ni;

    Gm.gx = gridx; Gm.gy = gridy; % grids
    Gm.x = xp; Gm.y = yp; % points to evaluate

    % storing summary
    update(m,uj1,Am,Nm,Gm);
end

print(m);
plot(m,uj1,Gm);
% plot_greeks(m,uj1,Gm,Am);
% disp(EuroRb(xp,yp,T))



