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
Uname = '';
ax = Smin; bx = Smax;
ay = ymin; by = ymax; % y dim
at = 0; bt = T; % t dim (IVP)
theta0 = 1/2; % Crank-Nicolson
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
    gridx = nugrid(gridxu,ax,bx,Gridxno(Gdno+1),K); %non-uniform spacing

    ny = nodey(ni); % for specific node counts
    ngridy = ny+1; %ninty(nn) = ny;
    neqy = ny+1; 
    numeqy = neqy;	% n-1 for D, n for DN and periodic, n+1 for N
    hy = (by-ay)/ny;
    gridy = ay + hy*[0:ny]; gridyu = gridy;
    gridy = nugrid(gridyu,ay,by,Gridyno(Gdno+1),bet); %non-uniform spacing

    %%%%%%%%%%%%
    % time dim %
    %%%%%%%%%%%%

    nt = nodet(ni);
    nit = 0; % total # of iterations
    ngridt = nt + 1;
    ht = (bt-at)/nt; htj = ht;
    gridt = at + ht*[0:nt*3]; % initialize larger
    nitv = zeros(size(gridt));
    tj = gridt(1); tj1 = tj; % for init
    dnorm = dnorm0(Penalty) / 2^(ni-1); % half each time
    % note - adaptive grid changed inside loop

    uj0 = DirechletBC(gridx,gridy,at); % IC
    uj1 = uj0; % for initialization
    uj00 = uj0;
    aux = zeros(size(uj0));

    [~, coefs00, ~] = rhscfd2(nx, ny, gridx, gridy, ...
                                0, 0);
    % P = spdiag(aux);
    %%%%%%%%%%%%%%%%%%%%%%%
    % Solve linear system %
    %%%%%%%%%%%%%%%%%%%%%%%

    for j = 1:nt*3
        tj = tj1;
        [htj, gridt] = nugridt(htj, j, tj, T, uj1, uj0,...
                                gridt, dnorm, 1);
        tj1 = tj + htj;

        uj0 = uj1; % after adaptive time step

        % assume time dependent coefficients
        [rhs0, coefs0, b0] = rhscfd2(nx, ny, gridx, gridy, ...
                                tj, coefs00);
        [rhs1, coefs1, b1] = rhscfd2(nx, ny, gridx, gridy, ...
                                tj1, coefs00);

        [A0,A0d,A0b,A0n] = cfd2(nx, ny, gridx, gridy, ...
                                coefs0,bet);
        [A1,A1d,A1b,A1n,Am] = cfd2(nx, ny, gridx, gridy, ...
                                coefs1,bet);
        Im = Am.Im;

        theta = max(ismember(j,1:3),theta0); % Rannacher Smoothing
        Aim = Im - theta*htj*(A1d+A1n) + A1b;
        Aex = Im + (1-theta)*htj*(A0d+A0n);
        rhs = htj*(theta*rhs1 + (1-theta)*rhs0) - b1;
        % Aim\(Aex*uj0 - rhs)

        [uj1,aux,nit,P,f,nitk] = t_step(uj0, rhs, Aim, Aex, ...
                    gridx, gridy, aux, htj, tj1, nit, Im);
        nitv(j) = nitk;

        if (tj1>=T)
            break
        end
    end
	nt = j;
    gridt = gridt(1:j+1);
    nitv = nitv(1:j+1);

    Nm.nx = nx+1; % grid points
    Nm.ny = ny+1;
    Nm.nt = nt; % steps
    Nm.nit = nit;
    Nm.ni = ni;

    Gm.gx = gridx; Gm.gy = gridy; Gm.gt = gridt; % grids
    Gm.x = xp; Gm.y = yp; % points to evaluate

    % storing summary
    update(m,uj1,Am,Nm,Gm,Aim);

    if TrackTime
        toc
    end
end

print(m,Display);
if (Display)
    plot(m,uj1,Gm);
    % plot(m,spdiags(P,0),Gm);
    % plot(m,aux,Gm);
    plot_greeks(m,uj1,Gm,Am);
    % plot_greeks_csx(m,uj1,Gm,Am,[Gm.gy(2),yp/3,2*yp/3,yp,3*yp/2]);
end
% disp(EuroRb(xp,yp,T))



