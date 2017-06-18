% New Script file

m = summary(opt.grid.ntimes); % summary class

for ni = 1:opt.grid.ntimes
    tic
	
    %%%%%%%%%%%%%%%%%%%
    % init space vars %
    %%%%%%%%%%%%%%%%%%%

    % physical grid
    % function grid(ni, opt_grid, x_center, y_center, penaltyno)
    pgrid = grid(ni, opt.grid, opt.var.K, opt.var.beta, opt.pen.no);

    [psoln] = solution(pgrid, opt.pde, opt.var, opt.sum.StoreU);
    % copy setting changes
    opt.pde = psoln.opt_pde; 

    %%%%%%%%%%%%%%%%%%%%%%%
    % Solve linear system %
    %%%%%%%%%%%%%%%%%%%%%%%

    for j = 1:pgrid.nt*3
        
        pgrid.stepj = j;
        % adaptive time step
        nugridt(pgrid, psoln, opt.grid);

        % copy after adaptive time step
        psoln.uj0 = psoln.uj1; 

        % set all matrix values
        setup_mat(psoln, pgrid, opt.pde, opt.var);

        % Aim\(Aex*uj0 - rhs)
        [opt_pen] = solve_mat(psoln, pgrid, opt.pde, opt.pen);

        % store values
        if opt.sum.StoreU; psoln.ucomp(:,j+1) = psoln.uj1; end;

        % reached terminal time
        if (pgrid.tj1 >= opt.grid.tmax); break; end;
    end

	pgrid.nt = j;
    pgrid.gridt = pgrid.gridt(1:j+1);
    pgrid.nitr_gt = pgrid.nitr_gt(1:j+1);

    if opt.sum.StoreU; psoln.ucomp = psoln.ucomp(:,1:j+1); end;

    % Nm.nx = nx+1; % grid points
    % Nm.ny = ny+1;
    % Nm.nt = nt; % steps
    % Nm.nit = nit;
    % Nm.ni = ni;

    % Gm.gx = gridx; Gm.gy = gridy; Gm.gt = gridt; % grids
    % Gm.x = xp; Gm.y = yp; % points to evaluate
    % Gm.xv = xv; Gm.yv = yv;

    % storing summary
    % update(m,uj1,Am,Nm,Gm,Aim);
    update(m, ni, psoln, pgrid, opt.pen.tol, opt.sum);

    if opt.sum.TrackTime; toc; end;
end

print(m, opt.sum.Display, opt);
if (opt.sum.Display)
    plot(m, psoln.uj1, pgrid);
    % plot(m,spdiags(P,0),Gm);
    % plot(m,aux,Gm);
    plot_greeks(m,psoln,pgrid,opt.grid,2*K,1.5);
    % plot_greeks_csx(m,uj1,Gm,Am,[Gm.gy(2),yp/3,2*yp/3,yp,3*yp/2]);
    plot_fb(m,psoln.uj1,pgrid,opt_pen,1.5*K,1);

    if StoreU; mesh_fb(m,ucomp,Gm); end;
end
% disp(EuroRb(xp,yp,T))



