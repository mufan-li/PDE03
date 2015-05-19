% format compact
% global Uno Uname BCno PDEno PDEname BC UxStep Unif Unift GridC Dim theta rho;
% global n0 ntimes convx convt IC_Offset dnormC Penalty nodex nodet;
% global T SigmaC Rf K Smin Smax;

% the following parameters are set outside of script.m
% -------------------------------------------------------- %
% Uno = -1;%set to zero (-1)
% BCno = 'European Call';%'PS2-4(iii)';
% PDEno = -1;%'PS2-4(iii)';
% rho = 0.5;
% theta = 0.5;    
% UxStep = 'forward';%['centered','forward','backward']
% n0 = 1;
% ntimes = 2;	% number of grid sizes to run
% convx = 1;% set to 0 to turn off hx size loop
% convt = 1;% set to 0 to turn off ht size loop
% GridC = 0;
% Unift = 1;
% dnormC = 2e-1;
% 
% T = 1;
% SigmaC = 0.2;
% Rf = log(1.03);
% K = 10;
% Smin = 0;
% Smax = 50;
% -------------------------------------------------------- %

% Changed less frequently
BC = 'DirechletBC';
Unif = 1; % 1 = non-unif
Dim = 2;
Uname = '';
ax = Smin; bx = Smax;
errg = zeros(1, ntimes);
tol = 1e-6; % tolerance for penalty

% global etaA etaB etaC nu mu R eee logSpaceC;
% global nint errg gridx gridxu; % for plotting
etaA = 1; etaB = 1;
mu = 0.5*(ax+bx); nu = 100;
etaC = 1e-2; R = -1/4;
eee = 1e-4;
% logSpaceC = 1;

% [udummy] = truevd(ax,0);
% disp(['U = ' Uname ' = {' num2str(Uno) '}'])
% disp(['domain [' num2str(ax)  ', ' num2str(bx) ']'])

conv_i = zeros(8,ntimes);

% plot fbt
figure;
hold on;

for ni = 1:ntimes
    
    nn = ni;
    % if convx==0
    %     nn = ntimes+n0+2;
    % end
    
%     n = 2^(nn+n0+1);
    n = nodex(ni); % for specific node counts
    ngrid = n+1; nint(nn) = n;
    neq = n-1; numeq = neq;	% n-1 for D, n for DN and periodic, n+1 for N
    hx = (bx-ax)/n;
    gridx = ax + hx*[0:n]; gridxu = gridx;
    gridx = nugrid(gridxu,ax,bx); %non-uniform spacing
    
    % set ht % non-unif adjusts to the initial setting of nt
    pt0 = 0;
    ptn = T;
%     ht = ((T)/2^(ni+n0))^(2^convx) * rho; % convx slows down ht size change
    ht = T/nodet(ni); % for specific node counts
    % if convt==0
    %     ht = T/2^(ntimes*1+n0+2);
    % end
    nt = ceil((ptn - pt0) / ht); % round to int
    ht = (T-pt0)/nt; % set again to line up with the end
    nt = nt*(1+Unift*2); % increase array size for non-unif gridt
    pt = pt0;

    % initial condition
    % changed to using BC
    uj0 = reshape(DirechletBC(gridx(2:numeq+1),pt),numeq,1);
    % dummy
    uj1 = zeros(numeq,1);
    % store
    ucomp = zeros(numeq,nt);
    
    % for non-unif grid
    gridt = zeros(1,nt);
    dnorm = dnormC*2^(-ni); % non-unif ht
    Rt = 1; % non-unif ht
    
    % for iterative methods
    nstep = 0;
    
    % saving the coef array to avoid calculation - indep of t
    coefs_in = zeros(neq,3);
    [~, coefs_in(:, 1), coefs_in(:, 2), coefs_in(:, 3)] = pde1(gridx(2:neq+1),0);
    
    % auxilary term for operator splitting method
    aux = 0; 
    auxt = zeros(neq,nt); % store aux
    pendiagt = zeros(neq,nt); % store diagonal of penalty matrix
    f = reshape(DirechletBC(gridx(2:neq+1),0),neq,1); % exercise value
    fbt = zeros(nt,1); % store fb
    
    for stepj = 1:nt
        % if not the initial step then assign uj0 = uj1
        if (pt ~= pt0)
            uj0=uj1;
        end
        
        % non-unif ht % note @stepj=1,pt=0,ht is used here although added to pt after        
        [ht, gridt, dnorm] = nugridt(ht, stepj, pt, T, ucomp, gridt, dnorm, Rt);
        
        if (stepj<= 2 + Unift*5) % Rannacher smoothing, theta = 1
            [rhs, coefs, coefs1] = rhscfd2(n, gridx, pt, ht, 1, coefs_in);
            [A, B, A2, A1, A0, B2, B1, B0] = cfd2(n, gridx, coefs, coefs1, pt, ht, 1);
        else
            [rhs, coefs, coefs1] = rhscfd2(n, gridx, pt, ht, theta, coefs_in);
            [A, B, A2, A1, A0, B2, B1, B0] = cfd2(n, gridx, coefs, coefs1, pt, ht, theta);
        end
        % function for penalty/other iterations
        [uj1, nstep, aux, P] = t_step(uj0, tol, A, B, rhs, nstep, ht, aux);

        auxt(:,stepj) = aux;
        pendiagt(:,stepj) = diag(full(P));
        ucomp(:,stepj) = uj1;
        fbnx = find(max(uj1 - (f+tol*1),0),1); % finds first non-zero
        fbt(stepj) = gridx(fbnx+1);

        if abs(pt + ht - T) < 1e-10
            break;
        else
            pt = pt + ht; % j+=1
        end
    end
    
%     disp([ni,stepj,nt,pt,dnorm])
    nt = stepj; % reset to make later calculations easier
    gridt = gridt(1:stepj+1);
    auxt = auxt(:,1:stepj);
    pendiagt = pendiagt(:,1:stepj);
    fbt = fbt(1:stepj);

    % plot
    plot(gridt(2:nt+1),fbt);
    lgv{nn} = num2str(nn);

%     figure;
%     plot(gridt(2:stepj)-gridt(1:stepj-1));
    
    ucomp = ucomp(1:numeq,1:nt);
    % calculate the average abs error at each time step
    error_i = 0;
    for stepi = nt
%         error_i = error_i + abs(ucomp(stepi,:)-truevd(gridxu(1+stepi),pt0+ht*(1:nt)));
        [EuroCall, EuroPut] = EuroBls(gridx(2:neq+1), K, Rf, T, 0, SigmaC);
        switch BCno
            case {'European Call'}
                TrueVal = EuroCall;
            case {'European Put'}
                TrueVal = EuroPut;
            case {'European Call + delta'}
                TrueVal = EuroBls(gridx(2:neq+1), K, Rf, T+IC_Offset, 0, SigmaC);
            case {'European Put + delta'}
                [ans,TrueVal] = EuroBls(gridx(2:neq+1), K, Rf, T+IC_Offset, 0, SigmaC);
            case {'American Put', 'American Call'}
                TrueVal = EuroCall;
            otherwise
                disp('True value default to Euro Call');
                TrueVal = EuroCall;
        end
        error_temp = abs(ucomp(:,nt)-reshape(TrueVal,[neq 1]));
        error_i = max([error_i; error_temp]); % truevd=0
    end
    
    %save/plot 4(ii)
    [uval, fb, delta, gamma] = rsummary(ucomp,gridx,neq,K,nt,tol);
    conv_i(:,ni) = [n;nt;max(error_i);nstep;uval;fb;delta;gamma];
    
%     % plot a few parts of the solution
%     figure;
%     stepj = [1 4 nt/4 nt/2 3*nt/4 15*nt/16 nt];
%     [null0, coefu, coefux, coefuxx, coefuxxx, coefuxxxx, null1, coefut] = pde1(0,0);
%     plot(gridx(2:n),ucomp(:,stepj(1)),'k',...
%         gridx(2:n),ucomp(:,stepj(2)),'b',...
%         gridx(2:n),ucomp(:,stepj(3)),'r',...
%         gridx(2:n),ucomp(:,stepj(4)),'o',...
%         gridx(2:n),ucomp(:,stepj(5)),'.',...
%         gridx(2:n),ucomp(:,stepj(6)),'k.',...
%         gridx(2:n),ucomp(:,stepj(7)),'r.'...
%         )
% %     axis([-1 1 -1.5 1.5]);
%     legend(strcat('step =',num2str(stepj(1))),...
%         strcat('step =',num2str(stepj(2))),...
%         strcat('step =',num2str(stepj(3))),...
%         strcat('step =',num2str(stepj(4))),...
%         strcat('step =',num2str(stepj(5))),...
%         strcat('step =',num2str(stepj(6))),...
%         strcat('step =',num2str(stepj(7))),...
%         'Location','NorthWest')
%     title(strcat('PDE: ', PDEname,', theta = ', num2str(theta), ', rho = ', num2str(ht/hx^2)));
%     xlabel(strcat('x, ht = ', num2str(ht), ', hx = ', num2str(hx), ', Ux Step = ',UxStep));
%     ylabel('u(x,t)');
    
    %plot 4(i)
%     figure;
%     plot((1:nt)*ht,error_i/neq);
%     [null0, coefu, coefux, coefuxx, coefuxxx, coefuxxxx, null1, coefut] = pde1(0,0);
%     title(strcat('Uno = ', num2str(Uno) ,', u = ', Uname, ...
%         ', theta = ', num2str(theta), ', rho = ', num2str(ht/hx^2)));
%     xlabel(strcat('t, ht = ', num2str(ht), ', hx = ', num2str(hx)));
%     ylabel('average squared error');
%     text((ptn+pt0)/5, max(error_i/neq)*0.9, ...
%         strcat( 'PDE: ', PDEname)...
%         , 'Color', 'k');

    % Note: pt is at j, results/error is at j+1 -> pt+ht
%     errg = errorfd(ngrid, gridx, pt, n, uj1, nn, errg);
end

% plot fbt
legend(lgv);
hold off;

% nint
% format short e
% disp(' ')
% disp('error on grid points')
% errg
% format short
% disp('order of convergence')
% errg(:, :) = max(errg(:, :),  0.222044604925e-15);
% LogNintRatio = log(nint(1, 2:ntimes)./nint(1, 1:ntimes-1));
% LogNintRatioMat = repmat(LogNintRatio, size(errg, 1), 1);
% if ntimes > 1
%     convg = log(errg(:, 1:ntimes-1)./errg(:, 2:ntimes))./LogNintRatioMat
% end
