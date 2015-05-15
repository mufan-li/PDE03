function [uj1, nstep, aux, P] = t_step(uj0, tol, A, B, rhs, nstep, ht, aux)
    global BCno Penalty gridx;
    m = length(uj0);
    P = spdiag(-1,m); % can't be equal

    switch BCno;
        case {'American Call', 'American Put'}
            switch Penalty
                case {1, 'Discrete'}
                    m = length(uj0);
                    ujk0 = uj0;
                    f = reshape(DirechletBC(gridx(2:m+1),0),m,1);
                    P0 = spdiag(-1,m); % can't be equal
                    P = spdiag((uj0 < (f-tol*0))) / tol; % initialize P
%                     sum(uj0 < (f-tol*-1))
                    for k = 1:10
                        ujk1 = (A+P)\(B*uj0 + rhs + P*f);

                        P0 = P; % construct P for stopping
                        P = spdiag((ujk1 < (f-tol*0))) / tol; 
                        % stop criterion - scaled
                        if (max(abs(ujk1 - ujk0)) / ...
                                max([1;ujk1])) < tol || isequal(P0,P)
                            break
                        end
                        ujk0 = ujk1; % next step
                    end
                    uj1 = ujk1;
                    nstep = nstep + k;
                    
%                 case {'Continuous'}
                    
                case {2, 'Operator Splitting'}
                    m = length(uj0);
                    aux0 = aux;
                    % using last step info
                    vj1 = A\(B*uj0 + rhs + aux0*ht); 
                    % handles call/put
                    f = reshape(DirechletBC(gridx(2:m+1),0),m,1); 
                    aux1 = max(aux0 - (vj1 - f)/ht,0);
                    aux = aux1;
                    uj1 = max(vj1 - ht*(aux0 - aux1), f);
                    nstep = nstep + 1;
                    
                otherwise
                    disp('No penalty');
                    uj1 = A\(B*uj0 + rhs);
                    nstep = nstep + 1;
            end
            
        otherwise
            uj1 = A\(B*uj0 + rhs);
            nstep = nstep + 1;
    end
end
