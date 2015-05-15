% function [rhs, coefs] = rhscfd2(n, gridx)
function [rhs, coefs, coefs1] = rhscfd2(n, gridx, pt, ht, theta, coefs_in)

% Evaluating Boundary Conditions
global BC BCno Unif Dim UxStep Uno;
switch Unif
    case{1} % non-uniform
        switch BC
            case {'DirechletBC'}
                
                % using BC without truevd             
                numeq = n-1;	% assume Dirichlet conditions
                m = numeq;	% simplicity has charm
                
                hx = gridx(2) - gridx(1);	% assume uniform grid
                % non-uniform grid
                h = gridx(2:n+1) - gridx(1:n); 
                h0 = h(1:m);
                h1 = h(2:m+1);
                
                rho = ht/hx^2;
                rhs = zeros(m,1);
                rhs0 = zeros(m,1);
                rhs1 = zeros(m,1);
                coefs = zeros(m,3);
                coefs1 = zeros(m,3);

                % switch BCno % save time when coef indep of t, truevd = 0
                %     case {-1}
                %         % calculate once before time loop
                %         coefs = coefs_in;
                %         coefs1 = coefs_in;
                %     otherwise
                        % still need coefs
                px = gridx(2:m+1);
                [rhs0, coefs(:, 1), coefs(:, 2), coefs(:, 3)] = pde1(px,pt);
                [rhs1, coefs1(:, 1), coefs1(:, 2), coefs1(:, 3)] = pde1(px,pt+ht);
                % end
                
                % combine the 2 vectors
                rhs = rhs0 * (1-theta) * ht + rhs1 * theta * ht;
                rhs = reshape(rhs, m, 1);
                
                % setup UxStep extra term % default to 'centered'
                rhs_ux_ex1 = - (1-theta)*coefs(1,2)*ht/(h0(1)*h1(1)*(h0(1)+h1(1)))*h1(1)^2;
                rhs_ux_im1 = - (theta)*coefs1(1,2)*ht/(h0(1)*h1(1)*(h0(1)+h1(1)))*h1(1)^2;
                rhs_ux_exm = + (1-theta)*coefs(m,2)*ht/(h0(m)*h1(m)*(h0(m)+h1(m)))*h0(m)^2;
                rhs_ux_imm = + (theta)*coefs1(m,2)*ht/(h0(m)*h1(m)*(h0(m)+h1(m)))*h0(m)^2;
                
                % Let pt be the j step, since IVP starts with pt=0 at j=0
                % Using boundaries directly
                px = gridx(1);
                % explicit, time step = j
                [true1] = DirechletBC(px,pt);
                rhs(1) = rhs(1) + (1-theta)*coefs(1,3)*ht/(h0(1)*h1(1)*(h0(1)+h1(1)))*2*h1(1)*true1 ...
                    + rhs_ux_ex1*true1;
                % implicit, time step = j+1
                [true1] = DirechletBC(px,pt+ht);
                rhs(1) = rhs(1) + (theta)*coefs1(1,3)*ht/(h0(1)*h1(1)*(h0(1)+h1(1)))*2*h1(1)*true1 ...
                    + rhs_ux_im1*true1;

                px = gridx(n+1);
                % explicit, time step = j
                [true1] = DirechletBC(px,pt);
                rhs(m) = rhs(m) + (1-theta)*coefs(m,3)*ht/(h0(m)*h1(m)*(h0(m)+h1(m)))*2*h0(m)*true1 ...
                    + rhs_ux_exm*true1; 
                % implicit, time step = j+1
                [true1] = DirechletBC(px,pt+ht);
                rhs(m) = rhs(m) + (theta)*coefs1(m,3)*ht/(h0(m)*h1(m)*(h0(m)+h1(m)))*2*h0(m)*true1 ...
                    + rhs_ux_imm*true1;
                
            otherwise
                numeq = n-1;	% assume Dirichlet conditions
                m = numeq;	% simplicity has charm
                h = gridx(2:n+1) - gridx(1:n);
                h0 = h(1:m);
                h1 = h(2:m+1);
                rhs = zeros(m,1);
                coefs = zeros(m,3);

                    for i = 1:m
                        px = gridx(i+1);
                        [rhs(i, 1), coefs(i, 1), coefs(i, 2), coefs(i, 3)] = pde1(px);
                    end

                rhs = (rhs - coefs(1:m,1).*reshape(truevd(gridx(2:m+1)),m,1)) .* ...
                        reshape((h0.*h1.*(h0+h1)),m,1) + ...
                        coefs(1:m,1).*reshape(truevd(gridx(2:m+1)),m,1);
                % Using boundaries directly

                px = gridx(1);
                [true1, true2, true3] = truevd(px);
                rhs(1) = rhs(1) - coefs(1, 3)*true1*2*h1(1) + coefs(1, 2)*true1*(h1(1)^2);
                px = gridx(n+1);
                [true1, true2, true3] = truevd(px);
                rhs(m) = rhs(m) - coefs(m, 3)*true1*2*h0(m) - coefs(m, 2)*true1*(h0(m)^2);
        end
    otherwise
        switch BC

            case{4}
                numeq = n-1;	% assume Dirichlet conditions
                m = numeq;	% simplicity has charm
                hx = gridx(2) - gridx(1);	% assume uniform grid
                rhs = zeros(m,1);
                coefs = zeros(m,3);

                    for i = 1:m
                        px = gridx(i+1);
                        [rhs(i, 1), coefs(i, 1), coefs(i, 2), coefs(i, 3)] = pde1(px);
                    end

                % Using boundaries directly
                px = gridx(1);
                [true1, true2, true3] = truevd(px);
                rhs(1) = rhs(1) - coefs(1, 3)*true1/hx^2*24/(65) + coefs(1, 2)*true1/(2*hx);
                px = gridx(2);
                [true1, true2, true3] = truevd(px);
                rhs(2) = rhs(2) + coefs(2, 3)*true1/hx^2/12;
                
                px = gridx(n);
                [true1, true2, true3] = truevd(px);
                rhs(m-1) = rhs(m-1) + coefs(m-1, 3)*true1/hx^2/12;
                px = gridx(n+1);
                [true1, true2, true3] = truevd(px);
                rhs(m) = rhs(m) + coefs(m, 3)*true1/hx^2*24/(65) - coefs(m, 2)*true1/(2*hx);
                
            case{3}
                % assume General Robin conditions
                numeq = n+1; % include both BC
                m = numeq;	% simplicity has charm
                hx = gridx(2) - gridx(1);	% assume uniform grid
                rhs = zeros(m,1);
                coefs = zeros(m,3);

                % modified to include BC
                for i = 1:m
                    px = gridx(i);
                    [rhs(i, 1), coefs(i, 1), coefs(i, 2), coefs(i, 3)] = pde1(px);
                end

                % First Order Differences
                px = gridx(1);
                [coef_a,coef_b,coef_f] = robinBC(px);
                rhs(1) = coef_f/hx^2;
                if (coefs(1,2) ~= 0)
                    rhs(1) = rhs(1) - coef_f/coef_b;
                end

                px = gridx(n+1);
                [coef_a,coef_b,coef_f] = robinBC(px);
                rhs(m) = coef_f/hx^2;
                if (coefs(m,2) ~= 0)
                    rhs(m) = rhs(m) - coef_f/coef_b;
                end

            case{2}
                % assume General Robin conditions
                numeq = n+1; % include both BC
                m = numeq;	% simplicity has charm
                hx = gridx(2) - gridx(1);	% assume uniform grid
                rhs = zeros(m,1);
                coefs = zeros(m,3);

                % modified to include BC
                for i = 1:m
                    px = gridx(i);
                    [rhs(i, 1), coefs(i, 1), coefs(i, 2), coefs(i, 3)] = pde1(px);
                end

                % Second Order Differences
                px = gridx(1);
                [coef_a,coef_b,coef_f] = robinBC(px);
                rhs(1) = coef_f/hx^2;
                if (coefs(1,2) ~= 0)
                    rhs(1) = rhs(1) - coef_f/coef_b;
                end

                px = gridx(n+1);
                [coef_a,coef_b,coef_f] = robinBC(px);
                rhs(m) = coef_f/hx^2;
                if (coefs(m,2) ~= 0)
                    rhs(m) = rhs(m) - coef_f/coef_b;
                end

            case{1}
%                 % assume General Robin conditions
%                 numeq = n+1; % include both BC
%                 m = numeq;	% simplicity has charm
%                 hx = gridx(2) - gridx(1);	% assume uniform grid
%                 rhs = zeros(m,1);
%                 coefs = zeros(m,3);
% 
%                 % modified to include BC
%                 for i = 1:m
%                     px = gridx(i);
%                     [rhs(i, 1), coefs(i, 1), coefs(i, 2), coefs(i, 3)] = pde1(px);
%                 end
% 
%                 % Centered Differences
%                 px = gridx(1);
%                 [coef_a,coef_b,coef_f] = robinBC(px);
%                 rhs(1) = rhs(1) + 2*coef_f/coef_b/hx;
%                 if (coefs(1,2) ~= 0)
%                     rhs(1) = rhs(1) - coef_f/coef_b;
%                 end
% 
%                 px = gridx(n+1);
%                 [coef_a,coef_b,coef_f] = robinBC(px);
%                 rhs(m) = rhs(m) - 2*coef_f/coef_b/hx;
%                 if (coefs(m,2) ~= 0)
%                     rhs(m) = rhs(m) - coef_f/coef_b;
%                 end

            case {'DirechletBC'}
                % using BC without truevd
                
                numeq = n-1;	% assume Dirichlet conditions
                m = numeq;	% simplicity has charm
                hx = gridx(2) - gridx(1);	% assume uniform grid
                rho = ht/hx^2;
                rhs = zeros(m,1);
                rhs0 = zeros(m,1);
                rhs1 = zeros(m,1);
                coefs = zeros(m,3);
                coefs1 = zeros(m,3);

                % still need coefs
                    for i = 1:m
                        px = gridx(i+1);
                        [rhs0(i), coefs(i, 1), coefs(i, 2), coefs(i, 3)] = pde1(px,pt);
                        [rhs1(i), coefs1(i, 1), coefs1(i, 2), coefs1(i, 3)] = pde1(px,pt+ht);
                    end

                % combine the 2 vectors
                rhs = rhs0 * (1-theta) * ht + rhs1 * theta * ht;
                
                % setup UxStep extra term
                switch UxStep
                    case{'forward'}
                        rhs_ux_ex1 = 0;
                        rhs_ux_im1 = 0;
                        rhs_ux_exm = + (1-theta)*coefs(m,2)*(ht/hx);
                        rhs_ux_imm = + (theta)*coefs1(m,2)*(ht/hx);
                    case{'backward'}
                        rhs_ux_ex1 = - (1-theta)*coefs(1,2)*(ht/hx);
                        rhs_ux_im1 = - (theta)*coefs1(1,2)*(ht/hx);
                        rhs_ux_exm = 0;
                        rhs_ux_imm = 0;
                    case{'centered'}
                        rhs_ux_ex1 = - (1-theta)*coefs(1,2)*(ht/2/hx);
                        rhs_ux_im1 = - (theta)*coefs1(1,2)*(ht/2/hx);
                        rhs_ux_exm = + (1-theta)*coefs(m,2)*(ht/2/hx);
                        rhs_ux_imm = + (theta)*coefs1(m,2)*(ht/2/hx);
                    otherwise
                        % default to 'centered'
                        rhs_ux_ex1 = - (1-theta)*coefs(1,2)*(ht/2/hx);
                        rhs_ux_im1 = - (theta)*coefs1(1,2)*(ht/2/hx);
                        rhs_ux_exm = + (1-theta)*coefs(m,2)*(ht/2/hx);
                        rhs_ux_imm = + (theta)*coefs1(m,2)*(ht/2/hx);
                end
                
                % Let pt be the j step, since IVP starts with pt=0 at j=0
                % Using boundaries directly
                px = gridx(1);
                % explicit, time step = j
                [true1] = DirechletBC(px,pt);
                rhs(1) = rhs(1) + (1-theta)*coefs(1,3)*rho*true1 + rhs_ux_ex1*true1;
                % implicit, time step = j+1
                [true1] = DirechletBC(px,pt+ht);
                rhs(1) = rhs(1) + (theta)*coefs1(1,3)*rho*true1 + rhs_ux_im1*true1;

                px = gridx(n+1);
                % explicit, time step = j
                [true1] = DirechletBC(px,pt);
                rhs(m) = rhs(m) + (1-theta)*coefs(m,3)*rho*true1 + rhs_ux_exm*true1; 
                % implicit, time step = j+1
                [true1] = DirechletBC(px,pt+ht);
                rhs(m) = rhs(m) + (theta)*coefs1(m,3)*rho*true1 + rhs_ux_imm*true1;

            otherwise

                numeq = n-1;	% assume Dirichlet conditions
                m = numeq;	% simplicity has charm
                hx = gridx(2) - gridx(1);	% assume uniform grid
                rho = ht/hx^2;
                rhs = zeros(m,1);
                rhs0 = zeros(m,1);
                rhs1 = zeros(m,1);
                coefs = zeros(m,3);
                coefs1 = zeros(m,3);

                    for i = 1:m
                        px = gridx(i+1);
                        [rhs0(i, 1), coefs(i, 1), coefs(i, 2), coefs(i, 3)] = pde1(px,pt);
                        [rhs1(i, 1), coefs1(i, 1), coefs1(i, 2), coefs1(i, 3)] = pde1(px,pt+ht);
                    end

                % combine the 2 vectors
                rhs = rhs0 * (1-theta) * ht + rhs1 * theta * ht;
                
                % Let pt be the j step, since IVP starts with pt=0 at j=0
                % Using boundaries directly
                px = gridx(1);
                % explicit, time step = j
                [true1, true2, true3, true4, true5, t6] = truevd(px,pt);
                rhs(1) = rhs(1) + (1-theta)*coefs(1,3)*rho*true1 - (1-theta)*coefs(1,2)*(ht/2/hx)*true1;
                % implicit, time step = j+1
                [true1, true2, true3, true4, true5, t6] = truevd(px,pt+ht);
                rhs(1) = rhs(1) + (theta)*coefs1(1,3)*rho*true1 - (theta)*coefs1(1,2)*(ht/2/hx)*true1;
%                 rhs(1) = rhs(1) - coefs(1, 3)*true1/hx^2 + coefs(1, 2)*true1/(2*hx);

                px = gridx(n+1);
                % explicit, time step = j
                [true1, true2, true3, true4, true5, t6] = truevd(px,pt);
                rhs(m) = rhs(m) + (1-theta)*coefs(m,3)*rho*true1 + (1-theta)*coefs(m,2)*(ht/2/hx)*true1; 
                % implicit, time step = j+1
                [true1, true2, true3, true4, true5, t6] = truevd(px,pt+ht);
                rhs(m) = rhs(m) + (theta)*coefs1(m,3)*rho*true1 + (theta)*coefs1(m,2)*(ht/2/hx)*true1;
%                 rhs(m) = rhs(m) - coefs(m, 3)*true1/hx^2 - coefs(m, 2)*true1/(2*hx);

        end
end