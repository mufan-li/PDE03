% function [A, A2, A1, A0, A1b, A1f] = cfd2(n, gridx, coefs);
% returns the matrix of second-order centered FD discretization
% of second-order DEs ()*u'' + ()*u' + ()*u = g, Dirichlet BCs
% if coefs does not exist then the DE is u'' = g
% optionally return A2, A1, A0 matrices corresponding to u'', u', u
% optionally return A1b, A1f matrices corresponding to u' backward, forward

function [A, B, A2, A1, A0, B2, B1, B0] = cfd2(n, gridx, coefs, coefs1, pt, ht, theta);

% Evaluating Boundary Conditions
global BC Unif UxStep;
switch Unif
    case{1}
        % non-uniform grid
        numeq = n-1;	% assume Dirichlet conditions
        m = numeq;	% simplicity has charm
        hx = gridx(2) - gridx(1);	% assume uniform grid
        rho = ht/hx^2;
        % assume non-uniform grid stay the same for all t
        h = gridx(2:n+1) - gridx(1:n);
        h0 = h(1:m);
        h1 = h(2:m+1);

        % centered Ux Step % originally intended for forward/backward Ux Step
        A10 = -1*spdiag(coefs1(:,2),m)*theta*ht*spdiag((h0.*h1.*(h0+h1)).^-1)*...
            spdiags(reshape([-h1(2:m).^2, 0, -h0.^2+h1.^2, 0, h0(1:m-1).^2],m,3), [-1,0,1], m, m);
        B10 = spdiag(coefs(:,2),m)*(1-theta)*ht*spdiag((h0.*h1.*(h0+h1)).^-1)*...
            spdiags(reshape([-h1(2:m).^2, 0, -h0.^2+h1.^2, 0, h0(1:m-1).^2],m,3), [-1,0,1], m, m);
        % implicit matrix
        A2 = -1*theta*spdiag(coefs(:, 3))*ht*spdiag((h0.*h1.*(h0+h1)).^-1)*...
            spdiags(reshape([2*h1(2:m), 0, -2*(h0+h1), 0, 2*h0(1:m-1)],m,3), [-1,0,1], m, m)+spdiag(1,m);
        A1 = A10;
        % note that r*u(i,j+1) is subtracted because it is not moved to the other side of the eqn
        A0 = -1*spdiag(coefs1(:,1),m)*theta*ht*spdiag(1,m);
        A = A2+A1+A0;

        % explicit matrix
        B2 = (1-theta)*spdiag(coefs(:, 3))*ht*spdiag((h0.*h1.*(h0+h1)).^-1)*...
            spdiags(reshape([2*h1(2:m), 0, -2*(h0+h1), 0, 2*h0(1:m-1)],m,3), [-1,0,1], m, m)+spdiag(1,m);
        B1 = B10;
        B0 = spdiag(coefs(:,1),m)*(1-theta)*ht*spdiag(1,m);
        B = B2+B1+B0;

    otherwise
        switch BC

        case{4}
            numeq = n-1;	% assume Dirichlet conditions
            m = numeq;	% simplicity has charm
            hx = gridx(2) - gridx(1);	% assume uniform grid
%             A = sptrid(1, -2, 1, m)/hx^2;
            A = spdiags([-1*ones(m, 1), 16*ones(m, 1), -30*ones(m, 1), ...
                16*ones(m, 1), -1*ones(m, 1)], [-2, -1, 0, 1, 2], m, m)/hx^2/12;
            A(1,1:5) = [12, -2, -17, 8, -1] /65/hx^2;
            A(m,m-4:m) = [1, -8, 17, 2, -12] /65/hx^2;
            
            if (nargin > 2)			% if coefs is present
                E = speye(m);
                A1 = sptrid(-1, 0, 1, m)/(2*hx);
                A = spdiag(coefs(:, 3))*A + spdiag(coefs(:, 2))*A1 + spdiag(coefs(:, 1))*E;
            end
            if (nargout > 1)
                A0 = speye(m);
                A1 = sptrid(-1, 0, 1, m)/(2*hx);
                A2 = sptrid(1, -2, 1, m)/hx^2;
                A1b= sptrid(-1, 1, 0, m)/hx;
                A1f= sptrid( 0,-1, 1, m)/hx;
            end
            
        case{3}
            % One sided first order
            % assume General Robin conditions

            numeq = n+1;	% assume General Robin conditions
            m = numeq;	% simplicity has charm
            hx = gridx(2) - gridx(1);	% assume uniform grid

            % creates n*n sparse matrix
            A = sptrid(1, -2, 1, m);

            % add Robin BC
            [coef_a,coef_b,coef_f] = robinBC(gridx(1));
            A(1,1) = coef_a-coef_b/hx;
            A(1,2) = coef_b/hx;
            [coef_a,coef_b,coef_f] = robinBC(gridx(m));
            A(m,m) = coef_a+coef_b/hx;
            A(m,m-1) = -coef_b/hx;
            A=A/hx^2;

            if (nargin > 2)			% if coefs is present
                E = speye(m);
                A1 = sptrid(-1, 0, 1, m)/(2*hx);
                % need to adjust for BC
                [coef_a,coef_b,coef_f] = robinBC(gridx(1));
                A1(1,1)=-coef_a/coef_b;
                A1(1,2)=0;
                [coef_a,coef_b,coef_f] = robinBC(gridx(m));
                A1(m,m)=-coef_a/coef_b;
                A1(m,m-1)=0;
                A = spdiag(coefs(:, 3))*A + spdiag(coefs(:, 2))*A1 + spdiag(coefs(:, 1))*E;
            end
            if (nargout > 1)
                A0 = speye(m);
                A1 = sptrid(-1, 0, 1, m)/(2*hx);
                A2 = sptrid(1, -2, 1, m)/hx^2;
                A1b= sptrid(-1, 1, 0, m)/hx;
                A1f= sptrid( 0,-1, 1, m)/hx;
            end


        case{2}
            % One sided second order
            % assume General Robin conditions

            numeq = n+1;	% assume General Robin conditions
            m = numeq;	% simplicity has charm
            hx = gridx(2) - gridx(1);	% assume uniform grid

            % creates n*n sparse matrix
            A = sptrid(1, -2, 1, m);

            % add Robin BC
            [coef_a,coef_b,coef_f] = robinBC(gridx(1));
            A(1,1) = coef_a-3*coef_b/(2*hx);
            A(1,2) = 4*coef_b/(2*hx);
            A(1,3) = -coef_b/(2*hx);
            [coef_a,coef_b,coef_f] = robinBC(gridx(m));
            A(m,m) = coef_a+3*coef_b/(2*hx);
            A(m,m-1) = -4*coef_b/(2*hx);
            A(m,m-2) = coef_b/(2*hx);
            A=A/hx^2;

            if (nargin > 2)			% if coefs is present
                E = speye(m);
                A1 = sptrid(-1, 0, 1, m)/(2*hx);
                % need to adjust for BC
                [coef_a,coef_b,coef_f] = robinBC(gridx(1));
                A1(1,1)=-coef_a/coef_b;
                A1(1,2)=0;
                [coef_a,coef_b,coef_f] = robinBC(gridx(m));
                A1(m,m)=-coef_a/coef_b;
                A1(m,m-1)=0;
                A = spdiag(coefs(:, 3))*A + spdiag(coefs(:, 2))*A1 + spdiag(coefs(:, 1))*E;
            end
            if (nargout > 1)
                A0 = speye(m);
                A1 = sptrid(-1, 0, 1, m)/(2*hx);
                A2 = sptrid(1, -2, 1, m)/hx^2;
                A1b= sptrid(-1, 1, 0, m)/hx;
                A1f= sptrid( 0,-1, 1, m)/hx;
            end

        case{1}
            % centered differences with ficticious point
            % assume General Robin conditions

            numeq = n+1;	% assume General Robin conditions
            m = numeq;	% simplicity has charm
            hx = gridx(2) - gridx(1);	% assume uniform grid

            % creates n*n sparse matrix
            A = sptrid(1, -2, 1, m);

            % add Robin BC
            [coef_a,coef_b,coef_f] = robinBC(gridx(1));
            A(1,1) = 2*hx*coef_a/coef_b - 2;
            A(1,2) = 2;
            [coef_a,coef_b,coef_f] = robinBC(gridx(m));
            A(m,m) = - 2*hx*coef_a/coef_b - 2;
            A(m,m-1) = 2;
            A=A/hx^2;

            if (nargin > 2)			% if coefs is present
                E = speye(m);
                A1 = sptrid(-1, 0, 1, m)/(2*hx);
                % need to adjust for BC
                [coef_a,coef_b,coef_f] = robinBC(gridx(1));
                A1(1,1)=-coef_a/coef_b;
                A1(1,2)=0;
                [coef_a,coef_b,coef_f] = robinBC(gridx(m));
                A1(m,m)=-coef_a/coef_b;
                A1(m,m-1)=0;
                A = spdiag(coefs(:, 3))*A + spdiag(coefs(:, 2))*A1 + spdiag(coefs(:, 1))*E;
            end
            if (nargout > 1)
                A0 = speye(m);
                A1 = sptrid(-1, 0, 1, m)/(2*hx);
                A2 = sptrid(1, -2, 1, m)/hx^2;
                A1b= sptrid(-1, 1, 0, m)/hx;
                A1f= sptrid( 0,-1, 1, m)/hx;
            end

        otherwise

            numeq = n-1;	% assume Dirichlet conditions
            m = numeq;	% simplicity has charm
            hx = gridx(2) - gridx(1);	% assume uniform grid
            rho = ht/hx^2;
            
            % setup Ux matrix
            switch UxStep
                case {'forward'}
                    A10 = -1*spdiag(coefs1(:,2),m)*theta*(ht/hx)*sptrid(0,-1,1,m);
                    B10 = spdiag(coefs(:,2),m)*(1-theta)*(ht/hx)*sptrid(0,-1,1,m);
                case {'backward'}
                    A10 = -1*spdiag(coefs1(:,2),m)*theta*(ht/hx)*sptrid(-1,1,0,m);
                    B10 = spdiag(coefs(:,2),m)*(1-theta)*(ht/hx)*sptrid(-1,1,0,m);
                case {'centered'}
                    A10 = -1*spdiag(coefs1(:,2),m)*theta*(ht/2/hx)*sptrid(-1,0,1,m);
                    B10 = spdiag(coefs(:,2),m)*(1-theta)*(ht/2/hx)*sptrid(-1,0,1,m);
                otherwise
                    % centered default
                    disp('default to center stepping');
                    A10 = -1*spdiag(coefs1(:,2),m)*theta*(ht/2/hx)*sptrid(-1,0,1,m);
                    B10 = spdiag(coefs(:,2),m)*(1-theta)*(ht/2/hx)*sptrid(-1,0,1,m);
            end
            
            % implicit matrix
            A2 = spdiag(coefs1(:,3),m)*theta*rho*sptrid(-1, 2, -1, m) + spdiag(1,m);
            A1 = A10;
            % note that r*u(i,j+1) is subtracted because it is not moved to the other side of the eqn
            A0 = -1*spdiag(coefs1(:,1),m)*theta*ht*spdiag(1,m);
            A = A2+A1+A0;
            
            % explicit matrix
            B2 = spdiag(coefs(:,3),m)*(1-theta)*rho*sptrid(1,-2,1,m)+spdiag(1,m);
            B1 = B10;
            B0 = spdiag(coefs(:,1),m)*(1-theta)*ht*spdiag(1,m);
            B = B2+B1+B0;
            
%             if (nargin > 2)			% if coefs is present
%                 E = speye(m);
%                 A1 = sptrid(-1, 0, 1, m);
%                 A = spdiag(coefs(:, 3))*A + spdiag(coefs(:, 2))*A1 + spdiag(coefs(:, 1))*E;
%             end
%             if (nargout > 1)
%                 A0 = speye(m);
%                 A1 = sptrid(-1, 0, 1, m)/(2*hx);
%                 A2 = sptrid(1, -2, 1, m)/hx^2;
%                 A1b= sptrid(-1, 1, 0, m)/hx;
%                 A1f= sptrid( 0,-1, 1, m)/hx;
%             end
        end
end