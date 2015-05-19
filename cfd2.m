% function [A, A2, A1, A0, A1b, A1f] = cfd2(n, gridx, coefs);
% returns the matrix of second-order centered FD discretization
% of second-order DEs ()*u'' + ()*u' + ()*u = g, Dirichlet BCs
% if coefs does not exist then the DE is u'' = g
% optionally return A2, A1, A0 matrices corresponding to u'', u', u
% optionally return A1b, A1f matrices corresponding to u' backward, forward

function [A] = cfd2(nx, ny, gridx, gridy, coefs)

% Evaluating Boundary Conditions
global BC;

% coefs = [coefu, coefux, coefuxx, coefuy, coefuyy, coefuxy]
m = (nx-1) * (ny-1);
Ix = speye(nx-1);
Iy = speye(ny-1);

% need to adjust these matrices for non-uniform grids
A2x = sptrid(1, -2, 1, nx-1);
A1x = sptrid(-1, 0, 1, nx-1);
A2y = sptrid(1, -2, 1, ny-1);
A1y = sptrid(-1, 0, 1, ny-1);

A = spdiag(coefs(:,1),m) * kron(Ix,Iy) + ... %u
    spdiag(coefs(:,2),m) * kron(A1x,Iy) + ... %ux
    spdiag(coefs(:,3),m) * kron(A2x,Iy) + ... %uxx
    spdiag(coefs(:,4),m) * kron(Ix,A1y) + ... %uy
    spdiag(coefs(:,5),m) * kron(Ix,A2y) + ... %uyy
    spdiag(coefs(:,6),m) * kron(A1x,A1y); %uxy

% A = spdiag(coefs(:, 3))*A + spdiag(coefs(:, 2))*A1 + spdiag(coefs(:, 1))*E;

end

%%%%%%%%%%%%%%%%%%%%%%
% % non-uniform grid %
%%%%%%%%%%%%%%%%%%%%%%

% numeq = n-1;  % assume Dirichlet conditions
% m = numeq;    % simplicity has charm
% hx = gridx(2) - gridx(1); % assume uniform grid
% rho = ht/hx^2;
% % assume non-uniform grid stay the same for all t
% h = gridx(2:n+1) - gridx(1:n);
% h0 = h(1:m);
% h1 = h(2:m+1);

% % centered Ux Step % originally intended for forward/backward Ux Step
% A10 = -1*spdiag(coefs1(:,2),m)*theta*ht*spdiag((h0.*h1.*(h0+h1)).^-1)*...
%     spdiags(reshape([-h1(2:m).^2, 0, -h0.^2+h1.^2, 0, h0(1:m-1).^2],m,3), [-1,0,1], m, m);
% B10 = spdiag(coefs(:,2),m)*(1-theta)*ht*spdiag((h0.*h1.*(h0+h1)).^-1)*...
%     spdiags(reshape([-h1(2:m).^2, 0, -h0.^2+h1.^2, 0, h0(1:m-1).^2],m,3), [-1,0,1], m, m);
% % implicit matrix
% A2 = -1*theta*spdiag(coefs(:, 3))*ht*spdiag((h0.*h1.*(h0+h1)).^-1)*...
%     spdiags(reshape([2*h1(2:m), 0, -2*(h0+h1), 0, 2*h0(1:m-1)],m,3), [-1,0,1], m, m)+spdiag(1,m);
% A1 = A10;
% % note that r*u(i,j+1) is subtracted because it is not moved to the other side of the eqn
% A0 = -1*spdiag(coefs1(:,1),m)*theta*ht*spdiag(1,m);
% A = A2+A1+A0;

% % explicit matrix
% B2 = (1-theta)*spdiag(coefs(:, 3))*ht*spdiag((h0.*h1.*(h0+h1)).^-1)*...
%     spdiags(reshape([2*h1(2:m), 0, -2*(h0+h1), 0, 2*h0(1:m-1)],m,3), [-1,0,1], m, m)+spdiag(1,m);
% B1 = B10;
% B0 = spdiag(coefs(:,1),m)*(1-theta)*ht*spdiag(1,m);
% B = B2+B1+B0;