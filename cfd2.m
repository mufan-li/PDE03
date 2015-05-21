% function [A, A2, A1, A0, A1b, A1f] = cfd2(n, gridx, coefs);

function [A] = cfd2(nx, ny, gridx, gridy, coefs)

% Evaluating Boundary Conditions
global BC;

% coefs = [coefu, coefux, coefuxx, coefuy, coefuyy, coefuxy]
m = (nx-1) * (ny-1);
Ix = speye(nx-1);
Iy = speye(ny-1);

% for non-uniform grids
% converted to column vectors
hx = gridx(2:nx+1) - gridx(1:nx);
hx0 = hx(1:nx-1)'; hx1 = hx(2:nx)';
hy = gridy(2:ny+1) - gridy(1:ny);
hy0 = hy(1:ny-1)'; hy1 = hy(2:ny)';

% adjusted matrices for non-uniform grids
% A2x = sptrid(1, -2, 1, nx-1);
A2x = spdiags([t(2./hx0./(hx0+hx1)), ...
        -2./hx0./hx1, ...
        h(2./hx1./(hx0+hx1))], [-1 0 1],nx-1,nx-1);
% A1x = sptrid(-1, 0, 1, nx-1);
A1x = spdiags([t(-hx1./hx0./(hx0+hx1)), ...
        (hx1-hx0)./hx0./hx1, ...
        h(hx0./hx1./(hx0+hx1))],[-1 0 1],nx-1,nx-1);

    % A10 = spdiag((h0.*h1.*(h0+h1)).^-1)*...
    % spdiags(reshape([-h1(2:m).^2, 0, 
    %     -h0.^2+h1.^2, 
    %     0, h0(1:m-1).^2],m,3), [-1,0,1], m, m);

% A2y = sptrid(1, -2, 1, ny-1);
A2y = spdiags([t(2./hy0./(hy0+hy1)), ...
        -2./hy0./hy1, ...
        h(2./hy1./(hy0+hy1))], [-1 0 1],ny-1,ny-1);
% A1y = sptrid(-1, 0, 1, ny-1);
A1y = spdiags([t(-hy1./hy0./(hy0+hy1)), ...
        (hy1-hy0)./hy0./hy1, ...
        h(hy0./hy1./(hy0+hy1))], [-1 0 1],ny-1,ny-1);

A = spdiag(coefs(:,1),m) * kron(Ix,Iy) + ... %u
    spdiag(coefs(:,2),m) * kron(A1x,Iy) + ... %ux
    spdiag(coefs(:,3),m) * kron(A2x,Iy) + ... %uxx
    spdiag(coefs(:,4),m) * kron(Ix,A1y) + ... %uy
    spdiag(coefs(:,5),m) * kron(Ix,A2y) + ... %uyy
    spdiag(coefs(:,6),m) * kron(A1x,A1y) ; %uxy

end

% removes the first element and append zero
function tail = t(vec)
    m = length(vec);
    tail = vec;
    tail(1:m-1) = vec(2:m);
    tail(m) = 0;
end

% removes the last element and append to zero
function head = h(vec)
   m = length(vec);
   head = vec;
   head(2:m) = vec(1:m-1);
   head(1) = 0; 
end

