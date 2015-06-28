% function [A, A2, A1, A0, A1b, A1f] = cfd2(n, gridx, coefs);

function [A,Ad,Ab,Am] = cfd2(nx, ny, gridx, gridy, coefs)

% Evaluating Boundary Conditions
global BCno;

m = (nx+1) * (ny+1);
Ix = speye(nx+1); Ix(1,1) = 0; Ix(nx+1,nx+1) = 0;
Iy = speye(ny+1); Iy(1,1) = 0; Iy(ny+1,ny+1) = 0;
Ic = speye(m);

Jx = 0 .* speye(nx+1); Jx(1,1) = 1; Jx(nx+1,nx+1) = 1;
Jx0 = Jx; Jx0(nx+1,nx+1) = 0; Jxn = Jx; Jxn(1,1) = 0;
Jy = 0 .* speye(ny+1); Jy(1,1) = 1; Jy(ny+1,ny+1) = 1;
Jy0 = Jy; Jy0(ny+1,ny+1) = 0; Jyn = Jy; Jyn(1,1) = 0;

Ix0 = Ix+Jx0; Iy0 = Iy+Jy0;
Ixn = Ix+Jxn; Iyn = Iy+Jyn;

% for non-uniform grids
% converted to column vectors
hx = gridx(2:nx+1) - gridx(1:nx);
hx0 = hx(1:nx-1)'; hx1 = hx(2:nx)';
hy = gridy(2:ny+1) - gridy(1:ny);
hy0 = hy(1:ny-1)'; hy1 = hy(2:ny)';

% adjusted matrices for non-uniform grids
A2x = spdiags([t0(2./hx0./(hx0+hx1)), ...
        b0(-2./hx0./hx1), ...
        h0(2./hx1./(hx0+hx1))], [-1 0 1],nx+1,nx+1);

A1x = spdiags([t0(-hx1./hx0./(hx0+hx1)), ...
        b0((hx1-hx0)./hx0./hx1), ...
        h0(hx0./hx1./(hx0+hx1))],[-1 0 1],nx+1,nx+1);

A2y = spdiags([t0(2./hy0./(hy0+hy1)), ...
        b0(-2./hy0./hy1), ...
        h0(2./hy1./(hy0+hy1))], [-1 0 1],ny+1,ny+1);

A1y = spdiags([t0(-hy1./hy0./(hy0+hy1)), ...
        b0((hy1-hy0)./hy0./hy1), ...
        h0(hy0./hy1./(hy0+hy1))], [-1 0 1],ny+1,ny+1);

% temporary solution, need general form
switch BCno
    case {3} % Spread Call
        Im = kron(Ixn,Iy0) - kron(Jxn,Jy0);
        Ixp = Ixn; Iyp = Iy0;

    case {2} % PDE BC for American Min Call
        Im = kron(Ixn,Iyn) - kron(Jxn,Jyn);
        Ixp = Ixn; Iyp = Iyn;

    case {1} % PDE BC for American Max Call
        Im = kron(Ix0,Ix0) - kron(Jx0,Jx0);
        Ixp = Ix0; Iyp = Iy0;

    otherwise
        Im = kron(Ix,Iy);
        Ixp = Ix; Iyp = Iy;
end

Ab = Ic - Im;

Ad = spdiag(coefs(:,1),m) * Im + ... %u
    spdiag(coefs(:,2),m) * kron(A1x,Iyp) + ... %ux
    spdiag(coefs(:,3),m) * kron(A2x,Iyp) + ... %uxx
    spdiag(coefs(:,4),m) * kron(Ixp,A1y) + ... %uy
    spdiag(coefs(:,5),m) * kron(Ixp,A2y) + ... %uyy
    spdiag(coefs(:,6),m) * kron(A1x,A1y) ; %uxy

A = Ab + Ad; % BC conditions

Am.Ix = Ix;
Am.Iy = Iy;
Am.Im = Im;
Am.A1x = A1x;
Am.A1y = A1y;
Am.A2x = A2x;
Am.A2y = A2y;
Am.Ab = Ab;

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

% append 2 zeros to the end
function tail = t0(vec)
    tail = [vec;0;0];
end

% add 2 zeros to the head
function head = h0(vec)
    head = [0;0;vec];
end

% add zeros in head and tail
function b = b0(vec)
    b = [0;vec;0];
end
