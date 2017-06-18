% function [A,Ad,Ab,An,Am] = cfd2(coefs, pgrid, opt_pde, opt_var)

function [A,Ad,Ab,An,Am] = cfd2(coefs, pgrid, opt_pde, opt_var)
    % cfd2(nx, ny, gridx, gridy, coefs, bet)

% Evaluating Boundary Conditions
% global BCno;

nx = pgrid.nx;
ny = pgrid.ny;
gridx = pgrid.gridx;
gridy = pgrid.gridy;
BCno = opt_pde.BCno;
bet = opt_var.beta;

m = (nx+1) * (ny+1);
Ixc = speye(nx+1); Ix = Ixc; Ix(1,1) = 0; Ix(nx+1,nx+1) = 0;
Iyc = speye(ny+1); Iy = Iyc; Iy(1,1) = 0; Iy(ny+1,ny+1) = 0;
Ic = speye(m); In = 0*Ic; % Neumann

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
A2x = A2fb(A2x,hx,nx);

A1x = spdiags([t0(-hx1./hx0./(hx0+hx1)), ...
        b0((hx1-hx0)./hx0./hx1), ...
        h0(hx0./hx1./(hx0+hx1))],[-1 0 1],nx+1,nx+1);
A1x = A1fb(A1x,hx,nx);

A2y = spdiags([t0(2./hy0./(hy0+hy1)), ...
        b0(-2./hy0./hy1), ...
        h0(2./hy1./(hy0+hy1))], [-1 0 1],ny+1,ny+1);
A2y = A2fb(A2y,hy,ny);

A1y = spdiags([t0(-hy1./hy0./(hy0+hy1)), ...
        b0((hy1-hy0)./hy0./hy1), ...
        h0(hy0./hy1./(hy0+hy1))], [-1 0 1],ny+1,ny+1);
% one-sided convection
A1y = A1yconv(A1y,gridy,hy,bet,opt_pde); % Heston
A1y = A1fb(A1y,hy,ny);

% Neumann matrices
An = 0 * Ic;

% general form for the PDE BC on sides
% - with discretizations of orthogonal directions
switch BCno
    case {10} % Heston
        Im = kron(Ix,Iy0);
        In = kron(Ix,Jyn) + kron(Jxn,Iyc); % Neumann
        T2yn = 0*Iy; T2yn(ny+1,ny:ny+1) = [2 -2] / hy(ny)^2;
        T2xn = 0*Ix; T2xn(nx+1,nx:nx+1) = [2 -2] / hx(nx)^2;
        T1y0 = Jy0 * A1y;

        Anu = spdiag(coefs(:,1),m) * In; %u
        Anx = spdiag(coefs(:,2),m) * kron(Ix*A1x,Jyn); %ux
        Anxx = spdiag(coefs(:,3),m) * ...
            (kron(Ix*A2x,Jyn) + kron(T2xn,Iyn) ) ; %uxx
        Any = spdiag(coefs(:,4),m) * ...
            (kron(Jxn,Iy*A1y) + kron(Jxn,T1y0) ) ; %uy
        Anyy = spdiag(coefs(:,5),m) * ...
            (kron(Ixn,T2yn) + kron(Jxn,Iy*A2y) ) ; %uyy

        An = Anu + Anx + Anxx + Any + Anyy;

    case {4} % Margrabe
        Im = kron(Ixn,Iyn);

    case {3} % Spread Call
        Im = kron(Ixn,Iyc) - kron(Jxn,Jy0);

    case {2} % PDE BC for American Min Call
        Im = kron(Ixn,Iyn) - kron(Jxn,Jyn);

    case {1} % PDE BC for American Max Call
        Im = kron(Ix0,Ix0) - kron(Jx0,Jx0);

    otherwise
        Im = kron(Ix,Iy);
end

Ab = Ic - Im - In; % Dirichlet BC

Ad0 = spdiag(coefs(:,1),m) * Im; %u
Adx = spdiag(coefs(:,2),m) * kron(A1x,Iyc); %ux
Adxx = spdiag(coefs(:,3),m) * kron(A2x,Iyc); %uxx
Ady = spdiag(coefs(:,4),m) * kron(Ixc,A1y); %uy
Adyy = spdiag(coefs(:,5),m) * kron(Ixc,A2y); %uyy
Adxy = spdiag(coefs(:,6),m) * kron(A1x,A1y); %uxy
% only apply PDE at interior points
Ad = Im * (Ad0+Adx+Adxx+Ady+Adyy+Adxy); 

A = Ab + Ad + An; % BC conditions

Am.Ix = Ix;
Am.Iy = Iy;
Am.Im = Im;
Am.A1x = A1x;
Am.A1y = A1y;
Am.A2x = A2x;
Am.A2y = A2y;
Am.Ab = Ab;

% Am.Ad0 = Ad0;
% Am.Adx = Adx;
% Am.Adxx = Adxx;
% Am.Ady = Ady;
% Am.Adyy = Adyy;
% Am.Adxy = Adxy;

Am.An = An;
Am.Anu = Anu;
Am.Anx = Anx;
Am.Anxx = Anxx;
Am.Any = Any;
Am.Anyy = Anyy;

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

% one-sided discretizations
function A2 = A2fb(A2,h,n)
    % second order
    a = -2*(2*h(1)+2*h(2)+h(3))/h(1)/h(2)/(h(1)+h(2));
    b = 2*(2*h(1)+h(2)+h(3))/h(2)/h(3)/(h(1)+h(2));
    c = -2*(2*h(1)+h(2))/h(3)/(h(2)+h(3))/(h(1)+h(2)+h(3));
    d = -(a+b+c);
    A2(1,1:4) = [d a b c];

    a = -2*(2*h(n)+2*h(n-1)+h(n-2))/h(n)/h(n-1)/(h(n)+h(n-1));
    b = 2*(2*h(n)+h(n-1)+h(n-2))/h(n-1)/h(n-2)/(h(n)+h(n-1));
    c = -2*(2*h(n)+h(n-1))/...
        h(n-2)/(h(n-1)+h(n-2))/(h(n)+h(n-1)+h(n-2));
    d = -(a+b+c);
    A2(n+1,n-2:n+1) = [c b a d];
end

% one-sided discretizations
function A1 = A1fb(A1,h,n)
    %% second order
    % A1(1,1:3) = [-(2*h(1) + h(2))/h(1)/(h(1)+h(2)), ...
    %     (h(1)+h(2))/h(1)/h(2), ...
    %     -h(1)/h(2)/(h(1)+h(2)) ];

    % A1(n+1,n-1:n+1) = [h(n)/h(n-1)/(h(n)+h(n-1)), ...
    %     -(h(n) + h(n-1))/h(n)/h(n-1), ...
    %     (2*h(n) + h(n-1))/h(n)/(h(n) + h(n-1)) ];

    %% first order
    A1(1,1:2) = [-1/h(1), 1/h(1)];

    A1(n+1,n:n+1) = [-1/h(1), 1/h(1)];
end

% one-sided convection for Heston
function A1yc = A1yconv(A1y,gridy,hy,bet,opt_pde)
    % global PDEno;
    PDEno = opt_pde.PDEno;
    switch PDEno
    case {101}
        my = length(gridy);
        nbet = sum(gridy<bet);
        A1yf = spdiag([1./hy,0]) * sptrid(0, -1, 1, my);
        A1yb = spdiag([0,1./hy]) * sptrid(-1, 1, 0, my);
        
        A1yc = A1yb;
        A1yc(1:nbet,:) = A1yf(1:nbet,:);

        % A1y = spdiags([t0(-hy1./hy0./(hy0+hy1)), ...
        % b0((hy1-hy0)./hy0./hy1), ...
        % h0(hy0./hy1./(hy0+hy1))], [-1 0 1],ny+1,ny+1);
    otherwise
        A1yc = A1y;
end
end


