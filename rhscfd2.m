% function [rhs, coefs] = rhscfd2(n, gridx)
function [rhs, coefs, b] = rhscfd2(nx, ny, gridx, gridy, tj)

% Evaluating Boundary Conditions
global BCno Rbno Unif Dim UxStep Uno;
global Rf q1;

mx = nx+1; my = ny+1;
numeq = mx * my;
m = numeq;
hxn = gridx(nx+1) - gridx(nx);

rhs = zeros(m,1);
coefs = zeros(m,7); % including y and cross derivatives

% find coefs - note iterate through y-indices first
% note - coefs are already in the order the diagonals
% specific x,y size
for i = 1:mx
    ind = (1:my)+(i-1)*my;
    [rhs(ind), coefs(ind,1), coefs(ind,2), coefs(ind,3), ...
        coefs(ind,4), coefs(ind,5), coefs(ind,6), coefs(ind,7)] =...
        pde2(gridx(i), gridy(:),tj);
end

% need to find rhs
vfx = zeros(mx,1); vfx(1) = 1;
vlx = zeros(mx,1); vlx(mx) = 1;
vfy = zeros(my,1); vfy(1) = 1;
vly = zeros(my,1); vly(my) = 1;
ox = ones(mx,1);
oy = ones(my,1);

% vectors including corners, size (nx+1) and (ny+1)
% converted to column vectors
ux0 = DirechletBC(gridx,gridy(1),tj);
uxn = DirechletBC(gridx,gridy(ny+1),tj);
u0y = DirechletBC(gridx(1),gridy,tj);
uny = DirechletBC(gridx(nx+1),gridy,tj);

% add the entire rhs
% removed and handled otherwise
switch BCno
    case {10}
        % Heston Model
        % PDE at y=0, D at x=0, Neumann at the rest
        b = kron(vfx,u0y);
        switch Rbno
        case {30} % call phix = 1;
            rhs = rhs - ...
                spdiag(coefs(:,2),m) * ...
                kron(vlx,oy) - ... 
                2/hxn * spdiag(coefs(:,3),m) * ...
                kron(vlx,ht0(oy)+vly);
        otherwise % put and otherwise phix = 0;
            % do nothing
        end
    case {4}
        % PDE BC, for Margrabe
        % x=0,x=max,y=max, and (max,max)
        b = kron(vfx,u0y) + kron(vlx,h(uny));
    case {3}
        % PDE BC, Dirich at x=0 and y=ymax
        % for American Spread Call (x-y)
        b = kron(vfx,ht0(u0y)) + kron(vlx,ht(uny)) + ...
            kron(ht0(uxn),vly) + kron(vfx,ht(u0y));
    case {2}
        % PDE BC, only dirichlet at corners and min 2 sides
        % for American Min Call
        b = kron(vfx,ht0(u0y)) + kron(vlx,ht(uny)) + ...
            kron(ht0(ux0),vfy) + kron(vfx,ht(u0y));
    case {1}
        % PDE BC, only dirichlet at corners and max 2 sides
        % for American Max Call
        b = kron(vfx,ht(u0y)) + kron(vlx,ht0(uny)) + ...
            kron(vlx,ht(uny)) + kron(ht0(uxn),vly);
    otherwise
        b = kron(vfx,ht0(u0y)) + kron(vlx,ht0(uny)) + ...
            kron(ux0,vfy) + kron(uxn,vly);
end

end % end function

% replace head and tail with zero
function b0 = ht0(vec)
    n = length(vec);
    b0 = vec;
    b0(1) = 0;
    b0(n) = 0;
end

% keep only the head and tail
function b0 = ht(vec)
    n = length(vec);
    b0 = vec;
    b0(2:n-1) = 0;
end

% keep only the head
function b0 = h(vec)
    n = length(vec);
    b0 = vec;
    b0(2:n) = 0;
end

% keep only the tail
function b0 = t(vec)
    n = length(vec);
    b0 = vec;
    b0(1:n-1) = 0;
end
