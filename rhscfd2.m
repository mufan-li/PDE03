% function [rhs, coefs] = rhscfd2(n, gridx)
function [rhs, coefs] = rhscfd2(nx, ny, gridx, gridy, tj)

% Evaluating Boundary Conditions
global BC BCno Unif Dim UxStep Uno;

mx = nx-1; my = ny-1;
numeq = mx * my;
m = numeq;

% note - need the beginning and end of the grid
hx = gridx(2:nx+1) - gridx(1:nx);
hx0 = hx(1:nx-1)'; hx1 = hx(2:nx)';
hx00 = hx(1); hx0n = hx(nx-1);
hx10 = hx(2); hx1n = hx(nx);

hy = gridy(2:ny+1) - gridy(1:ny);
hy0 = hy(1:ny-1)'; hy1 = hy(2:ny)';
hy00 = hy(1); hy0n = hy(ny-1);
hy10 = hy(2); hy1n = hy(ny);

rhs = zeros(m,1);
coefs = zeros(m,6); % including y and cross derivatives

% find coefs - note iterate through y-indices first
% note - coefs are already in the order the diagonals
% specific x,y size
for i = 1:mx
    ind = (1:my)+(i-1)*my;
    [rhs(ind), coefs(ind,1), coefs(ind,2), coefs(ind,3), ...
        coefs(ind,4), coefs(ind,5), coefs(ind,6)] = ...
        pde2(gridx(i+1), gridy(2:ny),tj);
end

% need to find rhs
vfx = zeros(mx,1); vfx(1) = 1;
vlx = zeros(mx,1); vlx(mx) = 1;
vfy = zeros(my,1); vfy(1) = 1;
vly = zeros(my,1); vly(my) = 1;

% vectors including corners, size (nx+1) and (ny+1)
% converted to column vectors
ux0 = DirechletBC(gridx,gridy(1),tj);
uxn = DirechletBC(gridx,gridy(ny+1),tj);
u0y = DirechletBC(gridx(1),gridy,tj);
uny = DirechletBC(gridx(nx+1),gridy,tj);

%%%%%%%%%%%%%%%%%%%%
% cross derivative %
%%%%%%%%%%%%%%%%%%%%

% first find the diagonal vectors
lx = -hx1./hx0./(hx0+hx1);
dx = (hx1-hx0)./hx0./hx1;
vx = hx0./hx1./(hx0+hx1);

ly = -hy1./hy0./(hy0+hy1);
dy = (hy1-hy0)./hy0./hy1;
vy = hy0./hy1./(hy0+hy1);

% then find the rhs for each border
v0y = lx(1) * (ly.*u0y(1:ny-1) + dy.*u0y(2:ny) + vy.*u0y(3:ny+1)); 
vny = vx(nx-1) * (ly.*uny(1:ny-1) + dy.*uny(2:ny) + vy.*uny(3:ny+1));
vx0 = ly(1) * (lx.*ux0(1:nx-1) + dx.*ux0(2:nx) + vx.*ux0(3:nx+1));
vxn = vy(ny-1) * (lx.*uxn(1:nx-1) + dx.*uxn(2:nx) + vx.*uxn(3:nx+1));

% corner units in the y direction size (ny-1)
% used to remove the extra corners in cross derivatives
uc0y = [u0y(1) *ly(1); zeros(ny-3,1); u0y(ny+1) *vy(ny-1)] *lx(1);
ucny = [uny(1) *ly(1); zeros(ny-3,1); uny(ny+1) *vy(ny-1)] *vx(nx-1);

% note - ux0 and uxn are terms used for uyy
rhs_add2x = - (kron(vfx,u0y(2:ny)*2/hx00/(hx00+hx10)) ...
            + kron(vlx,uny(2:ny)*2/hx1n/(hx0n+hx1n)));

rhs_add2y = - (kron(ux0(2:nx)*2/hy00/(hy00+hy10),vfy) ...
            + kron(uxn(2:nx)*2/hy1n/(hy0n+hy1n),vly));

rhs_add1x = - (kron(vfx,-u0y(2:ny)*hx10/hx00/(hx00+hx10)) ...
            + kron(vlx,uny(2:ny)*hx0n/hx1n/(hx0n+hx1n)));

rhs_add1y = - (kron(-ux0(2:nx)*hy10/hy00/(hy00+hy10),vfy) ...
            + kron(uxn(2:nx)*hy0n/hy1n/(hy0n+hy1n),vly));

rhs_addxy = - ( kron(vfx,v0y - uc0y) + kron(vlx,vny - ucny) + ...
                kron(vx0,vfy) + kron(vxn,vly) );

% add the entire rhs
rhs = rhs + coefs(:,2) .* rhs_add1x + coefs(:,3) .* rhs_add2x ...
    + coefs(:,4) .* rhs_add1y + coefs(:,5) .* rhs_add2y ...
    + coefs(:,6) .* rhs_addxy;

end












