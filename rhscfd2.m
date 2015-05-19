% function [rhs, coefs] = rhscfd2(n, gridx)
function [rhs, coefs] = rhscfd2(nx, ny, gridx, gridy)

% Evaluating Boundary Conditions
global BC BCno Unif Dim UxStep Uno;

mx = nx-1; my = ny-1;
numeq = mx * my;
m = numeq;
hx = gridx(2) - gridx(1);
hy = gridy(2) - gridy(1);
rhs = zeros(m,1);
coefs = zeros(m,6); % including y and cross derivatives

% find coefs - note iterate through y-indices first
% note - coefs are already in the order the diagonals
% specific x,y size
for i = 1:mx
    ind = (1:my)+(i-1)*my;
    [rhs(ind), coefs(ind,1), coefs(ind,2), coefs(ind,3), ...
        coefs(ind,4), coefs(ind,5), coefs(ind,6)] = ...
        pde1(gridx(i+1), gridy(2:ny));
end

% need to find rhs
vfx = zeros(mx,1); vfx(1) = 1;
vlx = zeros(mx,1); vlx(mx) = 1;
vfy = zeros(my,1); vfy(1) = 1;
vly = zeros(my,1); vly(my) = 1;

% vectors including corners, size (nx+1) and (ny+1)
% converted to column vectors - consider fixing in DirechletBC()
ux0 = DirechletBC(gridx,gridy(1))';
uxn = DirechletBC(gridx,gridy(ny+1))';
u0y = DirechletBC(gridx(1),gridy)';
uny = DirechletBC(gridx(nx+1),gridy)';

% corner units in the y direction size (ny-1)
% used to remove the extra corners in cross derivatives
uc0y = [u0y(1); zeros(ny-3,1); u0y(ny+1)];
ucny = [uny(1); zeros(ny-3,1); uny(ny+1)];

% need to add the other coefficients
% note - ux0 and uxn are terms used for uyy
rhs_add2x = - (kron(vfx,u0y(2:ny)) + kron(vlx,uny(2:ny))) ./ hx^2;
rhs_add2y = - (kron(ux0(2:nx),vfy) + kron(uxn(2:nx),vly)) ./ hy^2;

rhs_addxy = - (kron(vfx,u0y(1:ny-1)-u0y(3:ny+1)-uc0y) + ...
    kron(vlx,-uny(1:ny-1)+uny(3:ny+1)-ucny) + ...
    kron(ux0(1:nx-1)-ux0(3:nx+1),vfy) + ...
    kron(-uxn(1:nx-1)-ux0(3:nx+1),vly)) ./ (4*hx*hy);

rhs_add1x = - (kron(vfx,-u0y(2:ny)) + kron(uny(2:ny),vlx)) ./ (2*hx);
rhs_add1y = - (kron(-ux0(2:nx),vfy) + kron(uxn(2:nx),vly)) ./ (2*hy);

rhs = rhs + coefs(:,2) .* rhs_add1x + coefs(:,3) .* rhs_add2x ...
    + coefs(:,4) .* rhs_add1y + coefs(:,5) .* rhs_add2y ...
    + coefs(:,6) .* rhs_addxy;

end