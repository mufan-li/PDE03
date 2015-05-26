% function [t1] = EuroRb(x, y, t);
% - returns value of the European Rainbown options
function [t1] = EuroRb(x, y, t);
global T Sx Sy rho Rf K Smin Smax;

s = sqrt(Sx^2 + Sy^2 - rho * Sx * Sy);
d1 = ( log(x./y) + 1/2*s^2.*(T-t) ) ./ (s .* sqrt(T-t));
d1(isnan(d1)) = 0; % NaN when 0/0 - value is zero
d2 = d1 - s .* sqrt(T-t);

t1 = x .* normcdf(d1) - y .* normcdf(d2);

% consider adding delta, gamma etc.

end