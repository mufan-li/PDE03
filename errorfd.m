% function errg = errorfd(ngrid, gridx, n, uvct, nn, errg, truef);
function errg = errorfd(ngrid, gridx, pt, n, uvct, nn, errg, truef);

global BC;
switch BC

case{0}
% No BC
if (nargin < 8) truef = 'truevd'; end;
for i = 1:n-1
    px = gridx(i+1);
    [true(1) true(2) true(3)] = feval(truef, px, pt);
    errg(1, nn) = max(errg(1, nn), abs(true(1)-uvct(i)));
end

otherwise
% assume General Robin conditions

if (nargin < 8) truef = 'truevd'; end;
for i = 1:n
    px = gridx(i);
    [true(1) true(2) true(3)] = feval(truef, px, pt);
    errg(1, nn) = max(errg(1, nn), abs(true(1)-uvct(i)));
end


end