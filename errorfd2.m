% function errg = errorfd(ngrid, gridx, n, uvct, nn, errg, truef);
function [errg,trueval] = errorfd2(ngridx, ngridy, gridx, gridy, ...
						ni, uj1, errg, bt);

global PDEno Rbno K OptionType;

trueval = zeros(size(uj1));

switch OptionType

case {1}
	errg(1,ni) = intp(uj1,gridx,gridy,K,K);

otherwise
	switch PDEno
	case {100}
		% European Rainbow
		truef = 'EuroRb';
	otherwise
		truef = 'truevd2';
	end

	for i = 1:ngridx
		trueval((1:ngridy) + (i-1)*(ngridy)) = ...
			feval(truef,gridx(i),gridy(1:ngridy),bt);
	end

	errg(1,ni) = max(abs(trueval - uj1));
end

end

% assumes x,y fall within uj1
function [up] = intp(uj1,gx,gy,x,y)
	mx = length(gx); my = length(gy);
	uj1m = reshape(uj1,my,mx);
	ix = find(gx>x,1); iy = find(gy>y,1);

	um = uj1m(iy-1:iy,ix-1:ix);
	hx = gx(ix)-gx(ix-1); hy = gy(iy)-gy(iy-1);
	dx = 1-abs(gx(ix-1:ix)-x)/hx; dy = 1-abs(gy(iy-1:iy)-y)/hy;
	up = dy * um * dx';
end