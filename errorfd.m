% function errg = errorfd(ngrid, gridx, n, uvct, nn, errg, truef);
function errg = errorfd(ngridx, ngridy, gridx, gridy, ni, uj1, errg);

trueval = zeros(size(uj1));

for i = 2:ngridx-1
	trueval((1:ngridy-2) + (i-2)*(ngridy-2)) = ...
		truevd(gridx(i),gridy(2:ngridy-1));
end

errg(1,ni) = max(trueval - uj1);

end