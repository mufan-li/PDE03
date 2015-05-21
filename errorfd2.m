% function errg = errorfd(ngrid, gridx, n, uvct, nn, errg, truef);
function errg = errorfd2(ngridx, ngridy, gridx, gridy, ...
						ni, uj1, errg, bt);

trueval = zeros(size(uj1));

for i = 2:ngridx-1
	trueval((1:ngridy-2) + (i-2)*(ngridy-2)) = ...
		truevd2(gridx(i),gridy(2:ngridy-1),bt);
end

errg(1,ni) = max(trueval - uj1);

end