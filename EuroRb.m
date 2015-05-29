% function [t1] = EuroRb(x, y, t);
% - returns value of the European Rainbown options
function [t1] = EuroRb(x, y, t);
global Rbno T Sx Sy rho Rf K Smin Smax;

switch Rbno
	case {2}
		% Max Payoff
		s = sqrt(Sx^2 + Sy^2 - 2*rho * Sx * Sy);

		d12p = ( log(x./y) + 1/2*s^2.*(t) ) ./ (s .* sqrt(t));
		d12p(isnan(d12p)) = 0; % NaN when 0/0 - value is zero
		d12m = d12p - s .* sqrt(t);

		d21p = ( log(y./x) + 1/2*s^2.*(t) ) ./ (s .* sqrt(t));
		d21p(isnan(d21p)) = 0; % NaN when 0/0 - value is zero
		d21m = d21p - s .* sqrt(t);

		d1p = ( log(x./K) + Rf + 1/2*Sx^2.*t ) ./ (Sx .* sqrt(t));
		d1m = d1p - Sx .* sqrt(t);
		d2p = ( log(y./K) + Rf + 1/2*Sy^2.*t ) ./ (Sy .* sqrt(t));
		d2m = d2p - Sy .* sqrt(t);

		r12 = (Sy - rho*Sx)/s; %r13,2
		r21 = (Sx - rho*Sy)/s; %r23,1
		
		t1 = x .* N2cdf(-d21m,d1p,r21) + y .* N2cdf(-d12m,d2p,r12) ...
			- K * exp(-Rf*t) * (1 - N2cdf(-d1m,-d2m,rho));
	case {1}
		% Max Payoff
		s = sqrt(Sx^2 + Sy^2 - 2*rho * Sx * Sy);
		d12p = ( log(x./y) + 1/2*s^2.*(t) ) ./ (s .* sqrt(t));
		d12p(isnan(d12p)) = 0; % NaN when 0/0 - value is zero
		d12m = d12p - s .* sqrt(t);

		d21p = ( log(y./x) + 1/2*s^2.*(t) ) ./ (s .* sqrt(t));
		d21p(isnan(d21p)) = 0; % NaN when 0/0 - value is zero
		d21m = d21p - s .* sqrt(t);

		t1 = x .* normcdf(-d21m) + y .* normcdf(-d12m);

	case {0}
		% Margrabe
		s = sqrt(Sx^2 + Sy^2 - 2*rho * Sx * Sy);
		d1 = ( log(x./y) + 1/2*s^2.*(t) ) ./ (s .* sqrt(t));
		d1(isnan(d1)) = 0; % NaN when 0/0 - value is zero
		d2 = d1 - s .* sqrt(t);

		t1 = x .* normcdf(d1) - y .* normcdf(d2);

	otherwise
		% case 0
		% note when Uno = -1, this returns zero, which is desired
		t1 = truevd2(x, y, t);
end



% consider adding delta, gamma etc.

end

% 2D mvncdf
function p = N2cdf(x1,x2,rho)
	x0 = zeros(max(size(x1),size(x2)));
	x = [(x0+x1)', (x0+x2)'];
	mu = [0 0];
	SG = [1 rho;rho 1];
	p = mvncdf(x,mu,SG);
	p = p';
end
