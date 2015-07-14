% generates a summary of all the fields
% handle class
classdef summary < handle
	properties
		colnames
		cvals
		pvals
		value
	end

	methods
		% initialize
		function m = summary(ntimes)

			m.colnames = {'nx','ny','nt','nit',...
				'Price','Change','Ratio',...
				'Delta_x','Change','Ratio',...
				'Delta_y','Change','Ratio',...
				'Gamma_x','Change','Ratio',...
				'Gamma_y','Change','Ratio',...
				'Free_Boundary','Change','Ratio',...
				'Cond','Ratio',...
				'NormInv','Ratio'};

			% if additional tab needed
			m.cvals = { '','','','',...
				'\t','\t','',...
				'\t','\t','',...
				'\t','\t','',...
				'\t','\t','',...
				'\t','\t','',...
				'\t','\t','',...
				'','','',''};

			m.pvals = { '%i','%i','%i','%i',...
				'%6.6f','%6.6f','%2.2f',... % price
				'%6.6f','%6.6f','%2.2f',... % delta_x
				'%6.6f','%6.6f','%2.2f',... % delta_y
				'%6.6f','%6.6f','%2.2f',... % gamma_x
				'%6.6f','%6.6f','%2.2f',... % gamma_y
				'%6.6f','%6.6f','%2.2f',... % fb
				'%6.0f','%2.2f','%6.0f','%2.2f',... % cond,norm
				};
			m.value = zeros(ntimes,length(m.colnames));
		end

		% updates after each grid
		function update(m,uj1,Am,Nm,Gm,Aim)

			% determine value
			uv = intp(uj1,Gm);

			Adx = kron(Am.A1x,Am.Iy) + Am.Ab;
			dxv = intp(Adx*uj1,Gm);

			Ady = kron(Am.Ix,Am.A1y) + Am.Ab;
			dyv = intp(Ady*uj1,Gm);

			Agx = kron(Am.A2x,Am.Iy) + Am.Ab;
			gxv = intp(Agx*uj1,Gm);

			Agy = kron(Am.Ix,Am.A2y) + Am.Ab;
			gyv = intp(Agy*uj1,Gm);

			fbv = 0;

			if (Nm.ni<=3)
				condA = cond(full(Aim),Inf);
				normA = norm(full(Aim)^-1,Inf);
			else
				condA = 0;
				normA = 0;
			end

			% change/ratios calculated at the end
			m.value(Nm.ni,:) = [Nm.nx,Nm.ny,Nm.nt,Nm.nit,...
				uv,0,0,dxv,0,0,dyv,0,0,gxv,0,0,gyv,0,0,...
				fbv,0,0,condA,0,normA,0];
		end

		function print_cols(m,cols)
			% print column names
			for j = cols
				fprintf(['%s',cell2mat(m.cvals(j)),'\t'], ...
					cell2mat(m.colnames(j)));
			end
			fprintf('\n');

			% print values
			for i = 1:size(m.value,1)
				for j = cols
					fprintf([cell2mat(m.pvals(j)),'\t'], ...
						m.value(i,j));
				end
				fprintf('\n');
			end

			fprintf('\n');

		end

		% prints results
		function print(m)
			global PDEname Uname RbName PenaltyName Unift;
			global xp yp T Gridno Unift;
			disp([PDEname,', u = ',Uname, ', ',...
				 RbName,', ',PenaltyName]);
			disp(['Gridno: ',int2str(Gridno),...
				', Unift: ',int2str(Unift), ', xp: ', ...
				num2str(xp), ', yp: ', num2str(yp)]);

			% get change/ratio
			% omit first 4 and last 4 columns
			for j0 = 1:(length(m.colnames) - 8)/3
				j = (j0-1)*3 + 5;
				m.value(:,j+1:j+2) = chg( m.value(:,j), j0);
			end
			% get ratio of matrix cond and norm
			for j = length(m.colnames)-2:2:length(m.colnames)
				m.value(:,j) = ratio(m.value(:,j-1));
			end

			print_cols(m,1:7); % grid and price
			% print_cols(m,8:13); % deltas
			% print_cols(m,14:19); % gammas
			% print_cols(m,20:22); % fb
			print_cols(m,23:26); % cond,norm
		end

		% plot the final surface
		function plot(m,uj1,Gm)
			figure;
			mesh(Gm.gx,Gm.gy,...
				reshape(uj1,length(Gm.gy),length(Gm.gx)));
		end

		% plot the interior surface
		function plot_int(m,uj1,Gm)
			nx = length(Gm.gx);
			ny = length(Gm.gy);
			v = reshape(uj1,length(Gm.gy),length(Gm.gx));
			figure;
			mesh(Gm.gx(2:nx-1),Gm.gy(2:ny-1),v(2:ny-1,2:nx-1));
		end

		% plot the surface of greeks
		% only plot the interior points
		function plot_greeks(m,uj1,Gm,Am)
			Adx = kron(Am.A1x,Am.Iy) + Am.Ab;
			plot_int(m,Adx*uj1,Gm);

			Ady = kron(Am.Ix,Am.A1y) + Am.Ab;
			plot_int(m,Ady*uj1,Gm);

			Agx = kron(Am.A2x,Am.Iy) + Am.Ab;
			plot_int(m,Agx*uj1,Gm);

			Agy = kron(Am.Ix,Am.A2y) + Am.Ab;
			plot_int(m,Agy*uj1,Gm);
		end

	end % end methods
end % end class

% interpolation at points given grid
function [up] = intp(uj1,Gm)
	gx = Gm.gx;
	gy = Gm.gy;
	x = Gm.x;
	y = Gm.y;

	mx = length(gx); my = length(gy);
	uj1m = reshape(uj1,my,mx);
	ix = find(gx>x,1); iy = find(gy>y,1);

	um = uj1m(iy-1:iy,ix-1:ix);
	hx = gx(ix)-gx(ix-1); hy = gy(iy)-gy(iy-1);
	dx = 1-abs(gx(ix-1:ix)-x)/hx; dy = 1-abs(gy(iy-1:iy)-y)/hy;
	up = dy * um * dx';
end

% return change and ratio as 2 columns
% vec is a column vector
function [vec3] = chg(vec,j0)
	global xp yp T OptionType;

	n = length(vec);

	% only return trueval for price for now
	if (j0==1) & (OptionType==0)
		vec1 = abs(vec - EuroRb(xp,yp,T));
	else
		vec1 = [0;vec(2:n) - vec(1:n-1)];
	end

	vec2 = vec1(1:n-1) ./ vec1(2:n);
	vec2(isnan(vec2) | isinf(vec2)) = 0;
	vec3 = [vec1 , [0;vec2]];
end

% return the ratio as one vector
function [vec1] = ratio(vec)
	n = length(vec);
	vec1 = [0; vec(2:n) ./ vec(1:n-1)];
end
