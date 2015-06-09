% generates a summary of all the fields
% [nx,ny,nt, nit,uv,dxv, dyv,gxv,gyv, fbv]
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
				'Free_Boundary','Change','Ratio'};
			m.cvals = { '','','','',...
				'\t','\t','',...
				'\t','\t','',...
				'\t','\t','',...
				'\t','\t','',...
				'\t','\t','',...
				'\t','\t',''};
			m.pvals = { '%i','%i','%i','%i',...
				'%6.6f','%6.6f','%2.2f',... % price
				'%6.6f','%6.6f','%2.2f',... % delta_x
				'%6.6f','%6.6f','%2.2f',... % delta_y
				'%6.6f','%6.6f','%2.2f',... % gamma_x
				'%6.6f','%6.6f','%2.2f',... % gamma_y
				'%6.6f','%6.6f','%2.2f'};   % fb
			m.value = zeros(ntimes,length(m.colnames));
		end

		% updates after each grid
		function update(m,uj1,Am,Nm,Gm)

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

			% change/ratios calculated at the end
			m.value(Nm.ni,:) = [Nm.nx,Nm.ny,Nm.nt,Nm.nit,...
				uv,0,0,dxv,0,0,dyv,0,0,gxv,0,0,gyv,0,0,...
				fbv,0,0];
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
			global PDEname Uname RbName PenaltyName;
			disp(strcat([PDEname,', u = ',Uname]));
			disp(strcat([RbName,', ',PenaltyName]));

			% get change/ratio
			jn = 1:3:(length(m.colnames) - 4);
			for j = jn+4
				m.value(:,j+1:j+2) = chg( m.value(:,j) );
			end

			print_cols(m,1:7);
			print_cols(m,8:13);
			print_cols(m,14:19);
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

% change and ratio
% vec is a column vector
function [vec3] = chg(vec)
	n = length(vec);
	vec1 = vec(2:n) - vec(1:n-1);
	vec2 = vec1(1:n-2) ./ vec1(2:n-1);
	vec2(isnan(vec2)) = 0;
	vec3 = [[0;vec1] , [0;0;vec2]];
end


