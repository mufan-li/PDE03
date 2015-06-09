% generates a summary of all the fields
% [nx,ny,nt, nit,uv,dxv, dyv,gxv,gyv, fbv]
% handle class
classdef summary < handle
	properties
		colnames
		value
	end

	methods
		% initialize
		function m = summary(ntimes)
			m.colnames = {'nx','ny','nt','nit','Price',...
				'Delta_x','Delta_y','Gamma_x','Gamma_y',...
				'Free_Boundary'};
			m.value = zeros(ntimes,length(m.colnames));
		end

		% updates after each grid
		function update(m,uj1,Am,Nm,Gm)

			% determine values
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

			m.value(Nm.ni,:) = [Nm.nx,Nm.ny,Nm.nt,Nm.nit,...
				uv,dxv,dyv,gxv,gyv,fbv];
		end

		% prints results
		function print(obj)
			disp(obj.colnames);
			disp(obj.value);
		end
	end
end

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