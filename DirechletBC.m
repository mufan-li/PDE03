% Handles very general BCs given in the form of
% u(a,t) = ga(t)
% u(b,t) = gb(t)
% u(x,t0) = g0(x)

function [u_val] = DirechletBC(x,y,t)
    global PDEno;

    mx = length(x);
    my = length(y);
    u_val = zeros(mx*my,1);

    switch PDEno
        case {100}
            % base case
            for i = 1:mx
                % convert to column vector
            	u_val((1:my) + (i-1)*(my)) = EuroRb(x(i),y,t)';
            end
            
        otherwise
            % assumes (x,t) are only at IC and BC
            for i = 1:mx
                u_val((1:my) + (i-1)*(my)) = truevd2(x(i),y,t)';
            end
    end
end