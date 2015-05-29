% Handles very general BCs given in the form of
% u(a,t) = ga(t)
% u(b,t) = gb(t)
% u(x,t0) = g0(x)

function [u_val] = DirechletBC(x,y,t)
    global PDEno Rbno Rf K T;

    mx = length(x);
    my = length(y);
    u_val = zeros(mx*my,1);

    % only use true vals here
    if (PDEno == 100)
        % Margrabe
        for i = 1:mx
            u_val((1:my) + (i-1)*(my)) = ...
                EuroRb(x(i),y,t);
                % max(x(i)-y,0);
        end
    else
        for i = 1:mx
                u_val((1:my) + (i-1)*(my)) = truevd2(x(i),y,t)';
        end
    end

end