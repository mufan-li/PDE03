% Handles very general BCs given in the form of
% u(a,t) = ga(t)
% u(b,t) = gb(t)
% u(x,t0) = g0(x)

function [u_val] = DirechletBC(x,y,t)
    global PDEno Rbno Rf K T;

    mx = length(x);
    my = length(y);
    u_val = zeros(mx*my,1);

    if (PDEno == 100)
        switch Rbno
        case {10}
            % American Margrabe
            for i = 1:mx
                u_val((1:my) + (i-1)*(my)) = max(x(i)-y,0);
            end
        otherwise
            % All European Rainbow 
            for i = 1:mx
                u_val((1:my) + (i-1)*(my)) = EuroRb(x(i),y,t);
            end
        end
        
    else
        for i = 1:mx
                u_val((1:my) + (i-1)*(my)) = truevd2(x(i),y,t)';
        end
    end

end