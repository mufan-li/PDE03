% Handles very general BCs given in the form of
% u(a,t) = ga(t)
% u(b,t) = gb(t)
% u(x,t0) = g0(x)

function [u_val] = DirechletBC(x,y)
    global BCno;
    switch BCno
        case {1}
            % base case
            u_val = truevd(x,y);
        otherwise
            % assumes (x,t) are only at IC and BC
            u_val = truevd(x,y);
    end
end