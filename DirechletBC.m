% Handles very general BCs given in the form of
% u(a,t) = ga(t)
% u(b,t) = gb(t)
% u(x,t0) = g0(x)

function [u_val] = DirechletBC(x,y,t)
    global PDEno Rbno RbName Rf K T;

    mx = length(x);
    my = length(y);
    u_val = zeros(mx*my,1);

    if (PDEno == 100)
        switch Rbno
        case {15}
            % American Min Put
            for i = 1:mx
                u_val((1:my) + (i-1)*(my)) = max(min(K-x(i),K-y),0);
            end
        case {14}
            % American Max Put
            RbName = 'American Max Put';
            for i = 1:mx
                u_val((1:my) + (i-1)*(my)) = max(max(K-x(i),K-y),0);
            end
        case {13}
            % American Min Call
            RbName = 'American Min Call';
            for i = 1:mx
                u_val((1:my) + (i-1)*(my)) = max(min(x(i)-K,y-K),0);
            end
        case {12}
            % American Max Call
            RbName = 'American Max Call';
            for i = 1:mx
                u_val((1:my) + (i-1)*(my)) = max(max(x(i)-K,y-K),0);
            end
        case {11}
            % American Max Payoff
            RbName = 'American Max Payoff';
            for i = 1:mx
                u_val((1:my) + (i-1)*(my)) = max(x(i),y);
            end
        case {10}
            % American Margrabe
            RbName = 'American Margrabe';
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