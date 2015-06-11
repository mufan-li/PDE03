% Handles very general BCs given in the form of
% u(a,t) = ga(t)
% u(b,t) = gb(t)
% u(x,t0) = g0(x)

function [u_val] = DirechletBC(x,y,t)
    global PDEno Rbno RbName Rf K T q1 q2 OptionType;

    mx = length(x);
    my = length(y);
    u_val = zeros(mx*my,1);
    OptionType = 1;
    x0 = x; y0 = y;
    x = x0 * exp(-q1*t);
    y = y0 * exp(-q2*t);

    if (PDEno == 100)
        switch Rbno
        case {19}
            RbName = 'American Arithmetic Average Put';
            for i = 1:mx
                u_val((1:my) + (i-1)*(my)) = max(K-(x(i)+y)/2,0);
            end
        case {18}
            RbName = 'American Geometric Average Put';
            for i = 1:mx
                u_val((1:my) + (i-1)*(my)) = ...
                    max(K-real(sqrt(x(i).*y)),0);
            end
        case {17}
            RbName = 'American Spread Put';
            for i = 1:mx
                u_val((1:my) + (i-1)*(my)) = max(-x(i)+y+K,0);
            end
        case {16}
            RbName = 'American Spread Call';
            for i = 1:mx
                u_val((1:my) + (i-1)*(my)) = ...
                    max(x(i)-y-K*exp(-Rf*t),0);
            end
        case {15}
            RbName = 'American Min Put';
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
                u_val((1:my) + (i-1)*(my)) = ...
                    max(min(x(i)-K*exp(-Rf*t),y-K*exp(-Rf*t)),0);
            end
        case {12}
            % American Max Call
            RbName = 'American Max Call';
            for i = 1:mx
                u_val((1:my) + (i-1)*(my)) = ...
                    max(max(x(i)-K*exp(-Rf*t),y-K*exp(-Rf*t)),0);
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
            OptionType = 0;
            for i = 1:mx
                u_val((1:my) + (i-1)*(my)) = EuroRb(x0(i),y0,t);
            end
        end
        
    else
        OptionType = 0;
        for i = 1:mx
                u_val((1:my) + (i-1)*(my)) = ...
                    truevd2(x0(i),y0,t)';
        end
    end

end