% Handles very general BCs given in the form of
% u(a,t) = ga(t)
% u(b,t) = gb(t)
% u(x,t0) = g0(x)

function [u_val,opt_pde] = DirichletBC(x,y,t,opt_pde,opt_var)
    % opt.pde.Uno = -1;
    % opt.pde.BCno = 1;
    % opt.pde.PDEno = 101;
    % % IC/BC choices, Euro 0:2, Amer 10:19, Heston 30:31
    % % will be reassigned in DirichletBC.m
    % opt.pde.Amer = true;
    % opt.pde.RbName = 'Heston American Put';

    mx = length(x);
    my = length(y);
    u_val = zeros(mx*my,1);
    % OptionType = opt_pde.Amer;
    opt_pde.Amer = true;
    x0 = x; y0 = y;
    q1 = opt_var.q1; q2 = opt_var.q2;
    K = opt_var.K; Rf = opt_var.Rf;
    x = x0 * exp(-q1*t);
    y = y0 * exp(-q2*t);
    opt_pde.BCno = 0;

    if (ismember(opt_pde.PDEno, [100,101]))
        switch opt_pde.Rbno
        case {31}
            opt_pde.RbName = 'Heston American Put';
            opt_pde.BCno = 10;
            for i = 1:mx
                u_val((1:my) + (i-1)*(my)) = ...
                    max(K-x(i)*exp(-q1*t),0);
            end
        case {30}
            opt_pde.RbName = 'Heston American Call';
            opt_pde.BCno = 10;
            for i = 1:mx
                u_val((1:my) + (i-1)*(my)) = ...
                    max(x(i)-K*exp(-Rf*t),0);
            end
        case {29}
            opt_pde.RbName = 'Heston European Put';
            opt_pde.BCno = 10;
            opt_pde.Amer = false;
            for i = 1:mx
                u_val((1:my) + (i-1)*(my)) = ...
                    max(K*exp(-Rf*t)-x(i),0);
            end
        case {19}
            opt_pde.RbName = 'American Arithmetic Average Put';
            opt_pde.BCno = 1;
            for i = 1:mx
                u_val((1:my) + (i-1)*(my)) = max(K-(x(i)+y)/2,0);
            end
        case {18}
            opt_pde.RbName = 'American Geometric Average Put';
            opt_pde.BCno = 2;
            for i = 1:mx
                u_val((1:my) + (i-1)*(my)) = ...
                    max(K-real(sqrt(x(i).*y)),0);
            end
        case {17}
            opt_pde.RbName = 'American Spread Put';
            opt_pde.BCno = 3;
            for i = 1:mx
                u_val((1:my) + (i-1)*(my)) = max(-x(i)+y+K,0);
            end
        case {16}
            opt_pde.RbName = 'American Spread Call';
            opt_pde.BCno = 3;
            for i = 1:mx
                u_val((1:my) + (i-1)*(my)) = ...
                    max(x(i)-y-K*exp(-Rf*t),0);
            end
        case {15}
            opt_pde.RbName = 'American Min Put';
            opt_pde.BCno = 1;
            % American Min Put
            for i = 1:mx
                u_val((1:my) + (i-1)*(my)) = max(min(K-x(i),K-y),0);
            end
        case {14}
            % American Max Put
            opt_pde.RbName = 'American Max Put';
            opt_pde.BCno = 2;
            for i = 1:mx
                u_val((1:my) + (i-1)*(my)) = max(max(K-x(i),K-y),0);
            end
        case {13}
            % American Min Call
            opt_pde.RbName = 'American Min Call';
            opt_pde.BCno = 2;
            for i = 1:mx
                u_val((1:my) + (i-1)*(my)) = ...
                    max(min(x(i)-K*exp(-Rf*t),y-K*exp(-Rf*t)),0);
            end
        case {12}
            % American Max Call
            opt_pde.RbName = 'American Max Call';
            opt_pde.BCno = 1;
            for i = 1:mx
                u_val((1:my) + (i-1)*(my)) = ...
                    max(max(x(i)-K*exp(-Rf*t),y-K*exp(-Rf*t)),0);
            end
        case {11}
            % American Max Payoff
            opt_pde.RbName = 'American Max Payoff';
            for i = 1:mx
                u_val((1:my) + (i-1)*(my)) = max(x(i),y);
            end
        case {10}
            opt_pde.BCno = 4; % incorrect BC!
            % American Margrabe
            opt_pde.RbName = 'American Margrabe';
            % BC is NOT Dirichlet at (xmax,ymax)!
            for i = 1:mx
                u_val((1:my) + (i-1)*(my)) = max(x(i)-y,0);
            end
        otherwise
            % All European Rainbow 
            opt_pde.Amer = false;
            for i = 1:mx
                [u_val((1:my) + (i-1)*(my)),opt] = EuroRb(x0(i),y0,t,opt);
            end
        end
        
    else
        opt_pde.Amer = false;
        for i = 1:mx
                u_val((1:my) + (i-1)*(my)) = ...
                    truevd2(x0(i), y0, t, opt_pde)';
        end
    end

end