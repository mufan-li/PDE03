function [uj1] = t_step(uj0, rhs, Aim, Aex, gridx, gridy)
    global Penalty;
    m = length(uj0); lx = length(gridx); ly = length(gridy);
    f = DirechletBC(gridx(2:lx-1),gridy(2:ly-1),0);

    switch Penalty;
    case {2,'Penalty Iteration'}
        uj1 = Aim\(Aex*uj0 - rhs);

    case {1,'Operator Splitting'}
        uj1 = Aim\(Aex*uj0 - rhs);

    case {0,'Explicit Penalty'}
    	uj1 = Aim\(Aex*uj0 - rhs);
    	ind = uj1 < f;
    	uj1(ind) = f(ind);

    otherwise
        uj1 = Aim\(Aex*uj0 - rhs);
    end
end
