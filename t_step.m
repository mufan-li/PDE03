function [uj1,aux1,nit,P] = t_step(uj0, rhs, Aim, Aex, gridx, gridy,...
					 aux0, htj, tj1, nit)

    global Penalty PenaltyName OptionType;
    m = length(uj0); mx = length(gridx); my = length(gridy);
    f = DirechletBC(gridx,gridy,0);
    aux1 = aux0; % unassigned
    nit = nit+1; % only penalty is different

    P = 0*spdiag(f); % init

    switch OptionType
    case {1}
	    switch Penalty;
	    case {3,'Discrete Penalty'}
	    	PenaltyName = 'Discrete Penalty';
	    	tol = 1e-6;
	    	ujk0 = uj0;
	    	P = spdiag((uj0-f)<0)/tol;

	    	for k = 1:10
	        	ujk1 = (Aim+P)\(Aex*uj0 - rhs + P*f);
	        	P0 = P;
	        	P = spdiag((ujk1-f)<0)/tol;

	        	if (max( ... % component-wise division
        			abs(ujk1-ujk0) ./ max(ones(size(ujk1)),ujk1) ...
	        	) < tol ) ...
	        		|| isequal(P0,P)
	        		break
	        	end
	        	ujk0 = ujk1;
	        end
	        uj1 = ujk1;
	        nit = nit+k-1; % additional iterations

	    case {4,'Splitting HaHo13'}
	    	PenaltyName = 'Splitting HaHo13';
	        vj1 = Aim\(Aex*uj0 - rhs + htj*aux0);
	        aux1 = max( aux0 - 1/htj*(vj1 - f) ,0);
	        uj1 = max( vj1-htj*aux0 ,f);

	    case {2,'Operator Splitting'}
	    	PenaltyName = 'Operator Splitting';
	        vj1 = Aim\(Aex*uj0 - rhs + htj*aux0);
	        aux1 = max( aux0 - 1/htj*(vj1 - f) ,0);
	        uj1 = max( vj1+htj*(aux1-aux0) ,f);

	    case {1,'Explicit Penalty'}
	    	PenaltyName = 'Explicit Penalty';
	    	uj1 = Aim\(Aex*uj0 - rhs);
	    	ind = uj1 < f;
	    	uj1(ind) = f(ind);

	    otherwise
	    	PenaltyName = 'None';
	        uj1 = Aim\(Aex*uj0 - rhs);
	    end
	    
	otherwise
		PenaltyName = 'None';
        uj1 = Aim\(Aex*uj0 - rhs);
	end

end

