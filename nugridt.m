% variable ht step
function [ht, gridt, dnorm] = nugridt(ht, stepj, pt, T, ucomp,...
                                gridt, dnorm, Rt)
    global Unift Smin Smax;
    switch Unift
        case {1}
            % simple non-uniform
            if (stepj==1)
%                 ht = ht/10;
                ht = ht/1e3;
            elseif (stepj>2)
                ht1 = dnorm * max([Rt DirechletBC(Smin,pt) ...
                    DirechletBC(Smax,pt) DirechletBC(Smin,pt-ht) ...
                    DirechletBC(Smax,pt-ht)]) * ht / ...
                    max(abs(ucomp(:,stepj-1) - ucomp(:,stepj-2)));
%                 dnorm = dnorm / 2^(ht1>ht*2);
                ht = min(ht1,ht*1.5); % restricting growth
%                 ht = ht1;
            end
            
            % check T
            if (pt+ht>T)
                ht = T - pt;
            end
            
            gridt(stepj+1) = gridt(stepj) + ht;
            
        otherwise
            % uniform ht
            ht = ht; % for explicitness
            gridt(stepj+1) = gridt(stepj) + ht;
    end
end