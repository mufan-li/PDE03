% variable ht step
function [ht, gridt] = nugridt(ht, stepj, pt, T, uj1, uj0,...
                                gridt, dnorm, Rt)
    global Unift Smin Smax;
    switch Unift
        case {1}
            % simple non-uniform
            if (stepj==1)
                ht = ht/1e3;
            elseif (stepj>2)
                % element wise
                ht1 = dnorm * ht * ...
                min(max(max(Rt,uj0),uj1) ./ abs(uj0-uj1));
%                 dnorm = dnorm / 2^(ht1>ht*2);
                ht = min(ht1,ht*3); % restricting growth
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