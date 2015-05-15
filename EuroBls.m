% Black Scholes Formula

function [Call, Put, DeltaCall, DeltaPut, Gamma] = EuroBls(S, K, Rf, T, t, SigmaC)
    if (T==t)
        d1 = inf*((S>=K)-0.5);
    else
        d1 = (log(S./K) + (Rf + 0.5*SigmaC.^2).*(T-t)) ./ (SigmaC.*(T-t).^0.5);
    end
    d2 = d1 - SigmaC.*(T-t).^0.5;
    
    Call = S .* normcdf(d1) - K .* exp(-Rf.*(T-t)) .* normcdf(d2);
    Put = -S .* normcdf(-d1) + K .* exp(-Rf.*(T-t)) .* normcdf(-d2);
    
    DeltaCall = normcdf(d1);
    DeltaPut = normcdf(d1)-1;
    Gamma = normpdf(d1)./(S.*SigmaC.*(T-t).^0.5);
end
