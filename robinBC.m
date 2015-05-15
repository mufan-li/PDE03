% Handles General Robin BC Functions
% a(x) * u + b(x) * ux = f(x)

function [coef_a,coef_b,coef_f] = robinBC(x)
    
    % define a and b
    coef_a=1;
    coef_b=1;
    [true1,true2] = truevd(x);
    coef_f=coef_a*true1 + coef_b*true2;
end