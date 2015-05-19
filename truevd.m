% function [true1, true2, true3, true4, true5, t6, t7, t8, t9] = truevd(x)
%
% returns the values of the solution function, 1st, 2nd, 3rd and 4th
% derivatives on x
% Watch! Some functions are incomplete! (missing 3rd and 4th derivatives)
% Also some derivatives 5-8 are present.

function [true1] = truevd(x, y)

global Uno Uname;
% global etaA etaB etaC nu mu R eee;

% note only take/return 1d vectors
o = ones(max(size(x),size(y))); z = zeros(max(size(x),size(y)));

% UnoVec = [-1:7,10:13,30:31];

switch Uno
    case {31}
        Uname = 'sin(x.*y/1e-1)';
        true1 = sin(x.*y/1e-1);
    
    case {30}
        Uname = 'exp(x.*y)';
        true1 = exp(x.*y);
    
    case {13}
        Uname = 'x.^2 .* y.^2';
        true1 = x.^2 .* y.^2;

    case {12}
        Uname = 'x .* y.^2';
        true1 = x .* y.^2;

    case {11}
        Uname = 'x.^2 .* y';
        true1 = x.^2 .* y;
        
    case {10}
        Uname = 'x .* y';
        true1 = x .* y;
        
    case {7}
        Uname = 'x.^3 + y.^3';
        true1 = x.^3 + y.^3;
    
    case {6}
        Uname = 'x.^2 + y.^4';
        true1 = x.^2 + y.^4;
    
    case {5}
        Uname = 'x.^2 + y.^3';
        true1 = x.^2 + y.^3;
    
    case {4}
        Uname = 'x.^3 + y.^2';
        true1 = x.^3 + y.^2;
        
    case {3}
        Uname = 'x.^2 + y.^2';
        true1 = x.^2 + y.^2;
        
    case {2}
        Uname = 'x.^2 + y';
        true1 = x.^2 + y;
        
    case {1}
        Uname = 'x + y';
        true1 = x+y;
        
    case {0}
        Uname ='1';
        true1 = o;

    case {-1}
        Uname ='0';
        true1 = z;
end % switch

end % function








