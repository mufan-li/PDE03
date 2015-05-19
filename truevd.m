% function [t1,t2,t3,t4,t5,t6] = truevd(x, y)
% t1-6 are u,ux,uxx,uy,uyy,uxy respectively
function [t1,t2,t3,t4,t5,t6] = truevd(x, y)

global Uno Uname;
% global etaA etaB etaC nu mu R eee;

% note only take/return 1d vectors
o = ones(max(size(x),size(y))); z = zeros(max(size(x),size(y)));

% UnoVec = [-1:7,10:13,30:31];

switch Uno
    case {31}
        Uname = 'sin(x.*y/1e-1)';
        t1 = sin(x.*y/1e-1);
    
    case {30}
        Uname = 'exp(x.*y)';
        t1 = exp(x.*y);
    
    case {13}
        Uname = 'x.^2 .* y.^2';
        t1 = x.^2 .* y.^2;

    case {12}
        Uname = 'x .* y.^2';
        t1 = x .* y.^2;

    case {11}
        Uname = 'x.^2 .* y';
        t1 = x.^2 .* y;
        
    case {10}
        Uname = 'x .* y';
        t1 = x .* y;
        
    case {7}
        Uname = 'x.^3 + y.^3';
        t1 = x.^3 + y.^3;
    
    case {6}
        Uname = 'x.^2 + y.^4';
        t1 = x.^2 + y.^4;
    
    case {5}
        Uname = 'x.^2 + y.^3';
        t1 = x.^2 + y.^3;
    
    case {4}
        Uname = 'x.^3 + y.^2';
        t1 = x.^3 + y.^2;
        
    case {3}
        Uname = 'x.^2 + y.^2';
        t1 = x.^2 + y.^2;
        
    case {2}
        Uname = 'x.^2 + y';
        t1 = x.^2 + y;
        t2 = 2*x;
        t3 = 2*o;
        t4 = o;
        t5 = z;
        t6 = z;
        
    case {1}
        Uname = 'x + y';
        t1 = x+y;
        t2 = o;
        t3 = z;
        t4 = o;
        t5 = z;
        t6 = z;
        
    case {0}
        Uname ='1';
        t1 = o;
        t2 = z;
        t3 = z;
        t4 = z;
        t5 = z;
        t6 = z;

    case {-1}
        Uname ='0';
        t1 = z;
        t2 = z;
        t3 = z;
        t4 = z;
        t5 = z;
        t6 = z;

end % switch

end % function








