% function [t1,t2,t3,t4,t5,t6] = truevd(x, y)
% t1-6 are u,ux,uxx,uy,uyy,uxy respectively
function [t1,t2,t3,t4,t5,t6] = truevd(x, y)

global Uno Uname;
% global etaA etaB etaC nu mu R eee;

% note only take/return 1d vectors
o = ones(max(size(x),size(y)));
z = zeros(max(size(x),size(y)));

% UnoList = [-1:8,10:18,30:31];

switch Uno
    case {31}
        Uname = 'sin(x.*y*c)';
        c = 4*pi;
        t1 = sin(x.*y*c);
        t2 = c*y.*cos(x.*y*c);
        t3 = -c^2*y.^2.*sin(x.*y*c);
        t4 = c*x.*cos(x.*y*c);
        t5 = -c^2*x.^2.*sin(x.*y*c);
        t6 = -c^2*x.*y.*sin(x.*y*c) + c*cos(x.*y*c);
    
    case {30}
        Uname = 'exp(x.*y)';
        t1 = exp(x.*y);
        t2 = y.*exp(x.*y);
        t3 = y.^2.*exp(x.*y);
        t4 = x.*exp(x.*y);
        t5 = x.^2.*exp(x.*y);
        t6 = x.*y.*exp(x.*y) + exp(x.*y);

    case {18}
        Uname = 'x.^3 .* y.^4';
        t1 = x.^3 .* y.^4;
        t2 = 3*x.^2.*y.^4;
        t3 = 6*x.*y.^4;
        t4 = 4*x.^3.*y.^3;
        t5 = 12*x.^3.*y.^2;
        t6 = 12*x.^2.*y.^3;

    case {17}
        Uname = 'x.^4 .* y.^3';
        t1 = x.^4 .* y.^3;
        t2 = 4*x.^3.*y.^3;
        t3 = 12*x.^2.*y.^3;
        t4 = 3*x.^4.*y.^2;
        t5 = 6*x.^4.*y;
        t6 = 12*x.^3.*y.^2;

    case {16}
        Uname = 'x.^3 .* y.^3';
        t1 = x.^3 .* y.^3;
        t2 = 3*x.^2.*y.^3;
        t3 = 6*x.*y.^3;
        t4 = 3*x.^3.*y.^2;
        t5 = 6*x.^3.*y;
        t6 = 9*x.^2.*y.^2;

    case {15}
        Uname = 'x.^2 .* y.^3';
        t1 = x.^2 .* y.^3;
        t2 = 2*x.*y.^3;
        t3 = 2*y.^3;
        t4 = 3*x.^2.*y.^2;
        t5 = 6*x.^2.*y;
        t6 = 6*x.*y.^2;

    case {14}
        Uname = 'x.^3 .* y.^2';
        t1 = x.^3 .* y.^2;
        t2 = 3*x.^2.*y.^2;
        t3 = 6*x.*y.^2;
        t4 = 2*x.^3.*y;
        t5 = 2*x.^3;
        t6 = 6*x.^2.*y;

    case {13}
        Uname = 'x.^2 .* y.^2';
        t1 = x.^2 .* y.^2;
        t2 = 2*x.*y.^2;
        t3 = 2*y.^2;
        t4 = 2*x.^2.*y;
        t5 = 2*x.^2;
        t6 = 4*x.*y;

    case {12}
        Uname = 'x .* y.^2';
        t1 = x .* y.^2;
        t2 = y.^2;
        t3 = z;
        t4 = 2*x.*y;
        t5 = 2*x;
        t6 = 2*y;

    case {11}
        Uname = 'x.^2 .* y';
        t1 = x.^2 .* y;
        t2 = 2*x.*y;
        t3 = 2*y;
        t4 = x.^2;
        t5 = z;
        t6 = 2*x;
        
    case {10}
        Uname = 'x .* y';
        t1 = x .* y;
        t2 = y;
        t3 = z;
        t4 = x;
        t5 = z;
        t6 = o;

    case {8}
        Uname = 'x.^3 + y.^4';
        t1 = x.^3 + y.^4;
        t2 = 3*x.^2;
        t3 = 6*x;
        t4 = 4*y.^3;
        t5 = 12*y.^2;
        t6 = z;

    case {7}
        Uname = 'x.^4 + y.^3';
        t1 = x.^4 + y.^3;
        t2 = 4*x.^3;
        t3 = 12*x.^2;
        t4 = 3*y.^2;
        t5 = 6*y;
        t6 = z;

    case {6}
        Uname = 'x.^3 + y.^3';
        t1 = x.^3 + y.^3;
        t2 = 3*x.^2;
        t3 = 6*x;
        t4 = 3*y.^2;
        t5 = 6*y;
        t6 = z;
    
    case {5}
        Uname = 'x.^2 + y.^3';
        t1 = x.^2 + y.^3;
        t2 = 2*x;
        t3 = 2*o;
        t4 = 3*y.^2;
        t5 = 6*y;
        t6 = z;
    
    case {4}
        Uname = 'x.^3 + y.^2';
        t1 = x.^3 + y.^2;
        t2 = 3*x.^2;
        t3 = 6*x;
        t4 = 2*y;
        t5 = 2*o;
        t6 = z;
        
    case {3}
        Uname = 'x.^2 + y.^2';
        t1 = x.^2 + y.^2;
        t2 = 2*x;
        t3 = 2*o;
        t4 = 2*y;
        t5 = 2*o;
        t6 = z;
        
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








