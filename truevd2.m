% function [t1,t2,t3,t4,t5,t6,t7] = truevd2(x, y, t)
% t1-7 are u,ux,uxx,uy,uyy,uxy,ut respectively
function [t1,t2,t3,t4,t5,t6,t7] = truevd2(x, y, t)

global Uno Uname;
% global etaA etaB etaC nu mu R eee;

% note only take/return 1d vectors
s = max(size(x),max(size(y),size(t)));
o = ones(s); 
z = zeros(s);

% UnoList = [-1:8,10:18,30:31,40:43,50];

switch Uno

case {50}
    Uname = 'exp(x .* y .* t^2)';
    t1 = exp(x .* y .* t.^2);
    t2 = y .* t.^2 .* exp(x .* y .* t.^2);
    t3 = y.^2 .* t.^4 .* exp(x .* y .* t.^2);
    t4 = x .* t.^2 .* exp(x .* y .* t.^2);
    t5 = x.^2 .* t.^4 .* exp(x .* y .* t.^2);
    t6 = t.^2 .* exp(x .* y .* t.^2) + ...
        x .* y .* t.^4 .* exp(x .* y .* t.^2);
    t7 = 2 * x .* y .* t;

case {43}
    Uname = 'x .* y .* t^3';
    t1 = x .* y .* t.^3;
    t2 = y .* t.^3;
    t3 = z;
    t4 = x .* t.^3;
    t5 = z;
    t6 = t.^3;
    t7 = 3 * x .* y .* t.^2;

case {42}
    Uname = 'x .* y .* t^2';
    t1 = x .* y .* t.^2;
    t2 = y .* t.^2;
    t3 = z;
    t4 = x .* t.^2;
    t5 = z;
    t6 = t.^2;
    t7 = 2 * x .* y .* t;

case {41}
    Uname = 'x + y + t^3';
    t1 = x + y + t.^3;
    t2 = o;
    t3 = z;
    t4 = o;
    t5 = z;
    t6 = z;
    t7 = 3*t.^2;

case {40}
    Uname = 'x + y + t^2';
    t1 = x + y + t.^2;
    t2 = o;
    t3 = z;
    t4 = o;
    t5 = z;
    t6 = z;
    t7 = 2*t;

case {31}
    Uname = 'sin(x.*y*c) + t';
    c = 4*pi;
    t1 = sin(x.*y*c) + t;
    t2 = c*y.*cos(x.*y*c);
    t3 = -c^2*y.^2.*sin(x.*y*c);
    t4 = c*x.*cos(x.*y*c);
    t5 = -c^2*x.^2.*sin(x.*y*c);
    t6 = -c^2*x.*y.*sin(x.*y*c) + c*cos(x.*y*c);
    t7 = o;

case {30}
    Uname = 'exp(x.*y) + t';
    t1 = exp(x.*y) + t;
    t2 = y.*exp(x.*y);
    t3 = y.^2.*exp(x.*y);
    t4 = x.*exp(x.*y);
    t5 = x.^2.*exp(x.*y);
    t6 = x.*y.*exp(x.*y) + exp(x.*y);
    t7 = o;

case {18}
    Uname = 'x.^3 .* y.^4 + t';
    t1 = x.^3 .* y.^4 + t;
    t2 = 3*x.^2.*y.^4;
    t3 = 6*x.*y.^4;
    t4 = 4*x.^3.*y.^3;
    t5 = 12*x.^3.*y.^2;
    t6 = 12*x.^2.*y.^3;
    t7 = o;

case {17}
    Uname = 'x.^4 .* y.^3 + t';
    t1 = x.^4 .* y.^3 + t;
    t2 = 4*x.^3.*y.^3;
    t3 = 12*x.^2.*y.^3;
    t4 = 3*x.^4.*y.^2;
    t5 = 6*x.^4.*y;
    t6 = 12*x.^3.*y.^2;
    t7 = o;

case {16}
    Uname = 'x.^3 .* y.^3 + t';
    t1 = x.^3 .* y.^3 + t;
    t2 = 3*x.^2.*y.^3;
    t3 = 6*x.*y.^3;
    t4 = 3*x.^3.*y.^2;
    t5 = 6*x.^3.*y;
    t6 = 9*x.^2.*y.^2;
    t7 = o;

case {15}
    Uname = 'x.^2 .* y.^3 + t';
    t1 = x.^2 .* y.^3 + t;
    t2 = 2*x.*y.^3;
    t3 = 2*y.^3;
    t4 = 3*x.^2.*y.^2;
    t5 = 6*x.^2.*y;
    t6 = 6*x.*y.^2;
    t7 = o;

case {14}
    Uname = 'x.^3 .* y.^2 + t';
    t1 = x.^3 .* y.^2 + t;
    t2 = 3*x.^2.*y.^2;
    t3 = 6*x.*y.^2;
    t4 = 2*x.^3.*y;
    t5 = 2*x.^3;
    t6 = 6*x.^2.*y;
    t7 = o;

case {13}
    Uname = 'x.^2 .* y.^2 + t';
    t1 = x.^2 .* y.^2 + t;
    t2 = 2*x.*y.^2;
    t3 = 2*y.^2;
    t4 = 2*x.^2.*y;
    t5 = 2*x.^2;
    t6 = 4*x.*y;
    t7 = o;

case {12}
    Uname = 'x .* y.^2 + t';
    t1 = x .* y.^2 + t;
    t2 = y.^2;
    t3 = z;
    t4 = 2*x.*y;
    t5 = 2*x;
    t6 = 2*y;
    t7 = o;

case {11}
    Uname = 'x.^2 .* y + t';
    t1 = x.^2 .* y + t;
    t2 = 2*x.*y;
    t3 = 2*y;
    t4 = x.^2;
    t5 = z;
    t6 = 2*x;
    t7 = o;
    
case {10}
    Uname = 'x .* y + t';
    t1 = x .* y + t;
    t2 = y;
    t3 = z;
    t4 = x;
    t5 = z;
    t6 = o;
    t7 = o;

case {8}
    Uname = 'x.^3 + y.^4 + t';
    t1 = x.^3 + y.^4 + t;
    t2 = 3*x.^2;
    t3 = 6*x;
    t4 = 4*y.^3;
    t5 = 12*y.^2;
    t6 = z;
    t7 = o;

case {7}
    Uname = 'x.^4 + y.^3 + t';
    t1 = x.^4 + y.^3 + t;
    t2 = 4*x.^3;
    t3 = 12*x.^2;
    t4 = 3*y.^2;
    t5 = 6*y;
    t6 = z;
    t7 = o;

case {6}
    Uname = 'x.^3 + y.^3 + t';
    t1 = x.^3 + y.^3 + t;
    t2 = 3*x.^2;
    t3 = 6*x;
    t4 = 3*y.^2;
    t5 = 6*y;
    t6 = z;
    t7 = o;

case {5}
    Uname = 'x.^2 + y.^3 + t';
    t1 = x.^2 + y.^3 + t;
    t2 = 2*x;
    t3 = 2*o;
    t4 = 3*y.^2;
    t5 = 6*y;
    t6 = z;
    t7 = o;

case {4}
    Uname = 'x.^3 + y.^2 + t';
    t1 = x.^3 + y.^2 + t;
    t2 = 3*x.^2;
    t3 = 6*x;
    t4 = 2*y;
    t5 = 2*o;
    t6 = z;
    t7 = o;
    
case {3}
    Uname = 'x.^2 + y.^2 + t';
    t1 = x.^2 + y.^2 + t;
    t2 = 2*x;
    t3 = 2*o;
    t4 = 2*y;
    t5 = 2*o;
    t6 = z;
    t7 = o;
    
case {2}
    Uname = 'x.^2 + y + t';
    t1 = x.^2 + y + t;
    t2 = 2*x;
    t3 = 2*o;
    t4 = o;
    t5 = z;
    t6 = z;
    t7 = o;
    
case {1}
    Uname = 'x + y + t';
    t1 = x + y + t;
    t2 = o;
    t3 = z;
    t4 = o;
    t5 = z;
    t6 = z;
    t7 = o;
    
case {0}
    Uname ='1';
    t1 = o;
    t2 = z;
    t3 = z;
    t4 = z;
    t5 = z;
    t6 = z;
    t7 = z;

case {-1}
    Uname ='0';
    t1 = z;
    t2 = z;
    t3 = z;
    t4 = z;
    t5 = z;
    t6 = z;
    t7 = z;

end % switch

end % function








