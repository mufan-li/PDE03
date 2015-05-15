% function [true1, true2, true3, true4, true5, t6, t7, t8, t9] = truevd(x)
%
% returns the values of the solution function, 1st, 2nd, 3rd and 4th
% derivatives on x
% Watch! Some functions are incomplete! (missing 3rd and 4th derivatives)
% Also some derivatives 5-8 are present.

function [true1, true2, true3, true4, true5, t6, t7, t8, t9] = truevd(x, t)

global Uno Uname Dim;
global etaA etaB etaC nu mu R eee;

% note only take/return 1d vectors
o = ones(max(size(x),size(t))); z = zeros(max(size(x),size(t)));
% temporarily only !!!
true4 = z;
true5 = z;
t6 = z; t7 = z; t8 = z; t9 = z;

switch Dim
    case {2}
        % Let t6 = ut, t7 = utt
        switch Uno
            case {31}
                Uname = 'sin(x.*t/1e-1)';
                true1 = sin(x.*t/1e-1); true2 = (t/1e-1) .* cos(x.*t/1e-1); 
                true3 = -1 * (t/1e-1).^2 .* sin(x.*t/1e-1);
                t6 = (x/1e-1) .* cos(x.*t/1e-1); 
            
            case {30}
                Uname = 'exp(x.*t)';
                true1 = exp(x.*t); true2 = t .* exp(x.*t); true3 = t.^2 .* exp(x.*t);
                t6 = x .* exp(x.*t);
            
            case {12}
                Uname = 'x .* t.^2';
                true1 = x .* t; true2 = t; true3 = z;
                t6 = 2 * x .* t; t7 = 2 * x;

            case {11}
                Uname = 'x.^2 .* t';
                true1 = x.^2 .* t; true2 = 2 * x .* t; true3 = 2 * t;
                t6 = x.^2; t7 = z;
                
            case {10}
                Uname = 'x .* t';
                true1 = x .* t; true2 = t; true3 = z;
                t6 = x; t7 = z;
                
            case {7}
                Uname = 'x.^3 + t.^3';
                true1 = x.^3 + t.^3; true2 = 3 * x.^2; true3 = 6 * x;
                t6 = 3 * t.^2; t7 = 6 * t; t8 = 6 * o;
            
            case {6}
                Uname = 'x.^2 + t.^4';
                true1 = x.^2 + t.^4; true2 = 2 * x.^1; true3 = 2 * o;
                t6 = 4 * t.^3; t7 = 12 * t.^2; t8 = 24 * t; t9 = 24 * o;
            
            case {5}
                Uname = 'x.^2 + t.^3';
                true1 = x.^2 + t.^3; true2 = 2 * x.^1; true3 = 2 * o;
                t6 = 3 * t.^2; t7 = 6 * t; t8 = 6 * o;
            
            case {4}
                Uname = 'x.^3 + t.^2';
                true1 = x.^3 + t.^2; true2 = 3 * x.^2; true3 = 6 * x;
                t6 = 2 * t; t7 = 2 * o;
                
            case {3}
                Uname = 'x.^2 + t.^2';
                true1 = x.^2 + t.^2; true2 = 2 * x; true3 = 2;
                true4 = z; true5 = z;
                t6 = 2 * t; t7 = 2;
                
            case {2}
                Uname = 'x.^2 + t';
                true1 = x.^2 + t; true2 = 2 * x; true3 = 2;
                true4 = z; true5 = z;
                t6 = o; t7 = z;
                
            case {1}
                Uname = 'x + t';
                true1 = x+t; true2 = o; true3 = z;
                true4 = z; true5 = z;
                t6 = o;
                
            case {0}
                Uname ='1';
                true1 = o;         true2 = z;           true3 = z;
                true4 = z;         true5 = z;

            case {-1}
                Uname ='0';
                true1 = z;         true2 = z;           true3 = z;
                true4 = z;         true5 = z;
        end
        
    otherwise
        switch Uno
        case {1003}
            Uname ='SP pg 390 ex 10.3 -- boundary layer at 0';
        %from Ascher+M+R pg 390 ex 10.3 -- boundary layer at 0, when 0 <~ eee
        %-1, 2+cos(pi*x), eee, [0, 1]
        true1 = cos(pi*x) - exp(-3*x/eee); % + O(eee^2)
        true2 = -pi*sin(pi*x) + 3*exp(-3*x/eee)/eee;
        true3 = -pi*pi*cos(pi*x) - 9*exp(-3*x/eee)/eee^2;

        case {1012}
            Uname ='SP pg 442 ex 10.12 -- boundary layers at -1 and 1';
        %from Ascher+M+R pg 442 ex 10.12 -- boundary layers at -1 and 1, when 0 <~ eee
        %-0.5, -x, eee, [-1, 1]
        true1 = exp(-(x+1)/eee) + 2*exp((x-1)/eee); % + O(eee)
        true2 = (-exp(-(x+1)/eee) + 2*exp((x-1)/eee))/eee;
        true3 = true1/eee^2;

        case {1092}
            Uname ='SP pg 371 ex 9.2 -- interior schock layer at 0';
        %from Ascher+M+R pg 371 ex 9.2 -- interior schock layer at 0, when 0 <~ eee
        %0, x, eee, [-1, 1]
        t = erf(1/sqrt(2*eee));
        true1 = cos(pi*x) + erf(x/sqrt(2*eee))/t;
        w = exp(-x.^2/(2.*eee))*sqrt(2)/(sqrt(eee).*sqrt(pi).*t);
        true2 = -sin(pi.*x).*pi + w;
        true3 = -cos(pi.*x)*pi.^2 - w.*x/eee;

        case {1004}
            Uname ='SP pg 394 ex 10.4 -- interior schock layer at 0';
        %from Ascher+M+R pg 394 ex 10.4 -- interior schock layer at 0, when 0 <~ eee
        % 0, 2*x, eee, [-1, 1]
        true1 = erf(x/sqrt(eee));
        true2 = 2*exp(-x.^2/eee)/sqrt(pi*eee);
        true3 = -2*x.*true2/eee;

        case {913}
            Uname ='PB ?? boundary?? layer at 0';
        %from Carey+Dinh PB 3 -- layer at 0, when etaC > 0 and R < 0.
        true1 = -(x + etaC).^R + (etaC^R*(1-x) + (1+etaC)^R*x);
        true2 = -R*(x + etaC)^(R-1) - etaC^R + (1+etaC)^R;
        true3 = -R*(R-1)*(x + etaC)^(R-2);

        case {911}
            Uname ='PB ?? boundary layer at 1';
        %from Carey+Dinh PB 1 -- boundary layer at 1, when 0 < etaC << 1
        w = sinh(1/etaC);
        true1 = sinh(x/etaC)/w;
        true2 = cosh(x/etaC)/w/etaC;
        true3 = true1/etaC^2;

        case {906}
            Uname ='PB 6 -- nonlinear';
        %PB 6 in [320] from Russell and Shampine -- nonlinear
        c = 1.336055694906108149;
        true1 = -log(2) + 2*log(c*sec(c*(x-.5)/2));
        true2 = c*tan(c*(x-.5)/2);
        true3 = (c^2 + true2.^2)/2;
        %true3 = c^2/2*(1 + true2.^2); % wrong -- old

        case {905}
            Uname ='PB 5 -- interior layer at mu';
        %PB 5 in [320] from Carey and Dinh PB 2 -- interior layer at mu, when nu large
        %coefu = 0; coefux = -2*nu*(x-mu); coefuxx = -1/nu - nu*(x-mu)^2;
        w = 1 + nu^2*(x-mu).^2;
        r = atan(nu*(x-mu)) + atan(nu*mu);
        true1 = (1-x).*r;
        true2 = -r+(1-x)*nu./w;
        true3 = -2*nu./w -2*(1-x)*nu^3.*(x-mu)./w.^2;

        case {904}
            Uname ='PB 4 -- boundary layers at 0 and 1';
        %PB 4 in [320]
        %adapted from Celia+Gray -- boundary layers at 0 and 1, when etaB large
        %coefu =-1; coefux = 1; coefuxx = 1;
        r = etaB/etaA; p = 1+r*x; q = 1+r*(1-x); w = log(p); v = log(q);
        c = 1/(log(1+r))^2;
        true1 = c*w.*v;
        true2 = c*r.*(v./p - w./q);
        true3 = -c*r.^2.*(v./p.^2 + 2./(p.*q) + w./q.^2);

        case {903}
            Uname =['PB 3 -- boundary layer at 0, etaB = ' num2str(etaB)];
        %PB 3 in [320] from Celia+Gray pg 180 -- boundary layer at 0, when etaB large
        %coefu = 0; coefux = etaB; coefuxx = etaA+etaB*x;
        r = etaB/etaA; c = log(1+r); v = 1+r*x;
        true1 = log(v)/c; true2 = r./v/c; true3 = -r^2./v.^2/c;
        true4 = 2*r^3./v.^3/c; true5 = -6*r^4./v.^4/c;
        %true1 = log(1+etaB*x/etaA)/log(1+etaB/etaA);
        %true2 = etaB/(log(1+etaB/etaA)*(1+etaB*x/etaA)*etaA);
        %true3 = -(etaB/etaA)^2/(log(1+etaB/etaA)*(1+etaB*x/etaA)^2);

        case {202}
            Uname = '1./(1 + 10*x)';
            true1 =  1./(1 + 10*x);
            true2 =    -10./(1 + 10*x).^2;
            true3 =    200./(1 + 10*x).^3;
            true4 =  -6000./(1 + 10*x).^4;
            true5 =-240000./(1 + 10*x).^5;

        case {201}
            Uname = 'sqrt(1 - x.^2) -- Arc';
            v = 1 - x.^2;
            true1 = sqrt(v);
            true2 = -x./sqrt(v);
            true3 =   -x.^2 ./ v.^(3/2) -   1./sqrt(v);
            true4 = -3*x.^3 ./ v.^(5/2) - 3*x./ v.^(3/2);
            true5 =-15*x.^4 ./ v.^(7/2) -18*x.^2 ./ v.^(5/2) - 3*x./ v.^(3/2);

        case {200}
            Uname = '1 ./ (1 + 25*x.^2) -- Runge 25';
            true1 = 1 ./ (1 + 25*x.^2);
            true2 = -50*x ./ (1 + 25*x.^2).^2;
            true3 = 5000*x.^2 ./ (1 + 25*x.^2).^3 - 50 ./ (1 + 25*x.^2).^2;
            true4 = -750000*x.^3 ./ (1 + 25*x.^2).^4 + 15000*x ./ (1 + 25*x.^2).^3;
            true5 = 150000000*x.^4 ./ (1 + 25*x.^2).^5 ...
                  - 4500000*x.^2 ./ (1 + 25*x.^2).^4 + 15000 ./ (1 + 25*x.^2).^3;

        case {199}
            Uname = '1 ./ (1 + x.^2) -- Runge 1';
            true1 =  1 ./ (1 + x.^2);
            true2 = -2*x./(1 + x.^2).^2;
            true3 =  8*x.^2./(1+x.^2).^3 - 2/(1+x.^2).^2;
            true4 =-48*x.^3./(1+x.^2).^4 + 24*x./(1+x.^2).^3;
            true5 =384*x.^4./(1+x.^2).^5 -288*x.^2./(1+x.^2).^4 + 24/(1+x.^2).^3;

        case {196}
            Uname =   'x^(19/2) -       3 *x^(17/2) +      3 *x^(15/2) -      x^(13/2)';
        % u = 0, u' = 0, u'' = 0 BCs in (0, 1), u in C^6
        true1=         x.^(19/2)-       3 *x.^(17/2)+      3 *x.^(15/2)-      x.^(13/2);
        true2=   19/2 *x.^(17/2)-    51/2 *x.^(15/2)+   45/2 *x.^(13/2)-13/2 *x.^(11/2);
        true3=  323/4 *x.^(15/2)-   765/4 *x.^(13/2)+  585/4 *x.^(11/2)-143/4 *x.^(9/2);
        true4= 4845/8 *x.^(13/2)-  9945/8 *x.^(11/2)+ 6435/8 *x.^(9/2)-1287/8 *x.^(7/2);
        true5=62985/16*x.^(11/2)-109395/16*x.^(9/2) +57915/16*x.^(7/2)-9009/16*x.^(5/2);

        case {144}
            Uname =    'x^(13/2) -    2 * x^(11/2) +        x^(9/2)';
        % u = 0, u' = 0 BCs in (0, 1), u in C^4
        true1 =         x.^(13/2) -    2 * x.^(11/2) +        x.^(9/2);
        true2 =  13/2 * x.^(11/2) -   11 * x.^(9/2)  +  9/2 * x.^(7/2);
        true3 = 143/4 * x.^(9/2)  - 99/2 * x.^(7/2)  + 63/4 * x.^(5/2);
        true4 =1287/8 * x.^(7/2)  -693/4 * x.^(5/2)  +315/8 * x.^(3/2);
        true5 =9009/16* x.^(5/2) -3465/8 * x.^(3/2)  +945/16* x.^(1/2);

        %a = 3;
        %true1 = x^a;       true2 = a*x^(a-1);     true3 = a*(a-1)*x^(a-2);

        case {175}
            Uname ='x.^(15/2)';
        true1 = x.^(15/2);  true2 = 15/2*x.^(13/2); true3 = 195/4*x.^(11/2);
        true4 =2145/8*x.^(9/2); true5 =19305/16*x.^(7/2);

        case {165}
            Uname ='x.^(13/2)';
        true1 = x.^(13/2);  true2 = 13/2*x.^(11/2); true3 = 143/4*x.^(9/2);
        true4 =1287/8*x.^(7/2); true5 =9009/16*x.^(5/2);

        case {155}
            Uname ='x.^(11/2)';
        true1 = x.^(11/2);  true2 = 11/2*x.^(9/2);  true3 = 99/4*x.^(7/2);
        true4 =693/8*x.^(5/2); true5 =3465/16*x.^(3/2);

        case {145}
            Uname ='x.^(9/2)';
        true1 = x.^(9/2);   true2 = 9/2*x.^(7/2);   true3 = 63/4*x.^(5/2);
        true4 =315/8*x.^(3/2); true5 =945/16*x.^(1/2);

        case {135}
            Uname ='x.^(7/2)';
        true1 = x.^(7/2);   true2 = 7/2*x.^(5/2);   true3 = 35/4*x.^(3/2);
        true4 =105/8*x.^(1/2); true5 =105/16*x.^(-1/2);

        case {125}
            Uname ='x.^(5/2)';
        true1 = x.^(5/2);   true2 = 5/2*x.^(3/2);   true3 = 15/4*x.^(1/2);
        true4 =15/8*x.^(-1/2); true5 =-15/16*x.^(-3/2);

        case {115}
            Uname ='x.^(3/2)';
        true1 = x.^(3/2);   true2 = 3/2*x.^(1/2);   true3 =  3/4*x.^(-1/2);
        true4 =-3/8*x.^(-3/2); true5 = 9/16*x.^(-5/2);

        case {109}
            Uname ='exp(x)';
        true1 = exp(x);    true2 = exp(x);      true3 = exp(x);
        true4 = exp(x);    true5 = exp(x);
        t6 = exp(x); t7 = exp(x); t8 = exp(x); t9 = exp(x);

        case {108}
            Uname ='exp(-x)';
        true1 = exp(-x);   true2 =-exp(-x);     true3 = exp(-x);
        true4 =-exp(-x);   true5 = exp(-x);
        t6 =-exp(-x); t7 = exp(-x); t8 =-exp(-x); t9 = exp(-x);

        case {107}
            Uname ='x^2*(1-x)^2*exp(x)';
        % u = 0, u' = 0 on x = 0, 1
        true1 = x.^2.*(1-x).^2.*exp(x);
        true2 = 2*x.*(1-x).^2.*exp(x) - 2*x.^2.*(1-x).*exp(x) + x.^2.*(1-x).^2.*exp(x);
        true3 = 2*(1-x).^2.*exp(x) - 8*x.*(1-x).*exp(x) + 4*x.*(1-x).^2.*exp(x) ...
              + 2*x.^2.*exp(x) - 4*x.^2.*(1-x).*exp(x) + x.^2.*(1-x).^2.*exp(x);
        true4 =-12*(1-x).*exp(x) + 6*(1-x).^2.*exp(x) + 12*x.*exp(x) ...
              - 24*x.*(1-x).*exp(x) + 6*x.*(1-x).^2.*exp(x) + 6*x.^2.*exp(x) ...
              - 6*x.^2.*(1-x).*exp(x) + x.^2.*(1-x).^2.*exp(x);
        true5 = 24*exp(x) - 48*(1-x).*exp(x) + 12*(1-x).^2.*exp(x) + 48*x.*exp(x) ...
              - 48*x.*(1-x).*exp(x) + 8*x.*(1-x).^2.*exp(x) + 12*x.^2.*exp(x) ...
              - 8*x.^2.*(1-x).*exp(x) + x.^2.*(1-x).^2.*exp(x);

        case {106}
            Uname ='x*(1-x)*exp(x)';
        % u = 0 on x = 0, 1
        true1 = x.*(1-x).*exp(x);
        true2 = (1-2*x).*exp(x) + true1;
        true3 =-2*exp(x) + (1-2*x).*exp(x) + true2;
        true4 =-4*exp(x) + (1-2*x).*exp(x) + true3;
        true5 =-6*exp(x) + (1-2*x).*exp(x) + true4;

        case {105}
            Uname ='exp(2*x) + (1-exp(2))*x.^2'; % u(0) = u(1), u'(0) = u'(1);
        true1 = exp(2*x) + (1-exp(2))*x.^2;
        true2 = 2*exp(2*x) + 2*(1-exp(2))*x;
        true3 = 4*exp(2*x) + 2*(1-exp(2));
        true4 = 8*exp(2*x);
        true5 =16*exp(2*x);

        case {102}
            Uname ='sin(x*pi/2)';
        true1 = sin(x*pi/2); true2 = pi/2*cos(x*pi/2); true3 =-pi^2/4*sin(x*pi/2);
        true4 =-pi^3/8*cos(x*pi/2); true5 = pi^4/16*sin(x*pi/2);
        %true1 = sin(2*pi*x); true2 = 2*pi*cos(2*pi*x); true3 =-(2*pi)^2*sin(2*pi*x);

        case {101}
            Uname ='cos(x)';
        true1 = cos(x);    true2 =-sin(x);      true3 =-cos(x);
        true4 = sin(x);    true5 = cos(x);

        case {100}
            Uname ='sin(x)';
        true1 = sin(x);    true2 = cos(x);      true3 =-sin(x);
        true4 =-cos(x);    true5 = sin(x);
        t6 = cos(x); t7 =-sin(x); t8 =-cos(x); t9 = sin(x);

        case {71}
            Uname ='x^8';
        true1 = x.^8;       true2 = 8*x.^7;       true3 = 56*x.^6;
        true4 =336*x.^5;    true5 =1680*x.^4;
        true6 = 6720*x.^3;  true7 =20160*x.^2;    true8 =40320*x; true9 = 40320;
        t6 = true6; t7 = true7; t8 = true8; t9 = true9;

        case {66}
            Uname ='x^7/210 - x^6/60 + x^5/60';
        %u^(3) = 0, u^(4) = 0 at x = 0, x = 1
        true1 = x.^7/210 -   x.^6/60 + x.^5/60;
        true2 = x.^6/30  -   x.^5/10 + x.^4/12;
        true3 = x.^5/5   -   x.^4/2  + x.^3/3;
        true4 = x.^4     - 2*x.^3    + x.^2;
        case {61}
            Uname ='x^7';
        true1 = x.^7;       true2 = 7*x.^6;       true3 = 42*x.^5;
        true4 =210*x.^4;    true5 = 840*x.^3;
        true6 = 2520*x.^2;  true7 = 5040*x;       true8 = 5040;
        t6 = true6; t7 = true7; t8 = true8;

        case {56}
            Uname ='x^6 - 3*x^5';
        %u'''' = 0 on x = 0, 1
        true1 = x.^6 - 3*x.^5; true2 = 6*x.^5 - 15*x.^4; true3 = 30*x.^4 - 60*x.^3;
        true4 = 120*x.^3 - 180*x.^2; true5 = 360*x.^2 - 360*x;
        case {54}
            Uname ='x.^3.*(x-1).^3'; %x^6 - 3*x^5 + 3*x^4 - x^3
        %u = 0, u' = 0, u'' = 0 on x = 0, 1
        true1 =     x.^6 -  3*x.^5 +  3*x.^4 -   x.^3;
        true2 =   6*x.^5 - 15*x.^4 + 12*x.^3 - 3*x.^2;
        true3 =  30*x.^4 - 60*x.^3 + 36*x.^2 - 6*x;
        true4 = 120*x.^3 -180*x.^2 + 72*x    - 6;
        true5 = 360*x.^2 -360*x    + 72;
        case {53}
            Uname ='x^6 - 2*x^5 + x^4';
        %u = 0, u' = 0, u''' = 0 on x = 0, 1
        true1 =     x.^4.*(x-1).^2;
        true2 =   6*x.^5 - 10*x.^4 + 4*x.^3;
        true3 =  30*x.^4 - 40*x.^3 + 12*x.^2;
        true4 = 120*x.^3 - 120*x.^2;
        true5 = 360*x.^2;
        case {52}
            Uname ='x^6 - x';
        %u = 0 on x = 0, 1
        true1 = x.^6 - x;  true2 = 6*x.^5 - 1;    true3 = 30*x.^4;
        true4 = 120*x.^3;  true5 = 360*x.^2;
        case {51}
            Uname ='x^6';
        true1 = x.^6;       true2 = 6*x.^5;       true3 = 30*x.^4;
        true4 =120*x.^3;    true5 = 360*x.^2;     true6 = 720*x;

        case {42}
            Uname ='x^5 - x';
        %u = 0 on x = 0, 1
        true1 = x.^5 - x;   true2 = 5*x.^4 - 1;   true3 = 20*x.^3;
        true4 = 60*x.^2;    true5 = 120*x;
        case {41}
            Uname ='x^5';
        true1 = x.^5;       true2 = 5*x.^4;       true3 = 20*x.^3;
        true4 = 60*x.^2;    true5 = 120*x;
        true6 = 120;        true7 = 0;            true8 = 0;       true9 = 0;
        t6 = true6; t7 = true7; t8 = true8; t9 = true9;

        case {33}
            Uname ='x^4 - 2*x^3 + x^2';
        %u = 0, u' = 0 on x = 0, 1
        true1 = x.^2.*(x-1).^2; true2 = 4*x.^3 - 6*x.^2 + 2*x; true3 = 12*x.^2 -12*x +2;
        true4 = 24*x - 12;  true5 = 24*o;
        case {32}
            Uname ='x^4 - x';
        %u = 0 on x = 0, 1
        true1 = x.^4 - x;   true2 = 4*x.^3 - 1;   true3 = 12*x.*x;
        true4 = 24*x;       true5 = 24*o;
        case {31}
            Uname ='x^4';
        true1 = x.^4;       true2 = 4*x.^3;       true3 = 12*x.*x;
        true4 = 24*x;       true5 = 24*o;

        case {23}
            Uname ='x.^3/3 - x.*x/2';
        %u' = 0 on x = 0, 1
        true1 = x.^3/3 - x.*x/2; true2 = x.*x - x;  true3 = 2*x-1;
        true4 = 2*o;       true5 = z;
        case {22}
            Uname ='x.^3 - x.*x';
        %u = 0 on x = 0, 1
        true1 = x^3 - x*x; true2 = 3*x*x - 2*x; true3 = 6*x-2;
        true4 = 6*o;       true5 = z;
        case {21}
            Uname ='x.^3';
        true1 = x.^3;      true2 = 3*x.*x;       true3 = 6*x;
        true4 = 6*o;       true5 = z;

        case {12}
            Uname ='x.*x - x';
        true1 = x.*x - x;  true2 = 2*x - 1;     true3 = 2*o;
        true4 = z;         true5 = z;
        case {11}
            Uname ='x.*x';
        true1 = x.*x;      true2 = 2*x;         true3 = 2*o;
        true4 = z;         true5 = z;

        case {1}
            Uname ='x';
        true1 = x;         true2 = o;           true3 = z;
        true4 = z;         true5 = z;

        case {0}
            Uname ='1';
        true1 = o;         true2 = z;           true3 = z;
        true4 = z;         true5 = z;

        case {-1}
            Uname ='0';
        true1 = z;         true2 = z;           true3 = z;
        true4 = z;         true5 = z;
        
        otherwise
            error(['truevd: no such function ' num2str(Uno)])
        end
end








