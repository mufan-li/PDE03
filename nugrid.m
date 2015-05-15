% non-uniform grid for price
function [gridx] = nugrid(gridxu,ax,bx)
    global GridC etaB K Smax;
    m = length(gridxu)-1;
    
    switch GridC
        case {21}
            %Black Scholes Grid 2
            % a = 4; % a=4 sets grid 21 ~ grid 20
            % c1 = asinh((Smax-K)/a);
            % c2 = asinh((0-K)/a);
            % gridx = K + a*sinh(c1*gridxu/Smax+c2*(1-gridxu/Smax));

            % centered strike
            a0 = 20;
            a = grid21a(gridxu,K,Smax,a0);
            gridx = grid21(gridxu,K,Smax,a);
            
        case {20}
            %Black Scholes Grid 1
            a = (ceil(0.38*m)+0.5)/m;
            b = fzero(@(z) sinh(z*(1-a))/sinh(z*a)-(Smax-K)/K,10);
            gridx = (1 + sinh(b.*(gridxu/Smax-a))/sinh(b*a) )*K;
            
        case {10} % Inverse Function, normalized
            gridxu_temp = gridxu / bx * (log(1+etaB*bx)/log(1+etaB));
            gridx = ((1+etaB) .^ gridxu_temp - 1) / etaB;

        case {9} % test with log spacing
            logSpaceC = 0.001;
            gridx = exp(gridxu / bx * (log(bx + logSpaceC) - ...
                log(logSpaceC)) + log(logSpaceC)) - logSpaceC;
 
        case {8} % test with log spacing
            logSpaceC = 0.01;
            gridx = exp(gridxu / bx * (log(bx + logSpaceC) - ...
                log(logSpaceC)) + log(logSpaceC)) - logSpaceC;
           
        case {7} % test with log spacing
            logSpaceC = 0.1;
            gridx = exp(gridxu / bx * (log(bx + logSpaceC) - ...
                log(logSpaceC)) + log(logSpaceC)) - logSpaceC;

        case {6} % test with log spacing
            logSpaceC = 1000;
            gridx = exp(gridxu / bx * (log(bx + logSpaceC) - ...
                log(logSpaceC)) + log(logSpaceC)) - logSpaceC;

        case {5} % test with log spacing
            logSpaceC = 100;
            gridx = exp(gridxu / bx * (log(bx + logSpaceC) - ...
                log(logSpaceC)) + log(logSpaceC)) - logSpaceC;

        case {4} % test with log spacing
            logSpaceC = 10;
            gridx = exp(gridxu / bx * (log(bx + logSpaceC) - ...
                log(logSpaceC)) + log(logSpaceC)) - logSpaceC;

        case {3} % test with log spacing
            logSpaceC = 1;
            gridx = exp(gridxu / bx * (log(bx + logSpaceC) - ...
                log(logSpaceC)) + log(logSpaceC)) - logSpaceC;
            
        case {2}
            % gridx = ((1+etaB).^gridxu - 1)/etaB; 
            gridx = ((1+etaB).^gridxu - 1)*2/((1+etaB)^bx - 1);
            
        case {1}
            gridx = (exp(gridxu)-1) / (exp(bx)-1) * bx; % Non-Uniform grids
  
        otherwise
            gridx = gridxu;
    end
end

function [a] = grid21a(gridxu,K,Smax,a0)
    ng = 2e3; % # of guesses
    plk.x = zeros(ng+1,1);
    plk.y = zeros(ng+1,2);
    
    for i = 1:ng+1
        a = a0 * (1+(i-ng/2-1)/(ng/2*10));% restrict to 10%
        gridx = grid21(gridxu,K,Smax,a);
        nk = find(max(gridx-K,0),1);
        yk = abs(gridx(nk-1:nk)-K);
        plk.x(i) = a/a0;
        plk.y(i,:) = yk;
    end
    
    [~,ik] = sort(abs(plk.y(:,1)-plk.y(:,2)));
    [~,ik2] = min(abs(plk.x(ik(1:10))-1));
    a = plk.x(ik(ik2))*a0;
%     plot(plk.x,plk.y);
%     disp(a);
end

function [gridx] = grid21(gridxu,K,Smax,a)
%     a = 4; % a=4 sets grid 21 ~ grid 20
    c1 = asinh((Smax-K)/a);
    c2 = asinh((0-K)/a);
    gridx = K + a*sinh(c1*gridxu/Smax+c2*(1-gridxu/Smax));
end