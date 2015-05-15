% Handles very general BCs given in the form of
% u(a,t) = ga(t)
% u(b,t) = gb(t)
% u(x,t0) = g0(x)

function [u_val] = DirechletBC(x,t)
    global BCno K Smin Smax Rf SigmaC IC_Offset;
    switch BCno
        case {'American Call'}
            if (abs(t)<1e-5)
                u_val = (x-K>0).*(x-K); %to return vector format
            elseif (abs(x-Smin)<1e-5)
                u_val = 0.*x.*t;
            elseif (abs(x-Smax)<1e-5)
                u_val = (x-K).*sign(t);
            else
                u_val = 0;
                disp('Point not on BC.');
            end
            
        case {'American Put'}
            if (abs(t)<1e-5)
                u_val = (K-x>0).*(K-x); %to return vector format
            elseif (abs(x-Smin)<1e-5)
                u_val = (K-x).*sign(t);
            elseif (abs(x-Smax)<1e-5)
                u_val = 0.*x.*t;
            else
                u_val = 0;
                disp('Point not on BC.');
            end
            
        case {'European Call + delta'}
            if (abs(t)<1e-5)
                u_val = EuroBls(x, K, Rf, IC_Offset, 0, SigmaC); %to return vector format
            elseif (abs(x-Smin)<1e-5)
                u_val = 0.*x.*(t+IC_Offset);
            elseif (abs(x-Smax)<1e-5)
                u_val = x-K.*exp(-Rf.*(t+IC_Offset));
            else
                u_val = 0;
                disp('Point not on BC.');
            end
            
        case {'European Put + delta'}
            if (abs(t)<1e-5)
                [ans, u_val] = EuroBls(x, K, Rf, IC_Offset, 0, SigmaC); %to return vector format
            elseif (abs(x-Smin)<1e-5)
                % note in this case t is transformed as tau=T-t
                u_val = K.*exp(-Rf.*(t+IC_Offset));
            elseif (abs(x-Smax)<1e-5)
                u_val = 0.*x.*(t+IC_Offset);
            else
                u_val = 0;
                disp('Point not on BC.');
            end
            
        case {'European Call'}
            if (abs(t)<1e-5)
                u_val = (x-K>0).*(x-K); %to return vector format
            elseif (abs(x-Smin)<1e-5)
                u_val = 0.*x.*t;
            elseif (abs(x-Smax)<1e-5)
                u_val = x-K.*exp(-Rf.*t);
            else
                u_val = 0;
                disp('Point not on BC.');
            end
            
        case {'European Put'}
            if (abs(t)<1e-5)
                u_val = (K-x>0).*(K-x); %to return vector format
            elseif (abs(x-Smin)<1e-5)
                % note in this case t is transformed as tau=T-t
                u_val = K.*exp(-Rf.*t);
            elseif (abs(x-Smax)<1e-5)
                u_val = 0.*x.*t;
            else
                u_val = 0;
                disp('Point not on BC.');
            end
            
        case {'PS2-4(iii)'}
            % hack for part 4(iii)(iv)
            % can be built for more general cases with error checks
            u_val = x ./ abs(x);
            u_val(isnan(u_val)) = 1;
    
        case {1}
            u_val = (x-Smax/2).^2 + t;
            
        otherwise
            % assumes (x,t) are only at IC and BC
            u_val = truevd(x,t);
    end
end