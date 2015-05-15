% function [T] = spdiag(d, n)
% T = spdiags([d.*ones(n, 1)], [0], n, n);
% n optional, n = length(d)

function [T] = spdiag(d, n)

[nr, nc] = size(d); if nr == 1, d = d'; end;
mx = max(nr, nc); mn = min(nr, nc);
if (nargin < 2)
    if (nr > 1) & (nc > 1)
        T = spdiags(diag(d), [0], mn, mn);
    else
        T = spdiags(d.*ones(mx, 1), [0], mx, mx);
    end
else
    if (mx == 1)
        T = spdiags([d*ones(n, 1)], [0], n, n);
    elseif (nr ~= n) & (nc ~= n)
        error(['spdiag: cannot form diagonal of ' num2str(n) ' elements'])
    elseif (nr > 1) & (nc > 1)
        T = spdiags(diag(d), [0], mn, mn);
    else
        T = spdiags(d.*ones(mx, 1), [0], mx, mx);
    end
end
