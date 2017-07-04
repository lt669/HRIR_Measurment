syms phi

% A = [a, b, c; c, a, b; b, c, a]
% 
% 
% B = subs(A, a, 5)
K = 16;
Order = 3;

z90 = sym(zeros(K));
[z90(1, 1), diag, antidiag] = deal(1);
for m = 1:Order
    temp = sym(zeros(2*m + 1));
    n = size(temp,1);
    
    a = cos(m*phi);
    b = sin(m*phi);
    
    [temp(1:n+1:end  ), diag    ] = deal([ a, diag,     a]);
    [temp(n:n-1:end-1), antidiag] = deal([-b, antidiag, b]);
    
    z90(m^2 + 1:(m+1)^2, m^2 + 1:(m+1)^2) = temp;
end