syms fP_LM(u, l, m)

L = 4;
M = 4;
delta = 2*pi/100; % Angle between nodes
theta = -pi/2:delta:pi/2; % Elevation vector


M_L = abs(M) + L;
fP_LM(u, l, m) = (1 / (2^l * factorial(l))) * ((1 - u^2)^(m/2)) * diff((u^2 - 1)^l, u , M_L);

if L == 0
    
    vP_LM = 1;
    
else
    
    vP_LM = zeros(length(theta), 1);
    
    for i = 1:length(theta)
        if abs(sin(theta(i))) == 1
            vP_LM(i) = 0;
        else
            vP_LM(i) = fP_LM(0.99999999, L, abs(M));
        end
    end
end

test = double(vP_LM);


P_LallM = legendre(L,sin(theta));
P_LM = (-1)^M * P_LallM(abs(M)+1,:)';