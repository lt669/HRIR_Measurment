function rot_zyx = asxRotMatrix_zyx(Order)

    K = (Order + 1)^2;

    %% Load y(90) rotaion matricies
    currentFile = which('asxrotMatrix_zyx.m'); 
    folder = fileparts(currentFile);
    files = dir([folder '/Ry90_CGs']);
    for i=3:length(files)
        load(files(i).name)
    end

    %% Build y(-90) rotation matrix from data
    y_90 = zeros(K);
    y_90(1, 1) = 1;
    for m = 1:Order
       y_90(m^2 + 1:(m+1)^2, m^2 + 1:(m+1)^2) = eval(sprintf('Ry90_%02d', m));
    end

    %% Build symbolic z(phi) rotation matrix
    z = sym(zeros(K));
    syms phi;
    [z(1, 1), diag, antidiag] = deal(1);
    for m = 1:Order
        temp = sym(zeros(2*m + 1));
        n = size(temp,1);

        a = cos(m*phi);
        b = sin(m*phi);

        [temp(1:n+1:end  ), diag    ] = deal([ a, diag,     a]);
        [temp(n:n-1:end-1), antidiag] = deal([-b, antidiag, b]);

        z(m^2 + 1:(m+1)^2, m^2 + 1:(m+1)^2) = temp;
    end
    
    %% Calculate z(90)

    z90 = subs(z, phi, pi/2);
    
    %% Optimisation
    
    syms Z Y X
    
    modzy   = y_90 \ z90;
    modyx   = y_90 * (z90 \ y_90);
    modxend_ = y_90;
    
    zZ = subs(z, phi, Z);
    zY = subs(z, phi, Y);
    zX = subs(z, phi, X);
    
    rot_zyx = modxend_\ zX ... % Roll  (x)
              *modyx*   zY ... % Pitch (y)
              *modzy*   zZ;    % Yaw   (z)
    
end