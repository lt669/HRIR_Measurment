function    r = ShHarmCoef_v1_0(L, M, phi, theta, norm)
%
% Version 1.0 Dec 2016
%
% Coefficient = f(degree, index, azimuth, elivation, normalisation)
%
% norm:  N3D        (n3d)          - default
%        SN3D       (sn3d)
%        Orthonomal (orthonormal)
%
% Script calculates the spherical harmonic coefficient (r) for a given 
% degree (L), index (M) and angle (phi, theta) where L >= 0 and 
% |M| <= L. Coefficients are calculated for N3D normalisation by default   
% although Orthonomal and SN3D normalisations are available as options 
% (norm).
%
% NOTE: Angles should be given in radians!
%
% An ambisonic co-ordinate system is used where phi is azimuth; 0degrees 
% is straight ahead; positive is anticlockwise. Theta is elevation; 
% 0degrees is on the horizontal axis; positive is upwards. 
%
% Written by Calum Armstrong, Department of Electronics, The University 
% of York.
%
% Updates ---
% 05 / 12 / 2016: Initial coding
% 06 / 12 / 2016: Fpdates to comments / tidy up
% 13 / 12 / 2016: Finalisation of version 1.0

    % Set default normalisation if none is given
    if nargin < 5
        norm = 'n3d'; % Default normalisation is N3D 
    end

    % Check input variables are of viable values, else throw error
    if L < 0
        error('        L out of range: L > 0');
    elseif abs(M) > L
        error('        M out of range: |M| > L');
    elseif abs(theta) > (pi/2)
        error('        theta out of range: |theta| > (pi/2)');
    end

    % Declare Legendre function
    syms fP_LM(u, l, m)

    % Define Legendre Formula
    M_L = abs(M) + L; % Pre-calculation
    fP_LM(u, l, m) = (1 / (2^l * factorial(l))) *...
                     ((1 - u^2)^(m/2)) *...
                     diff((u^2 - 1)^l,...
                     u ,...
                     M_L); 

    % Calculate Legendre function for given degree, index, elevation
    if M ~= 0 && abs(sin(theta)) == 1 % Exception that otherwise throws  
                                      % a zero divide error
        P_LM = 0;
    else
        P_LM = fP_LM(sin(theta), L, abs(M)); % General calculation
    end

    % Calculate partial normalization factor (not including sqrt(2) for
    % m != 0)
    switch norm
        case 'n3d'
            N_LM = sqrt((2*L+1)*...
                        factorial(L-abs(M))/factorial(L+abs(M)));
        case 'sn3d'
            N_LM = sqrt(factorial(L-abs(M))/factorial(L+abs(M)));
        case 'orthonomal'
            N_LM = sqrt((2*L+1)/(4*pi)*...
                        factorial(L-abs(M))/factorial(L+abs(M)));
        otherwise
            error('unknown normalisation');
    end

    % Calculate sperical harmonic function coefficients
    if M>0 % For positive indexes...
        r = sqrt(2) * N_LM * P_LM .* cos(M*phi);
    elseif M==0 % For zeroth index...
        r = N_LM * P_LM;
    else % For negative indexes...
        r = sqrt(2) * N_LM * P_LM .* sin(abs(M)*phi);
    end

    % Return coefficient as double
    r = double(r);

end