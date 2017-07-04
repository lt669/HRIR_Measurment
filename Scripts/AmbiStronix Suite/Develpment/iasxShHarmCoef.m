function [ phiMesh,...
           thetaMesh,...
           rMesh ] = iasxShHarmCoef( m, i, phi, theta, norm )

% Version 1.0 Feb 2017
%
% [Azimuths, Elevations, Coefficients] From f( Degree,
%                                              Order,
%                                              Azimuths, 
%                                              Elevations
%                                              * Nomalisation )
%
% Script calculates the coefficients 'rMesh' of the Spherical Harmonic:
% degree 'm', index 'i' over the [theta x phi] mesh matrix generated 
% from the inputs 'phi' and 'theta'. Normalization of 'norm' is applied. 
% Script also outputs the matricies 'phiMesh' and 'thetaMesh' whose data 
% correspond to the coefficients calculated. 
% i.e. Y_mi^sigma(phiMesh(x,y), thetaMesh(x,y)) = rMesh(x,y)
%
% For a descrete set of points, 'phi'=[a,b,c,d,e] 'theta'=[f,g,h,i,j] the
% Spherical Harmonic coefficients for these points only are given by the
% diagonal coefficients of mesh matrix 'rMesh'.
%
% An ambisonic co-ordinate system is used where phi is azimuth; 0degrees 
% is straight ahead; positive is anticlockwise. Theta is elevation; 
% 0degrees is on the horizontal axis; positive is upwards.
%
% Written by Calum Armstrong, Department of Electronics, The University 
% of York.
%
% Updates ---
% 10 / 02 / 2017: Initial coding

%% DEFAULTS AND DECLARATIONS

    if nargin < 5 % If normalisation is not declared
        norm = 'n3d'; % Default to N3D
    end

    syms fP_LM(u, m_,i_) % Declare Assosiated Legendre Polynomial (ASP)
    
    [phiMesh,thetaMesh] = meshgrid(phi,theta); % Generate Node Mesh

%% CALCULATE COEFFICIENTS

% Calculate ASP for all 'theta' for degree 'm' and |index 'i'|
    sigma = sign(i); % Pre-calculation / Derivations
    i = abs(i);      % |
    m_i =  m + i;    % /

    fP_LM(u, m_, i_) = (1 / (2^m_ * factorial(m_))) *...
                       ((1 - u^2)^(i_/2)) *...
                       diff((u^2 - 1)^m_,...
                       u ,... 
                       m_i); % Define ALP formula

    if i ~= 0 % Exception that otherwise throws zero divide error *            
        P_LM = zeros(length(theta), 1); % Set output to zeros
        logic = abs(sin(theta)) == 1;
        P_LM(not(logic)) = double(fP_LM(sin(theta(not(logic))), m, i));                                   
    else
        P_LM = double(fP_LM(sin(theta(:)), m, i));
    end

    % * fP_LM(1, m, !0) throws a zero divide error. In simple terms this
    % is when the ALP should equal 0 on the Z-axis. To
    % solve this, these '0' results are entered manually whilst the
    % remaining results are calculated by MatLab.

% Duplicate ALP vector to agree dimentions with Node Mesh
    P_LM = repmat(P_LM, [1, length(phi)]);

% Calculate partial normalization factor (N3D)(not including sqrt(2) for 
% |i| > 0)
    switch norm
        case 'n3d'
            N_mi = sqrt((2*m+1)*factorial(m-i)/factorial(m+i));
        case 'sn3d'
            N_mi = sqrt(factorial(m-i)/factorial(m+i));
        case 'orthonormal'
            N_mi = sqrt((2*m+1)/(4*pi)*factorial(m-i)/factorial(m+i));
        otherwise
            error('unknown normalisation');
    end

% Calculate sperical harmonic function coefficients over mesh
    if sigma==1 % For positive index...
        rMesh = sqrt(2) * N_mi * P_LM .* cos(i*phiMesh);
    elseif sigma==0 % For zero index...
        rMesh = N_mi * P_LM;
    else % For negative index...
        rMesh = sqrt(2) * N_mi * P_LM .* sin(i*phiMesh);
    end
    
end