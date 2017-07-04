function plotShHarm_vx(MaxL, resolution)
%
% Function Plots HArmonics Ontop of each other - Badly!
%
% Red - Positive
% Blue - Negative
%
% Written by Calum Armstrong, Department of Electronics, The University of
% York.
% Origionaly adapted from a script written by Mengliu Zhao, School of 
% Computing Science, Simon Fraser University. 2014/Dec/03.

% DEFAULTS AND DECLARATIONS
if nargin < 1 % If there are no input arguments
	MaxL = 2;
	resolution = 500;
end
if nargin < 2 % If there is 1 input argument
	resolution = 500;
end

syms fP_LM(u, l, m) % legandre function declaration

% DISCRETIZE SPHERE SURFACE
delta = 2*pi/resolution; % Angle between nodes
theta = -pi/2:delta:pi/2; % Elevation vector
phi = 0:delta:2*pi; % Azimuth vector
[phi,theta] = meshgrid(phi,theta); % 2-part mesh of nodes
                                                       
% Set figure background to white
sphHarm = figure('Color',[1 1 1]);

% Calculate maximum axis size required (L = Max; M = 0)
ah = sqrt(2*MaxL+1) + 0.1;

% PlOT SPHERICAL HARMONICS
for L = 0:MaxL % For each degree
    disp(sprintf('Calculating %d Degree', L)); % Inform console
    for M = -L:L % For each index

        % Legendre function over all elevation angles for 
        % current degree and |current index|
        M_L = abs(M) + L; % Pre-calculation
        fP_LM(u, l, m) = (1 / (2^l * factorial(l))) *...
                         ((1 - u^2)^(m/2)) *...
                         diff((u^2 - 1)^l,...
                         u ,...
                         M_L); % Define Legendre Formula

        outLength = length(theta(:, 1)); % Calculate output vector length
        vP_LM = zeros(outLength, 1); % Instantiate output
        [vP_LM(1), vP_LM(outLength)] = deal(not(boolean(M))); 
                                          % Calculate Legendre function (1)
        vP_LM(2:outLength-1) = fP_LM(sin(theta(2:outLength-1, 1)),...
                                     L, abs(M)); 
                                          % Calculate Legendre function (2)
        
        % Duplicate Legendre function vector to agree dimentions with 
        % spherical node mesh
        P_LM = repmat(vP_LM, [1, size(theta, 2)]);

        % Calculate partial normalization factor (N3D)
        % (not including sqrt(2) for |m| > 0)
        N_LM = sqrt((2*L+1)*factorial(L-abs(M))/factorial(L+abs(M)));

        % Calculate sperical harmonic function coefficients over mesh
        if M>0 % For positive orders...
            r = sqrt(2) * N_LM * P_LM .* cos(M*phi);
        elseif M==0 % For zeroth order...
            r = N_LM * P_LM;
        else % For negative orders...
            r = sqrt(2) * N_LM * P_LM .* sin(abs(M)*phi);
        end
        
        % Plot spherical harmonic
        figure(sphHarm); % Select figure
        hold on;
        
        plotPolarPattern(r, theta, phi, 2, ah) % plot
        %grid off
        %axis off
        
    end
end

