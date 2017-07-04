function plotShHarm_v1_2(D, Li, resolution)
%
% Script plots seperately the real and imaginary parts of spherical 
% harmonics up to degree = L (Ambisonics - 'order' = L). 
%
% Ambisonic co-ordinates, x -> Forward
%                         y -> left
%                         z -> up
%                         phi -> 0Degree = forward; anticlockwise is
%                         positive
%                         theta -> 0Degree = horizontal axis; up is
%                         positive
%
% N3D normalization
%
% Red - Positive
% Blue - Negative
%
% Written by Calum Armstrong, Department of Electronics, The University of
% York.
% Origionaly adapted from a script written by Mengliu Zhao, School of 
% Computing Science, Simon Fraser University. 2014/Dec/03.
%
% Updates ---
% 24 / 11 / 2016: Initial coding
% 29 / 11 / 2016: Swap to using positive / negative elevation co-ordinate
% system (sin -> cos, cos -> sin, for all theta related sums); 
% update / correct some comments
% 02 / 12 / 2016: Switch to calculating Associated Legendre Equations from
% first principles
% 08 / 12 / 2016: Introduced functionality to plot virtual microphone polar
% patterns and refactored making use of external plotting function

% DEFAULTS AND DECLARATIONS
if nargin < 3 % If there are no input arguments
	resolution = 500;
end

syms fP_LM(u, l, m) % legandre function declaration
fP_LM(u, l, m) = 1;

% DISCRETIZE SPHERE SURFACE
delta = 2*pi/resolution; % Angle between nodes
theta = -pi/2:delta:pi/2; % Elevation vector
phi = 0:delta:2*pi; % Azimuth vector
[phi,theta] = meshgrid(phi,theta); % 2-part mesh of nodes
                                                       
% Instantiate variable for plotting virtual microphone pickup pattern
vMic = zeros(length(phi(:, 1)), length(phi(1, :)));

% Set figure background to white
sphHarm = figure('Color',[1 1 1]);
virtualMic = figure;

MaxL  = sqrt(length(D(1, :)))-1;

% Calculate size of subplot array needed to plot spherical harmonics
width = 2*MaxL + 1;
height = MaxL + 1;

% Calculate maximum axis size required (L = Max; M = 0)
ah = sqrt(2*MaxL+1) + 0.1;

% PlOT SPHERICAL HARMONICS
for L = 0:MaxL % For each degree
    fprintf('Calculating %d Degree\n', L); % Inform console
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
        
        if M ~= 0 % Exception that otherwise throws zero divide error
            [vP_LM(1), vP_LM(outLength)] = deal(0); 
                                          % Calculate Legendre function (1)
            vP_LM(2:outLength-1) = fP_LM(sin(theta(2:outLength-1, 1)),...
                                     L, abs(M)); 
                                          % Calculate Legendre function (2)
        else
            vP_LM(1:outLength) = fP_LM(sin(theta(:, 1)), L, abs(M));
        end
        
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
        
%         % FOR FU-MA FIRST ORDER
%         if L == 0
%             r = 1/sqrt(2) * ones(length(phi(:, 1)), length(phi(1, :)));
%         elseif M == -1
%             r = sin(phi).*cos(theta);
%         elseif M == 0
%             r = sin(theta);
%         elseif M == 1
%             r = cos(phi).*cos(theta);
%         end

        % Plot spherical harmonic
        figure(sphHarm); % Select figure
        subplot(height,width,(width*L) + (MaxL-L) + M+1+L) % Select subplot
                                                           % location                                                     
        plotPolarPattern(r, theta, phi, 2, ah) % plot
        grid off
        axis off
        
        % Scale and sum current harmonic for contribution to a potential 
        % virtual microphone polar pattern plot
        ACN = (L*(L+1)+M+1);                                
        scale = D(Li, ACN);                            
        
        vMic = vMic + r * double(scale); % Sum scaled harmonic
        figure(virtualMic); % Select figure
        plotPolarPattern(vMic, theta, phi, 1, max(max(abs(vMic)))) 
                                                            % Plot pattern
        title('Virtual Microphone Pickup Pattern'); % Set title
        
    end
end

