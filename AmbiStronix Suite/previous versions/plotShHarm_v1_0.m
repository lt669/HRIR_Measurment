function plotShHarm_v1_0(MaxL, resolution)

% Script plots seperately the real and imaginary parts of spherical harmonics up to degree = L
% (Ambisonics - 'order' = L). 
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
% Adapted from a script written by Mengliu Zhao, School of Computing 
% Science, Simon Fraser University. 2014/Dec/03.
%
% Updates ---
% 24 / 11 / 2016: Initial coding
% 29 / 11 / 2016: Swap to using positive / negative elevation co-ordinate
% system (sin -> cos, cos -> sin, for all theta related sums); 
% update / correct some comments
% 02 / 12 / 2016: Switch to calculating Associated Legendre Equations from
% first principles

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

total = zeros(length(phi(:, 1)), length(phi(1, :)));

% Set figure background to white
sphHarm = figure('Color',[1 1 1]);
virtualMic = figure;

% Calculate size of subplot array needed to plot spherical harmonics
width = 2*MaxL + 1;
height = MaxL + 1;

% Calculate maximum axis size required (L = Max; M = 0)
ah = sqrt(2*MaxL+1) + 0.1;

% PlOT SPHERICAL HARMONICS
for L = 0:MaxL % For each degree
    disp(sprintf('Calculating %d Degree', L));
    for M = -L:L % For each order

        % Legendre function over all elevation angles for 
        % current degree and |current index|
        M_L = abs(M) + L; % Pre-calculation
        fP_LM(u, l, m) = (1 / (2^l * factorial(l))) *...
                         ((1 - u^2)^(m/2)) *...
                         diff((u^2 - 1)^l,...
                         u ,...
                         M_L); % Define Legendre Formula

        outLength = length(theta(:, 1));             
                     
        vP_LM = zeros(outLength, 1); % Set output vector size

        [vP_LM(1), vP_LM(outLength)] = deal(not(boolean(M)));
        vP_LM(2:outLength-1) = fP_LM(sin(theta(2:outLength-1, 1)), L, abs(M)); % Calculate
        
        % Duplicate vector to agree dimentions with spherical node mesh
        P_LM = repmat(vP_LM, [1, size(theta, 2)]);

        % Calculate partial normalization factor (N3D)
        N_LM = sqrt((2*L+1)*factorial(L-abs(M))/factorial(L+abs(M)));

        % Calculate sperical harmonic function coefficients over mesh
        if M>0 % For positive orders...
            r = sqrt(2) * N_LM * P_LM .* cos(M*phi);
        elseif M==0 % For zeroth order...
            r = N_LM * P_LM;
        else % For negative orders...
            r = sqrt(2) * N_LM * P_LM .* sin(abs(M)*phi);
        end

%         if (L == 0) % FOR SUMMING HARMONICS
%         total = total + r.*r; % abs(r) isnt quite so neat...
%         end
%         
%         if (L == 1) && (M == 1) % FOR SUMMING HARMONICS
%         total = total + r.*r; % abs(r) isnt quite so neat...
%         end
        
        scale = r(resolution/4 + 1, 1);
        total = total + r * scale;
        figure(virtualMic);
        plotPolarPattern(total, theta, phi)
        
        % Convert to ambisonic-cartisian co-ordinate system
        x = abs(r) .* cos(theta) .* cos(phi); 
        y = abs(r) .* cos(theta) .* sin(phi);
        z = abs(r) .* sin(theta);
        colourCode = double(r>=0); % Make note of positive / negative

        % Plot spherical harmonic
        figure(sphHarm);
        
        subplot(height,width,(width*L) + (MaxL-L) + M+1+L) % Select subplot
                                                           % location
        h = surf(x,y,z,colourCode); % Plot

        % ADJUST CAMERA VIEW
            % IDENTICAL VIEWING ANGLE FOR EACH HARMONIC
                view(290,30)
                camzoom(2)

            % RELATIVE VIEWING ANGLE FROM CENTRAL VIEWING POINT
%                 hCenter = ceil((MaxL/2));
%                 hView = 10 * (L - hCenter) + 5;
%                 wCenter = MaxL + 1;
%                 wView = 270 + (10 * (wCenter - (MaxL-L + M+1+L))); 
%                 view(wView,hView);
%                 camzoom(1.5)

        % Apply light to image
        camlight right

        % Re-size axes
        axis([-ah ah -ah ah -ah ah])
        xlabel('x')

        % Colour plots
        colormap(  [0, 0, 1;  % Negative Colour
                    1, 0, 0]) % Positive Colour

        % Hide lines
        set(h, 'LineStyle','none')

        % Remove extras
        grid off
        axis off
        % text(0,-ah,-ah,[num2str(L), ', ', num2str(M)])
    end

end

