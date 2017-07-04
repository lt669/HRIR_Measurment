function plotShHarm_v1_4(M, plotSHflag, D, li, norm, resolution)

% Version 1.4 Feb 2017
%
% Graphical Output From f(Maximum Degree,
%                         Output Flag,
%                         Ambisonic Decode Matrix, 
%                         Speaker number,
%                         * Nomalisation,
%                         * Resolution)
%
%
% Script has two functions. Firstly, to calculate and plot (if 
% 'plotSHFlag' == 1) the individual Spherical Harmonic (SH) functions up 
% to order 'M' with normalisation 'norm' (default N3D). If 
% 'plotSHFlag' == 0 the plotting functionality is supressed, however, 
% the SHs are still calculated. SH coefficients are calculated at 
% resolution 'resolution' (default 500). 
% Secondly, the script will output the 3-Dimentional Virtual Microphone 
% pickup pattern generated for speaker 'li' given the Ambisonic decode 
% matrix 'D' and the SHs calculated. This will be complimented by an 
% aproximation of a 2-Dimentional polar plot depicting a plane of 
% symetry for the 3D  pattern.
%
% Plot Key:
% Red - Positive
% Blue - Negative
%
% An ambisonic co-ordinate system is used where phi is azimuth; 0degrees 
% is straight ahead; positive is anticlockwise. Theta is elevation; 
% 0degrees is on the horizontal axis; positive is upwards.
%
% Written by Calum Armstrong, Department of Electronics, The University 
% of York.
% Origionaly adapted from a script written by Mengliu Zhao, School of 
% Computing Science, Simon Fraser University. 2014/Dec/03.
%
% Updates ---
% 24 / 11 / 2016: Initial coding.
% 29 / 11 / 2016: Swap to using positive / negative elevation .
%   co-ordinate system (sin -> cos, cos -> sin, for all theta related 
%   sums); update / correct some comments.
% 02 / 12 / 2016: Switch to calculating Associated Legendre Equations 
%   from first principles.
% 08 / 12 / 2016: Introduced functionality to plot virtual microphone 
%   polar patterns and refactored making use of external plotting 
%   function.
% 14 / 02 / 2017: Include functionality to output 2D PanPot plot and
%   eddited input options for plotting virtual microphone pickup 
%   patterns (e.g. with the use of decode matrix 'D'); Included 
%   functionality to edit normalisation.

%% DEFAULTS AND DECLARATIONS

    if nargin < 5 % If no normalisation input is given
        norm = 'n3d';
    end
    if nargin < 6 % If no resolution input is given
        resolution = 500;
    end

%% SETUP

% Discretize sphere surface
    delta = 2*pi/resolution; % Angle between nodes
    theta = -pi/2:delta:pi/2; % Elevation vector
    phi = 0:delta:2*pi; % Azimuth vector
                                                       
% Instantiate variable for plotting virtual microphone pickup pattern
    vMic = zeros(length(theta), length(phi));
    panPot = zeros(361);

% Generate figure
    if plotSHflag
        figure('Color',[1 1 1]);
    end

% Calculate size of subplot array needed to plot spherical harmonics
    width = 2*M + 1;
    height = M + 1;

% Calculate maximum axis size required (L = Max; M = 0)
    ah = sqrt(2*M+1) + 0.1;

%% CALCULATE AND PlOT SPHERICAL HARMONICS AND VIRTUAL MIC. PATTERN

    for m = 0:M % For each degree
        fprintf('Calculating %d Degree\n', m); % Inform console
        for i = -m:m % For each index

            % Calculate Spherical Harmonics coefficients
            [phiMesh, thetaMesh, rMesh] = ...
                                 iasxShHarmCoef(m, i, phi, theta, norm);

            % Plot spherical harmonic if selected
            if plotSHflag
                subplot(height,width,(width*m) + (M-m) + i+1+m);
                                               % Select subplot location     
                plotPolarPattern(rMesh, phiMesh, thetaMesh, 2, ah); 
                                                                  % plot
                grid off;
                axis off;
            end

            % Scale and sum current harmonic for contribution to a 
            % virtual microphone polar pattern plot
            ACN = (m*(m+1)+i+1); % Calculate Ambisonic Channel Number
            %scale = rMesh(resolution/4 + 1, 1); % Projection Decode
            scale = D(li, ACN); % Extract scaling factor from Decode 
                                % Matrix
            vMic = vMic + rMesh * double(scale); % Sum scaled harmonic
        end
    end

% Plot 3D Virtual Microphone Pickup pattern
    figure; % New figure
    suptitle('Virtual Microphone Pickup Pattern'); % Set title of figure
    subplot(1,2,1);
    [maxi, idx] = max(abs(vMic(:))); % Find max value of vMic 
                                     % (=> direction)
    plotPolarPattern(vMic, phiMesh, thetaMesh, 1, maxi); % Plot pattern
    title('3D');


%% CALCULATE AND PLOT PANPOT

    fprintf('Calculating PanPot\n'); % Inform console

% Extract directional information of virtual microphone pattern
    pEle = thetaMesh(idx); % Plane Elevation
    Offset = phiMesh(idx); % Plane rotation (phi)

% Nasty geometrical maths that calculates the angles (phi, theta) that
% correspond to a full rotation around a plane of symetry of the virtual
% microphone pattern...
    rot = 0:2*pi/360:2*pi; % Rotation index
    phiPP(1:90)=atan(tan(rot(1:90))/cos(pEle))+Offset; % Phi Pan Pot
    thetaPP(1:90)=asin(sin(pEle)*cos(rot(1:90))); % Theta Pan Pot
    phiPP(92:181)=(pi)-atan(tan(pi-rot(92:181))/cos(pEle))+Offset;
    thetaPP(92:181)=-1*asin(sin(pEle)*cos(pi-rot(92:181)));
    phiPP(182:270)=(pi)+atan(tan(rot(182:270)-pi)/cos(pEle))+Offset;
    thetaPP(182:270)=-1*asin(sin(pEle)*cos(rot(182:270)-pi));
    phiPP(272:361)=(2*pi)-atan(tan(2*pi-rot(272:361))/cos(pEle))+Offset;
    thetaPP(272:361)=asin(sin(pEle)*cos(2*pi-rot(272:361)));
    [phiPP(91),phiPP(271)]=deal(pi/2+Offset,3*pi/2+Offset);
    thetaPP([91,271])=0;

% Calculate Spherical Harmonics coefficients at angles of interest
    for m = 0:M 
        for i = -m:m 
            [~, ~, rMesh] = iasxShHarmCoef(m, i, phiPP, thetaPP, norm);
            ACN = (m*(m+1)+i+1);
            %scale = rMesh(resolution/4 + 1, 1);
            scale = D(li, ACN);
            panPot = panPot + rMesh * double(scale);
        end
    end

% Extract relavent Spherical Harmonic coefficients from mesh matrix
    idx = sub2ind(size(panPot), 1:361, 1:361);
    panPot = panPot(idx);

% Plot PanPot
    subplot(1,2,2);
    logic = panPot >= 0;
    polar(rot(logic), abs(panPot(logic)), 'r'); % Plot positive values
    hold on;
    polar(rot(not(logic)), abs(panPot(not(logic))), 'b'); % Negatives
    title('2D');
    
end

