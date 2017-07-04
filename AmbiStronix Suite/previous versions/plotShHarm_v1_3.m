function plotShHarm_v1_3(M, plotSHflag, D, li, resolution)%, h)
%
% Script plots seperately the real and imaginary parts of spherical 
% harmonics up to degree = L (Ambisonics - 'order' = L). 
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
if nargin < 5 % If there are no input arguments
	resolution = 500;
end

% DISCRETIZE SPHERE SURFACE
delta = 2*pi/resolution; % Angle between nodes
theta = -pi/2:delta:pi/2; % Elevation vector
phi = 0:delta:2*pi; % Azimuth vector
                                                       
% Instantiate variable for plotting virtual microphone pickup pattern
vMic = zeros(length(theta), length(phi));
panPot = zeros(361);

% Set figures
if plotSHflag
    sphHarm = figure('Color',[1 1 1]);
end
virtualMic = figure;

% Calculate size of subplot array needed to plot spherical harmonics
width = 2*M + 1;
height = M + 1;

% Calculate maximum axis size required (L = Max; M = 0)
ah = sqrt(2*M+1) + 0.1;

% PlOT SPHERICAL HARMONICS
for m = 0:M % For each degree
    fprintf('Calculating %d Degree\n', m); % Inform console
    for i = -m:m % For each index

        [phiMesh, thetaMesh, rMesh] = calcShHarm(m, i, phi, theta);

        if plotSHflag
            % Plot spherical harmonic
            figure(sphHarm); % Select figure
            subplot(height,width,(width*m) + (M-m) + i+1+m) % Select subplot
                                                            % location     
            plotPolarPattern(rMesh, thetaMesh, phiMesh, 2, ah) % plot
            grid off
            axis off
        end
        
        % Scale and sum current harmonic for contribution to a potential 
        % virtual microphone polar pattern plot
        ACN = (m*(m+1)+i+1); 
        %scale = rMesh(resolution/4 + 1, 1); % Projection Decode (0,0)
        scale = D(li, ACN);
        
        vMic = vMic + rMesh * double(scale); % Sum scaled harmonic
    end
end

figure(virtualMic); % Select figure
[maxi, idx] = max(abs(vMic(:)));
plotPolarPattern(vMic, thetaMesh, phiMesh, 1, maxi) % Plot pattern
title('Virtual Microphone Pickup Pattern'); % Set title

planeEle = thetaMesh(idx);
panPotOffset = phiMesh(idx);
rot = 0:2*pi/360:2*pi;

phiPanPot(1:90)   = atan(tan(rot(1:90))/cos(planeEle)) + panPotOffset;
thetaPanPot(1:90) = asin(sin(planeEle)*cos(rot(1:90)));

phiPanPot(92:181)   = (pi) - atan(tan(pi - rot(92:181))/cos(planeEle)) + panPotOffset;
thetaPanPot(92:181) = -1 * asin(sin(planeEle)*cos(pi - rot(92:181)));

phiPanPot(182:270)   = (pi) + atan(tan(rot(182:270) - pi)/cos(planeEle)) + panPotOffset;
thetaPanPot(182:270) = -1 * asin(sin(planeEle)*cos(rot(182:270) - pi));

phiPanPot(272:361)   = (2*pi) - atan(tan(2*pi - rot(272:361))/cos(planeEle)) + panPotOffset;
thetaPanPot(272:361) = asin(sin(planeEle)*cos(2*pi - rot(272:361)));

[phiPanPot(91), phiPanPot(271)] = deal(pi/2 + panPotOffset, 3*pi/2 + panPotOffset);
thetaPanPot([91, 271]) = 0;

% phiPanPot = rem(phiPanPot, 2*pi);
% phiPanPot = round(phiPanPot / (2*pi) * (resolution)) + 1;
% thetaPanPot = round(((thetaPanPot / (pi/2)) * (resolution/4)) + (resolution/4)) + 1;
% 
% idx = sub2ind(size(vMic), thetaPanPot, phiPanPot);
% 
% panPot = vMic(idx);

% PlOT PANPOT
for m = 0:M % For each degree
    fprintf('Calculating %d Degree\n', m); % Inform console
    for i = -m:m % For each index

        [phiMesh, thetaMesh, rMesh] = calcShHarm(m, i, phiPanPot, thetaPanPot);
        
        % Scale and sum current harmonic for contribution to a potential 
        % virtual microphone polar pattern plot
        ACN = (m*(m+1)+i+1); 
        %scale = rMesh(resolution/4 + 1, 1); % Projection Decode (0,0)
        scale = D(li, ACN);
        
        panPot = panPot + rMesh * double(scale); % Sum scaled harmonic
    end
end

idx = sub2ind(size(panPot), 1:361, 1:361);
panPot = panPot(idx);

figure;
logic = panPot >= 0;
polar(rot(logic), abs(panPot(logic)), 'r');
hold on;
polar(rot(not(logic)), abs(panPot(not(logic))), 'b');
end

