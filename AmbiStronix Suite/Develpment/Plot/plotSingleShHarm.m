function [phi, theta, r] = plotSingleShHarm(m, i, res, h)
%
% Script plots a single Spherical Harmonic of degree m and index i (N3D 
% normalisation) onto handle h with resolution res. Function is plotted 
% as a surface in spherical coordinates with the radius and colour of 
% points on the surface representing the value of the funtion 
% f(phi, theta). 
%
% Red - Positive
% Blue - Negative
%
% Written by Calum Armstrong, Department of Electronics, The University of
% York.
%
% Updates ---
% 09 / 02 / 2017: Initial coding

% DEFAULTS AND DECLARATIONS
if nargin < 3 % If resolution is not declared
	res = 500;
end
if nargin < 4 % If resolution is not declared
	h = figure;
end

syms fP_LM(u, m_, i_) % legandre function declaration

% DISCRETIZE SPHERE SURFACE
delta = 2*pi/res; % Angle between nodes
phi = 0:delta:2*pi; % Azimuth vector
theta = -pi/2:delta:pi/2; % Elevation vector

% Set figure background to white
set(h, 'color', [1 1 1]);

% Calculate maximum axis size required (L = Max; M = 0)
ah = sqrt(2*m+1) + 0.1;

[phi, theta, r] = calcShHarm(m, i, phi, theta);

% Plot spherical harmonic                                                   
plotPolarPattern(r, phi, theta, 1.5, ah) % plot
grid off
axis off

end