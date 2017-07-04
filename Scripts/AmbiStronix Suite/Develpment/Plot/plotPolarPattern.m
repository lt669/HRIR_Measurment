function plotPolarPattern(r, phi, theta, zoom, size)

% Graphical Output From f(Coefficients,
%                         Azimuths, 
%                         Elevations,
%                         Camera zoom,
%                         Axes Size)
%
% Script plots a surface of coefficients given in spherical coordinates.
% Coefficients 'r' and coordinates: 'phi', 'theta' should be given as 
% mesh matricies as produced by the function meshgrid. Surface is 
% plotted on axes with range +/- 'size' and with camera zoom 'zoom'.
%
% An ambisonic co-ordinate system is used where phi is azimuth; 0degrees 
% is straight ahead; positive is anticlockwise. Theta is elevation; 
% 0degrees is on the horizontal axis; positive is upwards.
%
% Written by Calum Armstrong, Department of Electronics, The University 
% of York.
%
% Updates ---
% 08 / 12 / 2016: Initial coding

%% PLOT SURFACE

% Convert Ambisonic-Polar to Cartisian co-ordinates
    x = abs(r) .* cos(theta) .* cos(phi); 
    y = abs(r) .* cos(theta) .* sin(phi);
    z = abs(r) .* sin(theta);
    colourCode = double(r>=0); % Make note of polarity

    surface = surf(x,y,z,colourCode); % Plot surface

% Adjust camera view
    view(290,30)
    camzoom(zoom)

% Apply light to image
    camlight right

% Re-size axes
    axis([-size size -size size -size size])
    xlabel('x','FontWeight','bold')
    ylabel('y','FontWeight','bold')
    zlabel('z','FontWeight','bold')

% Colour plots
    colormap( [0, 0, 1;   % Negative Colour
               1, 0, 0] ) % Positive Colour

% Hide lines
    set(surface, 'LineStyle','none')

end