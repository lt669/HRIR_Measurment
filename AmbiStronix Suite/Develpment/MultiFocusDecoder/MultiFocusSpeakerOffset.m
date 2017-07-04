% Function calculates the virtual sample offset required to align
% Multi-focus virtual loudspeakers with equivilant Single-Focus speakers to
% preserve ITD

% Positive Anticlockwise - (phi)
% Positive Northward - (theta)
% Positive interaural distance - right
% Negative interaural distance - left

% Handle inputs
aR = 0.11; % InterauralDistance / 2
b = 1.5; % Radius

aL = -aR; % InterauralDistance / 2
phi = (2 * pi) - phi; % In Radians
% theta = theta; % In Radians


x = asin(sin(phi).*cos(theta)) + pi/2;
DelR = b - ( ( sqrt(b^2 + (aR^2 * (cos(x).^2 - 1))) ) + ( aR * cos(x) ) );
DelL = b - ( ( sqrt(b^2 + (aL^2 * (cos(x).^2 - 1))) ) + ( aL * cos(x) ) );

sampDelR = round(DelR / 343 * 48000);
sampDelL = round(DelL / 343 * 48000);