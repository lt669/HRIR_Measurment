function plotSingleShHarmFlat(phi, theta, r, rplot, density, h)

% DEFAULTS AND DECLARATIONS
if nargin < 4 % If resolution is not declared
	rplot = r;
end
if nargin < 5 % If grid density is not declared
	density = 50;
end
if nargin < 6 % If handle is not declared
	h = figure;
end

figure(h);
hold on;
view(35,20);

surface = surf(phi, theta, rplot, r);
set(surface, 'LineStyle','none', 'FaceColor', 'texturemap');
% Gradient Color
%     color(:, 3) = 1:-1/64:0;
%     color(:, 1) = 0: 1/64:1;
% Solid Color
    color =  [0, 0, 1;   % Negative Colour
              1, 0, 0];  % Positive Colour
colormap(color);

s1 = size(r, 1);
s2 = size(r, 2);

delta = round(size(phi, 1)/density);

Phi   = phi  (1:delta:s1, 1:delta:s2);
Theta = theta(1:delta:s1, 1:delta:s2);
R     = rplot(1:delta:s1, 1:delta:s2);

S1 = size(R, 1);
S2 = size(R, 2);

Phi(S1+1, :)    = phi(s1, 1:delta:s2);
Phi(1:S1, S2+1)    = phi(1:delta:s1, s2);
Phi(S1+1, S2+1) = phi(s1, s2);

Theta(S1+1, :)    = theta(s1, 1:delta:s2);
Theta(1:S1, S2+1)    = theta(1:delta:s1, s2);
Theta(S1+1, S2+1) = theta(s1, s2);

R(S1+1, :)    = rplot(s1, 1:delta:s2);
R(1:S1, S2+1)    = rplot(1:delta:s1, s2);
R(S1+1, S2+1) = rplot(s1, s2);

mesh = surf(Phi, Theta, R);
set(mesh, 'FaceAlpha',0);

caxis([min(r(:)), max(r(:))]);

end