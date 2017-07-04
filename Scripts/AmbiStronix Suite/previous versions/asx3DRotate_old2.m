% Ambisonic SoundField Rotation

function B_ = asx3DRotate_old2(B, X, Y, Z, sequence)

% B - Ambisonic Format signal representing soundfield to be rotated. ROWS
% should correspond to Ambisonic Channels in ACN format.
% X - Angle of rotation about the x axis
% Y - Angle of rotation about the y axis
% Y - Angle of rotation about the z axis
% sequence - Order / sequence of rotatins...
%           (xyz, xzy, yzx, yxz, zyx, zxy)
%
% COORDINATE SYSTEM !
%
%                     z                     x - Positive is Forward
%                     ^                     y - Positive is Left
%                     | <           x       z - Positive is Up
%                     |  \        ^         
%                    _|_/    __ /           Positive rotations about axes 
%                     |       / \           are in the clockwise 
%                     |     /   |           direction when facing the
%                     |   /  <-/            positive direction of that 
%         /\          | /                   axis
%   y <-----|---------|---------------
%           /       / |
%         </      /   |
%               /     |
%             /       |
%           /         |

K = 4;
Order = 1;
B = [1;1;1;1];
X = pi;
Y = pi;
Z = pi;

% %% Find Order of Ambisonics
% K = size(B, 1);
% order = sqrt(K) - 1;
% if not(isinteger(order))
%     error('Bad number of Ambisonic Channels');
% elseif order > 21 || order < 1
%     error('Order out of avaiable range');
% end

%% Load y(90) rotaion matricies
currentFile = which('asx3DRotate_old2.m'); 
folder = fileparts(currentFile);
files = dir([folder '\Ry90_CGs']);
for i=3:length(files)
    load(files(i).name)
end

%% Build y(-90) rotation matrix from data
y_90 = zeros(K);
y_90(1, 1) = 1;
for m = 1:Order
   y_90(m^2 + 1:(m+1)^2, m^2 + 1:(m+1)^2) = eval(sprintf('Ry90_%02d', m));
end

%% Build symbolic z(phi) rotation matrix
z = sym(zeros(K));
syms phi;
[z(1, 1), diag, antidiag] = deal(1);
for m = 1:Order
    temp = sym(zeros(2*m + 1));
    n = size(temp,1);
    
    a = cos(m*phi);
    b = sin(m*phi);
    
    [temp(1:n+1:end  ), diag    ] = deal([ a, diag,     a]);
    [temp(n:n-1:end-1), antidiag] = deal([-b, antidiag, b]);
    
    z(m^2 + 1:(m+1)^2, m^2 + 1:(m+1)^2) = temp;
end

%% Calculate z(90)

z90 = subs(z, phi, pi/2);

%% Optimisation

modzy = y_90 \ z90;
modyx = y_90 * (z90 \ y_90);
modxend_ = y_90;

%% zyx
zZ = subs(z, phi, Z);
zY = subs(z, phi, Y);
zX = subs(z, phi, X);
                   
B_ = modxend_\ zX ... % Roll  (x)
     *modyx*   zY ... % Pitch (y)
     *modzy*   zZ ... % Yaw   (z)   
     *         B;
 
%%
f = @() test(X, Y, Z, z, B, y_90, modzy, modyx, phi);
timeit(f)
end
%%
function out = test(X, Y, Z, z, B, y_90, mod1, mod2, phi)
        zZ = subs(z, phi, Z);
        zY = subs(z, phi, Y);
        zX = subs(z, phi, X);

        out = y_90 \ zX * mod2 *... % Roll
                     zY * mod1 *... % Pitch
                     zZ        *... % Yaw
              B;
end