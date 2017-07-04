% Ambisonic SoundField Rotation

function B_ = asx3DRotate(B, X, Y, Z, sequence)

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


%% Preperation

% Load y(90) rotaion matricies
currentFile = which('asx3DRotate.m'); 
folder = fileparts(currentFile);
files = dir([folder '\Ry90_CGs']);
for i=3:length(files)
    load(files(i).name)
end

%% Find Order of Ambisonics
K = size(B, 1);
order = sqrt(K) - 1;
if not(isinteger(order))
    error('Bad number of Ambisonic Channels');
elseif order > 21 || order < 1
    error('Order out of avaiable range');
end



%% Calculate y(90) and y(-90)
y_90 = zeros(K);
y_90(1, 1) = 1;
for m = 1:Order
   y_90(m^2 + 1:(m+1)^2, m^2 + 1:(m+1)^2) = eval(sprintf('Ry90_%02d', m));
end

y90 = inv(y_90);

%% Calculate z
z = sym(zeros(K));
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

%% Calculate z(90) and z(-90)
z90 = subs(z, phi, pi/2);

z_90 = inv(z90);

%%
% Calculate z(90) and z(-90)
% z90old = zeros(K);
% [z90old(1, 1), diag, antidiag] = deal(1);
% for m = 1:Order
%     temp = zeros(2*m + 1);
%     n = size(temp,1);
%     
%     a = cos(m*pi/2);
%     b = sin(m*pi/2);
%     
%     [temp(1:n+1:end  ), diag    ] = deal([ a, diag,     a]);
%     [temp(n:n-1:end-1), antidiag] = deal([-b, antidiag, b]);
%     
%     z90old(m^2 + 1:(m+1)^2, m^2 + 1:(m+1)^2) = temp;
% end
% 
% z_90 = inv(z90);

%% Optimisation

yStart = (y_90\z90);
yEnd   = (z90\y_90);

%% zyx
zZ = subs(z, phi, Z);
zY = subs(z, phi, Y);
zX = subs(z, phi, X);

% %    Roll          Pitch                   Yaw
% B_ = y90*zX*y_90 * z_90*y_90*zY*y90*z90 *  zZ * B
% 
% %    Roll             Pitch                       Yaw
% B_ = y_90\zX*y_90 * (z90\y_90)*zY*(y_90\z90) *  zZ * B

%    Roll           Pitch             Yaw
B_ = y_90\zX*y_90 * yEnd*zY*yStart *  zZ * B

end