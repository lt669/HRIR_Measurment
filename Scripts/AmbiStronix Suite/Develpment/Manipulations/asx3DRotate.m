% Ambisonic SoundField Rotation

function B_ = asx3DRotate(B, R, Z, Y, X)

% B - Ambisonic Format signal representing soundfield to be rotated. ROWS
% should correspond to Ambisonic Channels in ACN format.
% R - Sysbolic Rotation Matrix composed by asxRotMatrix_....m
% X - Angle of rotation about the x axis
% Y - Angle of rotation about the y axis
% Y - Angle of rotation about the z axis
%
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

%% zyx

R = double(subs(R));

B_ = R * B;

end