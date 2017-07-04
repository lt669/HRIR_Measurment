%%
clear all; 
clc; 

r = 0.075;  % Array radius
c = 343;    % Speed of sound
Fs = 44100; % Sample rate
Fstep = 10; % Frequency stepsize
%

%% Read A-Format file begin
filename='speechmix'; %A-Format file 
[sig, fs] = audioread([filename '.wav']);
% Read A-Format file end

%% Generate W_45, X_45, Y_45, Z_45 begin
% Note the '_' denotes the uncompensated form of the B-Format signals
% Note the '45' denotes a 45º (clockwise) horizontal offset of the B-Format Ambisonic Channels

% The configuration of the 4 Mics of Twirling720 recoding device
% signal -> (azimuth angle (+'ve clockwise), elevation angle)
% sig(:,1) 90º ,  35º Right Up
% sig(:,2) 180º, -35º Back Down
% sig(:,3) 0º  , -35º Front Down
% sig(:,4) 270º,  35º Left Up

% (Using horizontal / anticlockwise ordering of mic capsules - 3, 4, 2, 1)
W_45 =   sig(:,3) + sig(:,4) + sig(:,2) + sig(:,1); % (= W_)
X_45 =   sig(:,3) + sig(:,4) - sig(:,2) - sig(:,1);
Y_45 = - sig(:,3) + sig(:,4) + sig(:,2) - sig(:,1);
Z_45 = - sig(:,3) + sig(:,4) - sig(:,2) + sig(:,1); % (= Z_)

% generate W_45, X_45, Y_45, Z_45 end

%% Rotate / Re-align W_45, X_45, Y_45, Z_45 begin
% Define rotation matrix, R, to rotate soundfield 45 clockwise to re-align
% with B-Format Ambisonic Channels
R = [ 1,          0,           0,  0;... % W
      0,  2^(1/2)/2,  -2^(1/2)/2,  0;... % X
      0,  2^(1/2)/2,   2^(1/2)/2,  0;... % Y
      0,          0,           0,  1];   % Z

% Build B-Format signal matrix (ROWS represent Ambiasonic Channels -
% required for rotation matrix multiplication)
BFormat_45 = [W_45, X_45, Y_45, Z_45]' ;

% Rotate
BFormat_ = R * BFormat_45;

% Rotate / Re-align W_45, X_45, Y_45, Z_45 end

%% generate W1,X1,Y1,Z1 compansation filters begin
% i = 1;
% 
% for f = Fstep:Fstep:Fs/2 
%     w = 2*pi*f;
%    
%     % Compute W compensation response
% 
%     Num = 1 + (j*w*r)/c - (1/3)*(((w*r)/c)*((w*r)/c));
%     Den = 1 + (1/3)*(j*w*r)/c;
%     FW(i) = Num/Den;
%     
%     % Compute XYZ compensation response
%     
%     Num = sqrt(6)*(1 + (1/3)*(j*w*r)/c - (1/3)*(((w*r)/c)*((w*r)/c)));
%     Den = 1 + (1/3)*(j*w*r)/c;
%     FXYZ(i) = Num/Den;
%         
%     i = i+1;
% end
% 
% f = linspace(0, Fs/2, (Fs/2)/Fstep);
% % generate W1,X1,Y1,Z1 compansation filters end

%% compansation filter processing to W1,X1,Y1,Z1 begin
% 
% % compansation filter processing to W1,X1,Y1,Z1 end

%% write B-Format begin

audiowrite( sprintf('%s_B_Format.wav',filename),BFormat_', fs);

n3dHOA_ = [BFormat_(1, :) * sqrt(2); BFormat_(3, :); BFormat_(4, :); BFormat_(2, :)];
audiowrite( sprintf('%s_n3dHOA.wav',filename), n3dHOA_', fs);

% write W1,X1,Y1,Z1 separately end 