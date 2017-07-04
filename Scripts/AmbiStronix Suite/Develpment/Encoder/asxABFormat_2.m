%%
%clear all; 
%clc; 

%% Read A-Format file begin
filename='speechmix_Aformat'; %A-Format file 
[sig, fs] = audioread([filename '.wav']);
% Read A-Format file end

%% Generate W_, X_, Y_, Z_ begin
% Note the '_' denotes the uncompensated form of the B-Format signals

% The configuration of the 4 Mics of Twirling720 recoding device
% signal -> (azimuth angle (+'ve clockwise), elevation angle)
% sig(:,1) 90º ,  35º Right Up
% sig(:,2) 180º, -35º Back Down
% sig(:,3) 0º  , -35º Front Down
% sig(:,4) 270º,  35º Left Up

% (Using horizontal / anticlockwise ordering of mic capsules - 3, 4, 2, 1)
W_ =   sig(:,3) + sig(:,4) + sig(:,2) + sig(:,1);
X_ =   sqrt(2) * (sig(:,3) - sig(:,2));
Y_ =   sqrt(2) * (sig(:,4) - sig(:,1));
Z_ = - sig(:,3) + sig(:,4) - sig(:,2) + sig(:,1);

% generate W_, X_, Y_, Z_ end

%% compansation filter processing to W_,X_,Y_,Z_ begin

% SPATIAL COMPENSATION
%                  (r = Array radius)
%                  spatialCompensation(r,      fs)
[filtW, filtXYZ] = spatialCompensation(0.0175, fs);
W = conv(W_, filtW);
X = conv(X_, filtXYZ);
Y = conv(Y_, filtXYZ);
Z = conv(Z_, filtXYZ);

% NOISE COMPENSATION
filtNoise = noiseCompensation(fs);
W2 = filter(filtNoise, W);
X2 = filter(filtNoise, X);
Y2 = filter(filtNoise, Y);
Z2 = filter(filtNoise, Z);

% % compansation filter processing to W_,X_,Y_,Z_ end

%% write data begin

BFormat_ = [W_ / sqrt(2), X_, Y_, Z_];
audiowrite( sprintf('%s_B_Format_.wav',filename),BFormat_, fs);

BFormat = [W / sqrt(2), X, Y, Z];
audiowrite( sprintf('%s_B_Format.wav',filename),BFormat, fs);

BFormatNoiseSupress = [W2 / sqrt(2), X2, Y2, Z2];
audiowrite( sprintf('%s_B_Format_NoiseSupress.wav',filename), BFormatNoiseSupress, fs);

sn3dHOA_ = [W_, Y_, Z_, X_];
audiowrite( sprintf('%s_sn3dHOA_.wav',filename), sn3dHOA_, fs);

sn3dHOA = [W, Y, Z, X];
audiowrite( sprintf('%s_sn3dHOA.wav',filename), sn3dHOA, fs);

sn3dHOANoiseSupress = [W2, Y2, Z2, X2];
audiowrite( sprintf('%s_sn3dHOA_NoiseSupress.wav',filename), sn3dHOANoiseSupress, fs);

% write data end 

