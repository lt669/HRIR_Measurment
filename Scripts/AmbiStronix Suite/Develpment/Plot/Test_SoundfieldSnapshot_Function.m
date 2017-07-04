% Test script for SoundfieldSnapshot function. Plots stereo interference.

close all; clear all; clc;
fs = 44100;
r = 1.5; % Loudspeaker array radius
t = 0:1/fs: r*2/343; % Time interval

s = sin(2*pi*440*t); % Test source at 440Hz

% Uncomment for pulse
%s = zeros(size(t));
%s(1) = 1;

res = 0.005; % Resolution
plotlim = [1 1];
T = r*1/343; % Time of snapshot (note its done by multiples of array radius)
figure()
s = s*10^(1/r); % Normalize Data
zoom = 1.5; % Sets a plot zoom (higher number = zoom out)
clim = [-2 2]
P = SoundfieldSnapshot([30 -30], r, [s; s], res, plotlim, T, fs, 'Stereo 440Hz', clim, zoom);