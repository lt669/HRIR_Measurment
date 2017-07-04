function P = SoundfieldSnapshot(L, r, S, res, area, T, fs, PlotTitle, clim, zoom)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function plots a normalized horizontal sound field with  wave
% interference from N sources. It is useful for plotting Ambisonic or
% Wavefield Synthesis type soundfields. Uses time-domain superposition.
%
% Inputs:
% L: Source Angles, [1 x N] in degrees
% r: Source radius from centre
% S: Loudspeaker signals (sine waves or pulses are best for analysis)
% fs: Sample rate of loudspeaker signals
% res: Plot resolution in meters
% area: Listening area dimensions [Length Breadth]
% T: Snapshot time in seconds
% PlotTitle: Character string for title of plot
% Plot magnitude limits
%
% Example:
% Source Data
% close all; clear all; clc;
% fs = 44100;
% r = 20;
% t = 0:1/fs: r*2/343;
% s = sin(2*pi*1000*t);
% res = 0.01;
% plotlim = [1 1];
% clim = [-1 1]
% T = r*2/343;
% zoom = 1
% figure(2)
% P = SoundfieldSnapshot(45, r, S', res, plotlim, T, fs, 'Plane wave, 1000Hz', clim, zoom);
% Gavin Kearney, Trinity College Dublin 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 343;

L = L*pi/180;

x_coord = r*cos(L);
y_coord = r*sin(L);


l = 1;

% xval = -area(1)/2:res:area(1)/2;
% yval = -area(2)/2:res:area(2)/2;

xval = -(area(1)*zoom):res:(area(1)*zoom);
yval = -(area(2)*zoom):res:(area(2)*zoom);

for y = yval
    m = 1;
    for x = -xval
        if sqrt(x.^2 + y.^2)<0.99*r
            for i = 1:length(L)
                dxyS(i) = sqrt((x-x_coord(i))^2 + (y-y_coord(i))^2);
                t1(i) = dxyS(i)/c;
                DelayPad1 = round(t1(i)*fs);
                DelayPad2 = round((T-t1(i))*fs);
                A(i) = (1/dxyS(i));
                Padded = [zeros(1,DelayPad1), S(i,1:(length(S)-DelayPad1)), zeros(1,DelayPad2)]; % Recorded Impulses
                Pt(i) = A(i)*Padded(round(T*fs));
            end
            P(l,m) = sum(Pt);
        else
            P(l,m) = 1; % This will plot as white
        end
        m = m+1;
    end
    l = l + 1
end

% Plot

%clim = ([-1 1]);
%imagesc(xval,yval, flipud(fliplr(P./max(abs(P(:))'))), clim);
imagesc(xval,yval, flipud(fliplr(P)'), clim);
colormap(gray);
colorbar
title(PlotTitle);
xlabel('X-offset (m)');
ylabel('Y-offset (m)');


