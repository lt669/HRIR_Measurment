function [ output, FS ] = asxA2B ( lS, HRTFs, phiDeg, thetaDeg, aFS )

% Version 1.0 Feb 2017
%
% [ Binaural Signal,
%   Sampling Frequency ] From f( Loudspeaker Signals,
%                                Folder of HRTFs,
%                                Loudspeaker Azimuths,
%                                Loudspeaker Elevations,
%                                Loudspeaker signal sampling freq. )
%
% Script converts a set of loudspeaker signals 'lS' designed for 
% loudspeakers at locations ('phiDeg', 'thetaDeg') (in degrees) to a 
% binaural signal based on a set of HRTFs in the folder 'HRTFs'.
%
% Written by Calum Armstrong, Department of Electronics, The University 
% of York.
%
% Updates ---
% 15 / 02 / 2017: Initial coding
%%
    
% Instantiate output vector
    output = zeros(1257404, 1);
    figure;
    
% Loop through each loudspeaker signal, convolve with HRTF 
% of matching angle, and sum outputs
    for x = 22:29
       filename = dir(sprintf('%s/azi_%d_ele_%d*.wav',...
                                   HRTFs, phiDeg(x), thetaDeg(x)));

       [HTRF, FS] = audioread(sprintf('%s/%s',...
                                         HRTFs, filename(1).name));
       
       blah = conv(lS(:, x), HTRF(:, 1));
       
        
       subplot(4, 1, 1); plot(lS(:, x)); ylim([-2 2]);
       subplot(4, 1, 2); plot(HTRF(:, 1)); 
       subplot(4, 1, 3); plot(blah); ylim([-2 2]);
       
       output = output + blah;
       
       subplot(4, 1, 4); plot(output); ylim([-2 2]);
       
       %linkaxes(a, 'xy');
    end
    
end