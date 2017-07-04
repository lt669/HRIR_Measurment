function [ output, FS ] = iasxA2B ( lS, HRTFs, phiDeg, thetaDeg, aFS, normTo1 )

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

%% Decode Virtual Loudspeaker Signals to Binaural
    
% Read in dummy HRTF to calculate size needed for output vector
    filename = dir(sprintf('%s/azi_*.wav', HRTFs));
    [HRTF, hrtfFS] = audioread(sprintf('%s/%s',...
                                         HRTFs, filename(1).name));
    
% Throw error if sampling frequencies of ambisonic file and HRTFs
% dont match
    if aFS ~= hrtfFS
        error(['Sampling frequency mis-match between ambisonic '...
               'file and HRTFs']);
    end
    
% Instantiate output vector
    output = zeros(max([...
                   length(lS(:, 1))+length(HRTF(:, 1))-1,...   
                   length(lS(:, 1)),...
                   length(HRTF(:, 1))]), 2);

% Loop through each loudspeaker signal, convolve with HRTF 
% of matching angle, and sum outputs
    for x = 1:length(lS(1, :))
       filename = dir(sprintf('%s/azi_%d_ele_%d*.wav',...
                                   HRTFs, phiDeg(x), thetaDeg(x)));
       if length(filename) > 1
           error('Multiple HRTF files for direction: azi %d ele %d', ...
                                                phiDeg(x), thetaDeg(x));
       end
       [HRTF, FS] = audioread(sprintf('%s/%s',...
                                         HRTFs, filename(1).name));
       output(:, 1) = output(:, 1) + ...
                                conv(lS(:, x),  HRTF(:, 1));   
       output(:, 2) = output(:, 2) + ...
                                conv(lS(:, x),  HRTF(:, 2));
    end

% Normalise output
    if normTo1
        output = output / max(abs(output(:)));
    end
end