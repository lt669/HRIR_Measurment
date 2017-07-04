% Script will output a binaural file which localises a White
% Noise burst at 10 degree increments around the horizontal plane by 
% direct convolution with KEMAR HRTFs taken from the SADIE online 
% database. This provides a 'best case' rendering to which
% Ambisonic -> binaural decoders should be measured against.
%
% Written by Calum Armstrong, Department of Electronics, The University 
% of York.

%% Code

% Loop through azimuth directions and convolve a White Noise burst with 
% the relevant HRTF at each point
first = 1;
for i = 0:10:350

    [HRTF, HRTFfs] = audioread(sprintf(...
                  'HRTFs/44K_16bit_KEMAR_DFC/azi_%d_ele_0_DFC.wav', i));
    [File, Filefs] = audioread('WN_48.wav');

% Convolve
    outputL = conv(File(:, 1), HRTF(:, 1));
    outputR = conv(File(:, 1), HRTF(:, 2));
    matrix = [outputL, outputR];

% Concatinate outputs
   if first == 1
       concMatrix = matrix;
       first = 0;
   else
       concMatrix = cat(1, concMatrix, matrix);
   end
   
% Play output
%    soundsc(output(:, :, i), 44100);
%    pause(0.6);

end

% Normalize output
concMatrix = concMatrix / max((abs(concMatrix(:))));

% Output .wav file 
audiowrite('HRTF_test.wav', concMatrix, HRTFfs);