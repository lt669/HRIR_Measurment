% Script will loop through a set of filenames stored in 'filenames',
% convole each with a White Noise burst in turn and play the result. It 
% has been used to test entire sets of HRTFs for locational accuracy.
%
% Written by Calum Armstrong, Department of Electronics, The University 
% of York.

%% Code

% For each file in 'filename'
for i = 1:50

    [HRTF, HRTFfs] = audioread(sprintf('48K_24bit_KEMAR_DFC/%s', filenames{i}));
    [File, Filefs] = audioread('WN_48.wav');

% Convolve
    outputL = conv(File(:, 1), HRTF(:, 1));
    outputR = conv(File(:, 1), HRTF(:, 2));
    output(:, :, i) = [outputL, outputR];

% Plot
%     figure;
%     plot(output);

% Display and play
    disp(filenames{i});
    soundsc(output(:, :, i), 44100);
    pause(0.6);

end