[HRTF, HRTFfs] = audioread('C:\Users\ca718\Google Drive\PhD\Large_Data\HRTF_Data_Sets\RIG\PinkCapsules\Calum\HRTFs\Calum_Centre_trimmed\azi_0_ele_90_FFC.wav');
[File, Filefs] = audioread('front_48.wav');

% Convolve
    outputL = conv(File(:, 1), HRTF(:, 1));
    outputR = conv(File(:, 1), HRTF(:, 2));
    output(:, :) = [outputL, outputR];
% Display and play
    soundsc(output, 48000);