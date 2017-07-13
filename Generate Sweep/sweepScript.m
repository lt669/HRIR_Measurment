% Generate sweep script
clc;
% path = 'Audio/Sweeps';
path = '/Users/Lewis 1/Documents/MATLAB/HRIR_MEASUREMENT_V2/Audio/GlobalAudio/Sweeps';
[x,invX] = generatesweep(path,20,22050,3,44100,0);
