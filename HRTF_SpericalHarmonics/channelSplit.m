

left = '/Users/Lewis 1/Documents/MATLAB/HRIR_MEASURMENT/Audio/Analysis/HRIR_SphericalHarmonic/KU_100/48/SP_Left.wav';
right = '/Users/Lewis 1/Documents/MATLAB/HRIR_MEASURMENT/Audio/Analysis/HRIR_SphericalHarmonic/KU_100/48/SP_Right.wav';

output = {                     'W',... 
                           'Y',  'Z',  'X',... 
                       'V',  'T',  'R',  'S', 'U',... 
                   'Q',  'O',  'M',  'K',  'L',  'N',  'P',... 
              'a4', 'b4', 'c4', 'd4', 'e4', 'f4', 'g4', 'h4', 'i4',... 
          'a5', 'b5', 'c5', 'd5', 'e5', 'f5', 'g5', 'h5', 'i5', 'j5', 'k5'};



x = audioread(right);      
mkdir('/Users/Lewis 1/Documents/MATLAB/HRIR_MEASURMENT/Audio/Analysis/HRIR_SphericalHarmonic/KU_100/48/split/left');
mkdir('/Users/Lewis 1/Documents/MATLAB/HRIR_MEASURMENT/Audio/Analysis/HRIR_SphericalHarmonic/KU_100/48/split/right');
for k=1:36
    y = x(:,k);
    name = strcat('/Users/Lewis 1/Documents/MATLAB/HRIR_MEASURMENT/Audio/Analysis/HRIR_SphericalHarmonic/KU_100/48/split/right/',char(output(k)),'.wav');
    audiowrite(name,y,48);
end