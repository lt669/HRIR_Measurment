function asxMultiWavEncode(degree, varargin)
%


%% asxWavEncode has changed!!!!! (save folder...)


% Version 1.0 Dec 2016
%
% Concatinated Ambisonic File Created From... 
%   f(order of ambisonics required, options, 
%     list of files to encode with angles e.g.: file1, phi, theta,
%                                               file2, phi, theta,
%                                               file3, phi, theta,)
%
% Options:       Angular unit: Radians    (rad)          - default
%                              Degrees    (deg) 
%
%               Normalisation: N3D        (n3d)          - default
%                              SN3D       (sn3d)
%                              Orthonomal (orthonormal)
%
% Script simplifies the process of encoding multiple mono audio files 
% into a single ambisonic file. The function automatically renders a 
% multi-channel ambisonic .wav file encoding input source audios 
% (varargin - file1, file2 etc.) in the spatial locations 
% (varargin - phi, theta etc.). Audio is encoded in degree (degree) with 
% options (varargin - options). One of three normalisation schemes can 
% be chosen and angles can be given in either radians or degrees. If a 
% multichannel / stereo file is input, the first / Left channel is 
% selected as the input. 
%
% NOTE: Options must be listed before audio files / angles!
%
% An ambisonic co-ordinate system is used where phi is azimuth; 0degrees 
% is straight ahead; positive is anticlockwise. Theta is elevation; 
% 0degrees is on the horizontal axis; positive is upwards.
%
% Written by Calum Armstrong, Department of Electronics, The University 
% of York.
%
% asxMultiWavEncode(1, 'deg', 'n3d',...
% 'C:\Users\ca718\Google Drive\PhD\Large_Data\Test_Files\Mono\48\front_48.wav', 4.19, 0,...
% 'C:\Users\ca718\Google Drive\PhD\Large_Data\Test_Files\Mono\48\back_48.wav' , 175.81, 0,...
% 'C:\Users\ca718\Google Drive\PhD\Large_Data\Test_Files\Mono\48\left_48.wav' , 90, 0,...
% 'C:\Users\ca718\Google Drive\PhD\Large_Data\Test_Files\Mono\48\right_48.wav', 270, 0,...
% 'C:\Users\ca718\Google Drive\PhD\Large_Data\Test_Files\Mono\48\up_48.wav'   , 90, 85.81,...
% 'C:\Users\ca718\Google Drive\PhD\Large_Data\Test_Files\Mono\48\down_48.wav' , 90, -85.81);
% 
% Updates ---
% 12 / 12 / 2016: Initial coding
% 13 / 12 / 2016: Finalisation of version 1.0

    % Set default normalisation
    norm = 'n3d'; 
    % Set default angular units
    units = 'rad';
    
    % Register and handle input options
    k = 1;
    first = 1;
    while k <= length(varargin) % For each input option given...
       switch varargin{k}
           case 'deg' % Angular units in degrees
               units = 'deg';
               k = k+1;
           case 'rad' % Angular units in radians
               units = 'rad';
               k = k+1;
           case 'n3d' % N3D normalisation
               norm = 'n3d';
               k = k+1;
           case 'sn3d' % SN3D normalisation
               norm = 'sn3d';
               k = k+1;
           case 'orthonormal' % Orthnormal normalisation
               norm = 'orthonormal';
               k = k+1;
           otherwise % List of files / angles
               % Encode file
               [matrix, FS] = asxWavEncode(varargin{k},...
                                                   varargin{k+1},...
                                                   varargin{k+2},...
                                                   degree,...
                                                   units,...
                                                   norm,...
                                                   'supOut');
               
               % Concatinate encoded matricies
               if first == 1
                   oldMatrix = matrix;
                   first = 0;
               else
                   oldMatrix = cat(1, oldMatrix, matrix);
               end
               
               k = k+3;
       end
    end

    % Generate filename
    if isa(degree, 'numeric')
        degree = num2str(degree);
    end
    filename = sprintf('Ambi_concat_degree%s_3d_%s_%d.wav',...
                        degree, norm, FS);
                    
    % Output .wav file from concatinated matrix
    audiowrite(filename, oldMatrix, FS);

end