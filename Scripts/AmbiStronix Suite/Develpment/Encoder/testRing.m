% Script outputs an Ambisonic Format test file which localises a White
% Noise burst at 10 degree increments around the horizontal plane. An
% example decode can be computed as follows
%
% asxDecoder('Ambi_ring_test_1_sn3d.wav', '48K_24bit_KEMAR_DFC',...
%            'example', 'sn3d' );
%
% Written by Calum Armstrong, Department of Electronics, The University 
% of York.

%% Code

degree = 1;
norm = 'sn3d';

% Loop through azimuth directions and encode a White Noise burst at each
% point
first = 1;
for az = 0:10:360
% Encode
    [matrix, FS] = asxWavEncode('WN_48.wav',...
                                az,...
                                0,...
                                degree,...
                                'deg',...
                                norm,...
                                'supOut',...
                                'raw');

% Concatinate encoded matricies
    if first == 1
        concMatrix = matrix;
        first = 0;
    else
        concMatrix = cat(1, concMatrix, matrix);
    end
end

% Normalize output
concMatrix = concMatrix / max((abs(concMatrix(:))));

% Output .wav file 
audiowrite(sprintf('Ambi_ring_test_%d_%s_48.wav', degree, norm),...
           concMatrix, FS);
