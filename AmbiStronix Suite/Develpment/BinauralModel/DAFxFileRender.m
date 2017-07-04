%% ENCODE AMBISONIC FILES
%
% for degree = 1:2:5
%    for phi = 0:5:180
%       for theta = 0:10:80 
%         asxWavEncode(...
%         'WN_48.wav',...
%         phi, theta, degree, 'deg', 'n3d');
%         fprintf('degree: %d; phi: %d; theta: %d\n', degree, phi, theta);
%       end
%    end
% end
%
% Manual for 0, 90
% asxWavEncode('WN_48.wav',...
% 0, 90, 5, 'deg', 'n3d');

%% RENDER PLAIN HRTF BINAURTAL FILES
%
% [HRTF, ~] = audioread('../HRTFs/48K_24bit_KEMAR_DFC/azi_0_ele_0_DFC.wav');
% [File, Filefs] = audioread('WN_48.wav');
% 
% i = 1;
% output = zeros(length(File(:, 1))+length(HRTF(:, 1))-1, 2, 334);
%            
% for phi = 0:5:180
%   for theta = 0:10:80 
%     [HRTF, ~] = audioread(sprintf('../HRTFs/48K_24bit_KEMAR_DFC/azi_%d_ele_%d_DFC.wav', phi, theta));
% 
%     % Convolve
%     outputL = conv(File(:, 1), HRTF(:, 1));
%     outputR = conv(File(:, 1), HRTF(:, 2));
%     output(:,:, i) = [outputL, outputR];
%     i = i + 1;
% 
%     fprintf('CALCULATE: phi: %d; theta: %d\n', phi, theta);
%   end
% end
% 
% [HRTF, HRTFfs] = audioread('../HRTFs/48K_24bit_KEMAR_DFC/azi_0_ele_90_DFC.wav'); 
% outputL = conv(File(:, 1), HRTF(:, 1));
% outputR = conv(File(:, 1), HRTF(:, 2));
% output(:,:, 334) = [outputL, outputR];
% i = i + 1;
% fprintf('CALCULATE: phi: 0; theta: 90\n'); 
% 
% 
% outputNorm = output / max(abs(output(:)));
% 
% 
% for index = 1:333
%     phi = (ceil((index)/9)-1) * 5; 
%     theta = (rem(index-1, 9)) * 10;
%     
%     filename = sprintf('HRTF_Render/Bi_HRTF_%d_%d.wav',...
%                        phi,...
%                        theta);
%                        
%     audiowrite(filename, outputNorm(:, :, index), 48000);
%     
%     fprintf('RENDER: %d of 334\n', index);
% end
% 
% filename = sprintf('HRTF_Render/Bi_HRTF_0_90.wav');   
% audiowrite(filename, outputNorm(:, :, 334), 48000);
% fprintf('RENDER: %d of 334\n', 334);


%% SINGLE FOCUS AMBISONIC DECODE

[HRTF, ~] = audioread('C:\Users\ca718\Google Drive\PhD\Large_Data\HRTF_Data_Sets\RIG\PinkCapsules\Calum\HRTFs\Calum_MF_ITDComp_trimmed\azi_0_ele_0_FFC.wav');
[File, ~] = audioread('C:\Users\ca718\Google Drive\PhD\Large_Data\Test_Files\Mono\White_Noise\0-1s\WN_48.wav');
i = 1;
output = zeros(length(File(:, 1))+length(HRTF(:, 1))-1, 2, 1002);

    config = 'octohedron'; %
    decode = 'pinv'; %
    saveFolder = 'C:\Users\ca718\Google Drive\PhD\Publications\!Conferences\AES 2017\AudioFiles\multifocus';
    
    weighting = 'maxre'; %
    normTo1 = 0; %
    HRTFFolder = 'C:\Users\ca718\Google Drive\PhD\Large_Data\HRTF_Data_Sets\RIG\PinkCapsules\Calum\HRTFs\Calum_MF_ITDComp_trimmed'; %
    dim = '3d'; %
    norm = 'n3d'; %
    
    [LSphi, LStheta, LSphiDeg, LSthetaDeg] = iasxLDir(config);
    
for degree = 1:2:5
    M = degree; %
    K = (M + 1)^2; %
    
    [D, ~] = iasxD(LSphi, LStheta, decode, M, dim, norm);
    D = iasxDWeighting(D, weighting, M);
    
    [ambisonic, aFS] = audioread(sprintf('C:\\Users\\ca718\\Google Drive\\PhD\\Large_Data\\Test_Files\\Ambisonic_Format\\WhiteNoise\\Ambi_WN_48_%d_0_90_n3d.wav', degree));
    lS = (D * ambisonic(:, 1:K)')';

    [output(:,:, i), ~] = iasxA2B  (lS,...
                                    HRTFFolder,...
                                    LSphiDeg,...
                                    LSthetaDeg,...
                                    aFS,...
                                    normTo1);
                                
    i = i + 1;
    
    fprintf('CALCULATE: degree: %d; phi: 0; theta: 90\n', degree);
    
   for phi = 0:5:180
      for theta = 0:10:80 
        
        [ambisonic, aFS] = audioread(sprintf('C:\\Users\\ca718\\Google Drive\\PhD\\Large_Data\\Test_Files\\Ambisonic_Format\\WhiteNoise\\Ambi_WN_48_%d_%d_%d_n3d.wav', degree, phi, theta));
        lS = (D * ambisonic(:, 1:K)')';
        
        [output(:,:, i), ~] = iasxA2B  (lS,...
                                        HRTFFolder,...
                                        LSphiDeg,...
                                        LSthetaDeg,...
                                        aFS,...
                                        normTo1);
        
        i = i + 1;
        
        %fprintf('CALCULATE: degree: %d; phi: %d; theta: %d\n', degree, phi, theta);
      end
   end
end

outputNorm = output / max(abs(output(:)));

for index = 1:1002
    degree = (ceil((index)/334)-1) * 2 + 1;
    phi = (ceil((rem(index-1, 334))/9)-1) * 5;
    if phi == -5, phi = 0; end
    theta = (rem(rem(index-1, 334)-1, 9)) * 10;
    if theta == -10, theta = 90; end
    
    filename = sprintf('%s\\Bi_test_%d_%d_%d.wav',...
                       saveFolder,...
                       degree,...
                       phi,...
                       theta);
                       
    audiowrite(filename, outputNorm(:, :, index), 48000);
    
    %fprintf('RENDER: %d of 1002\n', index);
end

%% NEED TO ROTATE SOUNDFIELD FOR LEFT / RITE EAR
%% Multi-FOCUS AMBISONIC DECODE

[HRTF, ~] = audioread('C:\Users\ca718\Google Drive\PhD\Large_Data\HRTF_Data_Sets\RIG\PinkCapsules\Calum\HRTFs\Calum_MF_ITDComp_trimmed\azi_0_ele_0_FFC.wav');
[File, ~] = audioread('C:\Users\ca718\Google Drive\PhD\Large_Data\Test_Files\Mono\White_Noise\0-1s\WN_48.wav');
i = 1;
output = zeros(length(File(:, 1))+length(HRTF(:, 1))-1, 2, 1002);

    config = 'octohedron'; %
    decode = 'pinv'; %
    saveFolder = 'C:\Users\ca718\Google Drive\PhD\Publications\!Conferences\AES 2017\AudioFiles\multifocus';
    
    weighting = 'maxre'; %
    normTo1 = 0; %
    HRTFFolder = 'C:\Users\ca718\Google Drive\PhD\Large_Data\HRTF_Data_Sets\RIG\PinkCapsules\Calum\HRTFs\Calum_MF_ITDComp_trimmed'; %
    dim = '3d'; %
    norm = 'n3d'; %
    
    [LSphi, LStheta, LSphiDeg, LSthetaDeg] = iasxLDir(config);
    
for degree = 1:2:5
    M = degree; %
    K = (M + 1)^2; %
    
    [D, ~] = iasxD(LSphi, LStheta, decode, M, dim, norm);
    D = iasxDWeighting(D, weighting, M);
    
    [ambisonic, aFS] = audioread(sprintf('C:\\Users\\ca718\\Google Drive\\PhD\\Large_Data\\Test_Files\\Ambisonic_Format\\WhiteNoise\\Ambi_WN_48_%d_0_90_n3d.wav', degree));
    lS = (D * ambisonic(:, 1:K)')';

    [output(:,:, i), ~] = iasxA2B  (lS,...
                                    HRTFFolder,...
                                    LSphiDeg,...
                                    LSthetaDeg,...
                                    aFS,...
                                    normTo1);
                                
    i = i + 1;
    
    fprintf('CALCULATE: degree: %d; phi: 0; theta: 90\n', degree);
    
   for phi = 0:5:180
      for theta = 0:10:80 
        
        [ambisonic, aFS] = audioread(sprintf('C:\\Users\\ca718\\Google Drive\\PhD\\Large_Data\\Test_Files\\Ambisonic_Format\\WhiteNoise\\Ambi_WN_48_%d_%d_%d_n3d.wav', degree, phi, theta));
        lS = (D * ambisonic(:, 1:K)')';
        
        [output(:,:, i), ~] = iasxA2B  (lS,...
                                        HRTFFolder,...
                                        LSphiDeg,...
                                        LSthetaDeg,...
                                        aFS,...
                                        normTo1);
        
        i = i + 1;
        
        %fprintf('CALCULATE: degree: %d; phi: %d; theta: %d\n', degree, phi, theta);
      end
   end
end

outputNorm = output / max(abs(output(:)));

for index = 1:1002
    degree = (ceil((index)/334)-1) * 2 + 1;
    phi = (ceil((rem(index-1, 334))/9)-1) * 5;
    if phi == -5, phi = 0; end
    theta = (rem(rem(index-1, 334)-1, 9)) * 10;
    if theta == -10, theta = 90; end
    
    filename = sprintf('%s\\Bi_test_%d_%d_%d.wav',...
                       saveFolder,...
                       degree,...
                       phi,...
                       theta);
                       
    audiowrite(filename, outputNorm(:, :, index), 48000);
    
    %fprintf('RENDER: %d of 1002\n', index);
end