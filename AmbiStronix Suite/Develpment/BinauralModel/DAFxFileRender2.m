%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RENDER DIRECT HRTF BINAURTAL FILES (Top Left 334 points)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RENDER DIRECT HRTF BINAURTAL FILES (lebedev 50 points)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% HRTFFolder = 'C:\Users\ca718\Google Drive\PhD\Large_Data\HRTF_Data_Sets\RIG\PinkCapsules\Calum\HRTFs\Calum_Center_trimmed';
% fileFolder = 'C:\Users\ca718\Google Drive\PhD\Large_Data\Test_Files\Mono\White_Noise\0-1s';
% saveFolder = 'C:\Users\ca718\Google Drive\PhD\Large_Data\Test_Files\Binaural\HRTF_Render_Calum_(LebDev50pt)';
% 
% [HRTF, ~] = audioread(sprintf('%s/azi_0_ele_0_FFC.wav', HRTFFolder));
% [File, Filefs] = audioread(sprintf('%s/WN_48.wav', fileFolder));
% output = zeros(length(File(:, 1))+length(HRTF(:, 1))-1, 2, 50);
% 
% [~, ~, LSphiDeg, LSthetaDeg] = iasxLDir('50ptLebedev');
% 
% for i = 1:50 
%         
%     phi = LSphiDeg(i);
%     theta = LSthetaDeg(i);
%         
%     [HRTF, ~] = audioread(sprintf('%s/azi_%d_ele_%d_FFC.wav', HRTFFolder, phi, theta));
% 
%     % Convolve
%     outputL = conv(File(:, 1), HRTF(:, 1));
%     outputR = conv(File(:, 1), HRTF(:, 2));
%     output(:,:, i) = [outputL, outputR];
% 
%     fprintf('CALCULATE: phi: %d; theta: %d\n', phi, theta);
%     
% end
% 
% outputNorm = output / max(abs(output(:)));
% 
% for i = 1:50
%     
%     phi = LSphiDeg(i);
%     theta = LSthetaDeg(i);
%     
%     filename = sprintf('%s/Bi_HRTF_%d_%d.wav',...
%                        saveFolder,...
%                        phi,...
%                        theta);
%                        
%     audiowrite(filename, outputNorm(:, :, i), 48000);
%     
%     fprintf('RENDER: %d of 50\n', i);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ENCODE AMBISONIC FILES (Top Left 334 points)(1, 3, 5 order) BAD NORMALISATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ENCODE AMBISONIC FILES (lebedev 50 points)(1, 3, 5 order)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%       file = 'C:\Users\ca718\Google Drive\PhD\Large_Data\Test_Files\Mono\White_Noise\0-1s\WN_48.wav';
% saveFolder = 'C:\Users\ca718\Google Drive\PhD\Large_Data\Test_Files\Ambisonic_Format\WhiteNoise(Lebedev50points)';
% 
% [~, ~, LSphiDeg, LSthetaDeg] = iasxLDir('50ptLebedev');
% 
% output = zeros(48000, 36, 250);
% 
% for order = 1:2:5
%    for i = 1:50
%        
%        phi = LSphiDeg(i);
%        theta = LSthetaDeg(i);
%        
%        [output(:, 1:(order+1)^2, (50*order-1)+i), fs] = asxWavEncode(file, phi, theta, order, saveFolder, 'deg', 'n3d', 'supOut', 'raw');
%                 
%        fprintf('Calculate degree: %d; phi: %d; theta: %d\n', order, phi, theta);
% 
%    end
% end
% 
% outputNorm = output / max(abs(output(:)));
% 
% [~, sourceFile, ~] = fileparts(sourceFile);
% 
% for order = 1:2:5
%    for i = 1:50
% 
%        phi = LSphiDeg(i);
%        theta = LSthetaDeg(i);
%     
%        filename = sprintf('%s/Ambi_%s_%d_%d_%d_%s.wav',...
%                           saveFolder,...
%                           sourceFile,...
%                           order,...
%                           phi,...
%                           theta,...
%                           'n3d'); % Construct filename
%                                
%        audiowrite(filename, outputNorm(:, 1:(order+1)^2, (50*order-1)+i), fs); % Save file
%          
%        fprintf('degree: %d; phi: %d; theta: %d\n', order, phi, theta);
% 
%    end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SINGLE FOCUS AMBISONIC DECODE (Top Left 334 points)(single layout - 1, 3, 5 order)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% [HRTF, ~] = audioread('C:/Users/ca718/Google Drive/PhD/Large_Data/HRTF_Data_Sets/RIG/PinkCapsules/Calum/HRTFs/Calum_Center_trimmed/azi_0_ele_0_FFC.wav');
% [File, ~] = audioread('C:/Users/ca718/Google Drive/PhD/Large_Data/Test_Files/Mono/White_Noise/0-1s/WN_48.wav');
% i = 1;
% output = zeros(length(File(:, 1))+length(HRTF(:, 1))-1, 2, 1002);
% 
%     config = 'octohedron'; %
%     decode = 'pinv'; %
%     saveFolder = 'C:\Users\ca718\Google Drive/PhD/Publications/!Conferences/AES 2017/AudioFiles/center';
%     
%     weighting = 'maxre'; %
%     normTo1 = 0; %
%     HRTFFolder = 'C:\Users\ca718\Google Drive/PhD/Large_Data/HRTF_Data_Sets/RIG/PinkCapsules/Calum/HRTFs/Calum_Center_trimmed'; %
%     dim = '3d'; %
%     norm = 'n3d'; %
%     
%     [LSphi, LStheta, LSphiDeg, LSthetaDeg] = iasxLDir(config);
%     
% for degree = 1:2:5
%     M = degree; %
%     K = (M + 1)^2; %
%     
%     [D, ~] = iasxD(LSphi, LStheta, decode, M, dim, norm);
%     D = iasxDWeighting(D, weighting, M);
%     
%     [ambisonic, aFS] = audioread(sprintf('C:\\Users\\ca718\\Google Drive//PhD//Large_Data//Test_Files//Ambisonic_Format//WhiteNoise//Ambi_WN_48_%d_0_90_n3d.wav', degree));
%     lS = (D * ambisonic(:, 1:K)')';
% 
%     [output(:,:, i), ~] = iasxA2B  (lS,...
%                                     HRTFFolder,...
%                                     LSphiDeg,...
%                                     LSthetaDeg,...
%                                     aFS,...
%                                     normTo1);
%                                 
%     i = i + 1;
%     
%     fprintf('CALCULATE: degree: %d; phi: 0; theta: 90\n', degree);
%     
%    for phi = 0:5:180
%       for theta = 0:10:80 
%         
%         [ambisonic, aFS] = audioread(sprintf('C:\\Users\\ca718\\Google Drive//PhD//Large_Data//Test_Files//Ambisonic_Format//WhiteNoise//Ambi_WN_48_%d_%d_%d_n3d.wav', degree, phi, theta));
%         lS = (D * ambisonic(:, 1:K)')';
%         
%         [output(:,:, i), ~] = iasxA2B  (lS,...
%                                         HRTFFolder,...
%                                         LSphiDeg,...
%                                         LSthetaDeg,...
%                                         aFS,...
%                                         normTo1);
%         
%         i = i + 1;
%         
%         fprintf('CALCULATE: degree: %d; phi: %d; theta: %d\n', degree, phi, theta);
%       end
%    end
% end
% 
% outputNorm = output / max(abs(output(:)));
% 
% for index = 1:1002
%     degree = (ceil((index)/334)-1) * 2 + 1;
%     phi = (ceil((rem(index-1, 334))/9)-1) * 5;
%     if phi == -5, phi = 0; end
%     theta = (rem(rem(index-1, 334)-1, 9)) * 10;
%     if theta == -10, theta = 90; end
%     
%     filename = sprintf('%s//Bi_test_%d_%d_%d.wav',...
%                        saveFolder,...
%                        degree,...
%                        phi,...
%                        theta);
%                        
%     audiowrite(filename, outputNorm(:, :, index), 48000);
%     
%     %fprintf('RENDER: %d of 1002\n', index);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SINGLE FOCUS AMBISONIC DECODE (lebedev 50 points)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     HRTFFolder = 'C:/Users/ca718/Google Drive/PhD/Large_Data/HRTF_Data_Sets/RIG/PinkCapsules/Calum/HRTFs/Calum_Center_trimmed'; %
%     ambiFolder = 'C:\Users\ca718\Google Drive\PhD\Large_Data\Test_Files\Ambisonic_Format\WhiteNoise(Lebedev50points)';
%           file = 'C:/Users/ca718/Google Drive/PhD/Large_Data/Test_Files/Mono/White_Noise/0-1s/WN_48.wav';
%     saveFolder = 'C:/Users/ca718/Google Drive/PhD/Publications/!Conferences/AES 2017/AudioFiles/center';
%     
%             M = 1; % degree
%        config = 'octohedron'; 
%        decode = 'pinv'; 
%     weighting = 'basic'; 
%     
%     dim = '3d'; %
%     norm = 'n3d'; %
%     normTo1 = 0; %
% 
%     [HRTF, ~] = audioread(sprintf('%s/azi_0_ele_0_FFC.wav', HRTFFolder));
%     [File, ~] = audioread(file);
%     output = zeros(length(File(:, 1))+length(HRTF(:, 1))-1, 2, 50);
% 
%     [LSphi, LStheta, LSphiDeg, LSthetaDeg] = iasxLDir(config);
%     [~, ~, AmbiphiDeg, AmbithetaDeg] = iasxLDir('50ptLebedev');
%     
%     [D, ~] = iasxD(LSphi, LStheta, decode, M, dim, norm);
%     D = iasxDWeighting(D, weighting, M);
%     
%    for i = 1:50 
%         
%         phi = AmbiphiDeg(i);
%         theta = AmbithetaDeg(i);
%        
%         [ambisonic, aFS] = audioread(sprintf('%s/Ambi_WN_48_%d_%d_%d_n3d.wav', ambiFolder, M, phi, theta));
%         lS = (D * ambisonic')';
%         
%         [output(:,:, i), ~] = iasxA2B  (lS,...
%                                         HRTFFolder,...
%                                         LSphiDeg,...
%                                         LSthetaDeg,...
%                                         aFS,...
%                                         normTo1);
%         
%         fprintf('CALCULATE: degree: %d; phi: %d; theta: %d\n', M, phi, theta);
%    end
% 
% outputNorm = output / max(abs(output(:)));
% 
% for i = 1:50
%     phi = AmbiphiDeg(i);
%     theta = AmbithetaDeg(i);
%     
%     filename = sprintf('%s//Bi_test_%d_%d_%d.wav',...
%                        saveFolder,...
%                        M,...
%                        phi,...
%                        theta);
%                        
%     audiowrite(filename, outputNorm(:, :, i), 48000);
%     
%     %fprintf('RENDER: %d of 50\n', index);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Multi-FOCUS AMBISONIC DECODE (lebedev 50 points)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    HRTFFolder = 'C:/Users/ca718/Google Drive/PhD/Large_Data/HRTF_Data_Sets/RIG/PinkCapsules/Calum/HRTFs/Calum_MF_ITDComp_trimmed'; %
    ambiFolder = 'C:\Users\ca718\Google Drive\PhD\Large_Data\Test_Files\Ambisonic_Format\WhiteNoise(Lebedev50points)';
          file = 'C:/Users/ca718/Google Drive/PhD/Large_Data/Test_Files/Mono/White_Noise/0-1s/WN_48.wav';
    saveFolder = 'C:\Users\ca718\Google Drive\PhD\Publications\!Conferences\AES 2017\AudioFiles\multifocus';
    
            M = 5; % degree
       config = '50ptLebedev'; 
       decode = 'pinv'; 
    weighting = 'basic'; 
    
    dim = '3d'; %
    norm = 'n3d'; %
    normTo1 = 0; %

    [HRTF, ~] = audioread(sprintf('%s/azi_0_ele_0_FFC.wav', HRTFFolder));
    [File, ~] = audioread(file);
    output = zeros(length(File(:, 1))+length(HRTF(:, 1))-1, 2, 50);
 
    [LSphi, LStheta, LSphiDeg, LSthetaDeg] = iasxLDir(config);
    [~, ~, AmbiphiDeg, AmbithetaDeg] = iasxLDir('50ptLebedev');
    
    rot = asxRotMatrix_zyx(M);
    
    [D, ~] = iasxD(LSphi, LStheta, decode, M, dim, norm);
    D = iasxDWeighting(D, weighting, M);
    
   for i = 1:50
        
        phi = AmbiphiDeg(i);
        theta = AmbithetaDeg(i);
       
        [ambisonic, aFS] = audioread(sprintf('%s/Ambi_WN_48_%d_%d_%d_n3d.wav', ambiFolder, M, phi, theta));
        ambisonicL = asx3DRotate(ambisonic', rot, -4.19*pi/180, 0, 0);
        lS(:, :, 1) = (D * ambisonicL)';
        ambisonicR = asx3DRotate(ambisonic', rot,  4.19*pi/180, 0, 0);
        lS(:, :, 2) = (D * ambisonicR)';
        
        [output(:,:, i), ~] = iasxA2BMF  (lS,...
                                          HRTFFolder,...
                                          LSphiDeg,...
                                          LSthetaDeg,...
                                          aFS,...
                                          normTo1);
        
        fprintf('CALCULATE: degree: %d; phi: %d; theta: %d\n', M, phi, theta);

   end

outputNorm = output / max(abs(output(:)));

for i = 1:50
    phi = AmbiphiDeg(i);
    theta = AmbithetaDeg(i);
    
    filename = sprintf('%s/Bi_test_%d_%d_%d.wav',...
                       saveFolder,...
                       M,...
                       phi,...
                       theta);
                       
    audiowrite(filename, outputNorm(:, :, i), 48000);
    
    %fprintf('RENDER: %d of 1002\n', index);
end
