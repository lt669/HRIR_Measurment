HRTF = 'C:/Users/ca718/Google Drive/PhD/Large_Data/Test_Files/Binaural/HRTF_Render';
compar = 'C:\Users\ca718\Google Drive\PhD\Publications\!Conferences\AES 2017\AudioFiles\center';
titleextra = 'FB';

[phiMesh,thetaMesh] = meshgrid(0:5:180,0:10:90);
phiMesh = phiMesh * pi / 180;
thetaMesh = thetaMesh * pi / 180;
x = cos(thetaMesh) .* cos(phiMesh); 
y = cos(thetaMesh) .* sin(phiMesh);
z = sin(thetaMesh);
                
[rL, rR, sL, sR, ILDD, ITDD] = deal(zeros(10,37,3));

% Filt = designfilt('lowpassfir', 'PassbandFrequency', 516.8,...
%          'StopbandFrequency', 742.9, 'PassbandRipple', 0.5, ...
%          'StopbandAttenuation', 75, 'SampleRate', 48000,...
%          'DesignMethod', 'kaiserwin');

% Filt = designfilt('highpassfir', 'PassbandFrequency', 4000,...
%          'StopbandFrequency', 3000, 'PassbandRipple', 0.5, ...
%          'StopbandAttenuation', 75, 'SampleRate', 48000,...
%          'DesignMethod', 'kaiserwin');     

Filt = 'N/A';

for method = 1:1 % decode method
    for d = 1 % order
        for p = 1:37 % azi
            for t = 1:9 % phi
                
                degree = ((d-1)*2) + 1;
                phi = (p-1) * 5;
                theta = (t-1) * 10;
                
                [rL(t,p,d),...
                 rR(t,p,d),...
                 sL(t,p,d),...
                 sR(t,p,d),...
                 ITDD(t,p,d),...
                 ILDD(t,p,d)] = DAFxFileComparison_func(HRTF, compar, degree, phi, theta, Filt);
                
                i = i + 1;
                
            end
        end
    end
end

for d = 1
    [rL(10,:,d),...
     rR(10,:,d),...
     sL(10,:,d),...
     sR(10,:,d),...
     ITDD(10,:,d),...
     ILDD(10,:,d)] = DAFxFileComparison_func(HRTF, compar, ((d-1)*2) + 1, 0, 90, Filt);
end

% norm_sL = 1 - abs(1 - sL);
% norm_sR = 1 - abs(1 - sL);

norm_sL = min(sL, 1/sL);
norm_sR = min(sR, 1/sR);

ILDD = min(ILDD, 1/ILDD);

errL = rL .* norm_sL;
errR = rR .* norm_sR;

%% SPEC. RMS ERROR / ITDD / ILDD 1ST ORDER
figure('name', sprintf('%s - %s',compar,titleextra), 'Position', [20 30 900 900]);
subplot('Position', [0.037 0.2 0.3 0.6]); specErr1L = surf(x,y,z,errL(:,:,1)); % Plot surface
%text(1.05,  0.69, 0, {'Spectral Error','(L)'},'HorizontalAlignment', 'left', 'FontSize',25, 'FontName', 'Times New Roman');
text(1.05,  0.84, 0, {'Spec. RMS Error','(L)'},'HorizontalAlignment', 'left', 'FontSize',25, 'FontName', 'Times New Roman');
formatDAFxGraph(specErr1L);

subplot('Position', [0.343 0.2 0.3 0.6]); specErr1R = surf(x,y,z,errR(:,:,1)); % Plot surface
%text(1.05,  0.69, 0, {'Spectral Error','(R)'},'HorizontalAlignment', 'left', 'FontSize',25, 'FontName', 'Times New Roman');
text(1.05,  0.84, 0, {'Spec. RMS Error','(R)'},'HorizontalAlignment', 'left', 'FontSize',25, 'FontName', 'Times New Roman');
formatDAFxGraph(specErr1R);

subplot('Position', [0.658 0.2 0.15 0.3]); ITDD1 = surf(x,y,z,ITDD(:,:,1)); % Plot surface
text(1.1, 0, 0, {'ITDD',''},'HorizontalAlignment', 'right','FontSize',25, 'FontName', 'Times New Roman'); 
formatDAFxGraph(ITDD1); 

subplot('Position', [0.812 0.2 0.15 0.3]); ILDD1 = surf(x,y,z,ILDD(:,:,1)); % Plot surface
text(1.1, 0, 0, {'ILDD',''},'HorizontalAlignment', 'right','FontSize',25, 'FontName', 'Times New Roman'); 
formatDAFxGraph(ILDD1);

colormap(bone);

%% ITDD / ILDD 3RD ORDER
% figure('name', sprintf('%s - %s',compar,titleextra), 'Position', [20 30 900 900]);
% subplot('Position', [0.097 0.1 0.4 0.8]); specErr1L = surf(x,y,z,ITDD(:,:,3)); % Plot surface
% %text(1.05,  0.52, 0, {'Spectral Error','(L)'},'HorizontalAlignment', 'left', 'FontSize',25, 'FontName', 'Times New Roman');
% %text(1.05,  0.65, 0, {'RMS Level Error','(L)'},'HorizontalAlignment', 'left', 'FontSize',25, 'FontName', 'Times New Roman');
% text(1.05,  0.22, 0, {'ITDD',''},'HorizontalAlignment', 'left', 'FontSize',25, 'FontName', 'Times New Roman');
% formatDAFxGraph(specErr1L);
% 
% subplot('Position', [0.503 0.1 0.4 0.8]); specErr1R = surf(x,y,z,ILDD(:,:,3)); % Plot surface
% %text(1.05,  0.52, 0, {'Spectral Error','(R)'},'HorizontalAlignment', 'left', 'FontSize',25, 'FontName', 'Times New Roman');
% %text(1.05,  0.65, 0, {'RMS Level Error','(R)'},'HorizontalAlignment', 'left', 'FontSize',25, 'FontName', 'Times New Roman');
% text(1.05,  0.22, 0, {'ILDD',''},'HorizontalAlignment', 'left', 'FontSize',25, 'FontName', 'Times New Roman');
% formatDAFxGraph(specErr1R);
% 
% colormap(bone);

%% SPECTRAL / RMS ERROR SEPERATE 1ST 3RD AND 5TH ORDER
figure('name', sprintf('%s - %s',compar,titleextra), 'Position', [20 100 800 600]);
subplot('Position', [0.039 0.543 0.14 0.373]); specErr1L = surf(x,y,z,rL(:,:,1)); % Plot surface
title({'Spectral','Err. (1)','(L)'},'HorizontalAlignment', 'left', 'Position', [0.7 1 0]); formatDAFxGraph(specErr1L);
subplot('Position', [0.181 0.543 0.14 0.373]); specErr1R = surf(x,-y,z,rR(:,:,1)); % Plot surface
title({'Spectral','Err. (1)','(R)'},'HorizontalAlignment', 'right', 'Position', [0.7 -1 0]); formatDAFxGraph(specErr1R);
subplot('Position', [0.359 0.543 0.14 0.373]); specErr3L = surf(x,y,z,rL(:,:,2)); % Plot surface
title({'Spectral','Err. (3)','(L)'},'HorizontalAlignment', 'left', 'Position', [0.7 1 0]); formatDAFxGraph(specErr3L);
subplot('Position', [0.501 0.543 0.14 0.373]); specErr3R = surf(x,-y,z,rR(:,:,2)); % Plot surface
title({'Spectral','Err. (3)','(R)'},'HorizontalAlignment', 'right', 'Position', [0.7 -1 0]); formatDAFxGraph(specErr3R);
subplot('Position', [0.679 0.543 0.14 0.373]); specErr5L = surf(x,y,z,rL(:,:,3)); % Plot surface
title({'Spectral','Err. (5)','(L)'},'HorizontalAlignment', 'left', 'Position', [0.7 1 0]); formatDAFxGraph(specErr5L);
subplot('Position', [0.821 0.543 0.14 0.373]); specErr5R = surf(x,-y,z,rR(:,:,3)); % Plot surface
title({'Spectral','Err. (5)','(R)'},'HorizontalAlignment', 'right', 'Position', [0.7 -1 0]); formatDAFxGraph(specErr5R);

subplot('Position', [0.039 0.085 0.14 0.373]); rmsErr1L = surf(x,y,z,norm_sL(:,:,1)); % Plot surface
title({'RMS','Err. (1)','(L)'},'HorizontalAlignment', 'left', 'Position', [0.7 1 0]); formatDAFxGraph(rmsErr1L);
subplot('Position', [0.181 0.085 0.14 0.373]); rmsErr1R = surf(x,-y,z,norm_sR(:,:,1)); % Plot surface
title({'RMS','Err. (1)','(R)'},'HorizontalAlignment', 'right', 'Position', [0.7 -1 0]); formatDAFxGraph(rmsErr1R);
subplot('Position', [0.359 0.085 0.14 0.373]); rmsErr3L = surf(x,y,z,norm_sL(:,:,2)); % Plot surface
title({'RMS','Err. (3)','(L)'},'HorizontalAlignment', 'left', 'Position', [0.7 1 0]); formatDAFxGraph(rmsErr3L);
subplot('Position', [0.501 0.085 0.14 0.373]); rmsErr3R = surf(x,-y,z,norm_sR(:,:,2)); % Plot surface
title({'RMS','Err. (3)','(R)'},'HorizontalAlignment', 'right', 'Position', [0.7 -1 0]); formatDAFxGraph(rmsErr3R);
subplot('Position', [0.679 0.085 0.14 0.373]); rmsErr5L = surf(x,y,z,norm_sL(:,:,3)); % Plot surface
title({'RMS','Err. (5)','(L)'},'HorizontalAlignment', 'left', 'Position', [0.7 1 0]); formatDAFxGraph(rmsErr5L);
subplot('Position', [0.821 0.085 0.14 0.373]); rmsErr5R = surf(x,-y,z,norm_sR(:,:,3)); % Plot surface
title({'RMS','Err. (5)','(R)'},'HorizontalAlignment', 'right', 'Position', [0.7 -1 0]); formatDAFxGraph(rmsErr5R);

colormap(bone);

%% SPEC. RMS Errors 1ST 3RD AND 5TH ORDER
% figure('name', compar, 'Position', [20 100 1500 500]);
% subplot('Position', [0.048 0.05 0.14 0.84]); totErr1L = surf(x,y,z,errL(:,:,1)); % Plot surface
% title('Total Spec. Err. (1) (L)'); formatDAFxGraph(totErr1L); 
% subplot('Position', [0.192 0.05 0.14 0.84]); totErr1R = surf(x,-y,z,errR(:,:,1)); % Plot surface
% title('Total Spec. Err. (1) (R)'); formatDAFxGraph(totErr1R);
% 
% subplot('Position', [0.358 0.05 0.14 0.84]); totErr3L = surf(x,y,z,errL(:,:,2)); % Plot surface
% title('Total Spec. Err. (3) (L)'); formatDAFxGraph(totErr3L);
% subplot('Position', [0.502 0.05 0.14 0.84]); totErr3R = surf(x,-y,z,errR(:,:,2)); % Plot surface
% title('Total Spec. Err. (3) (R)'); formatDAFxGraph(totErr3R);
% 
% subplot('Position', [0.668 0.05 0.14 0.84]); totErr5L = surf(x,y,z,errL(:,:,3)); % Plot surface
% title('Total Spec. Err. (5) (L)'); formatDAFxGraph(totErr5L);
% subplot('Position', [0.812 0.05 0.14 0.84]); totErr5R = surf(x,-y,z,errR(:,:,3)); % Plot surface
% title('Total Spec. Err. (5) (R)'); formatDAFxGraph(totErr5R);
% 
% colormap(hot);

%% ITDD / ILDD 1ST 3RD AND 5TH ORDER
% figure('name', compar, 'Position', [20 100 1500 500]);
% 
% subplot('Position', [0.048 0.05 0.14 0.84]); ITDD1 = surf(x,y,z,ITDD(:,:,1)); % Plot surface
% title('ITDD (1)'); formatDAFxGraph(ITDD1); 
% subplot('Position', [0.192 0.05 0.14 0.84]); ILDD1 = surf(x,-y,z,ILDD(:,:,1)); % Plot surface
% title('ILDD (1)'); formatDAFxGraph(ILDD1);
% 
% subplot('Position', [0.358 0.05 0.14 0.84]); ITDD3 = surf(x,y,z,ITDD(:,:,2)); % Plot surface
% title('ITDD (3)'); formatDAFxGraph(ITDD3);
% subplot('Position', [0.502 0.05 0.14 0.84]); ILDD3 = surf(x,-y,z,ILDD(:,:,2)); % Plot surface
% title('ILDD (3)'); formatDAFxGraph(ILDD3);
% 
% subplot('Position', [0.668 0.05 0.14 0.84]); ITDD5 = surf(x,y,z,ITDD(:,:,3)); % Plot surface
% title('ITDD (5)'); formatDAFxGraph(ITDD5);
% subplot('Position', [0.812 0.05 0.14 0.84]); ILDD5 = surf(x,-y,z,ILDD(:,:,3)); % Plot surface
% title('ILDD (5)'); formatDAFxGraph(ILDD5);
% 
% colormap(hot);

%% SPEC. RMS ERROR / ITDD / ILDD 1ST 3RD AND 5TH ORDER
% figure('name', sprintf('%s - %s',compar,titleextra), 'Position', [20 100 1500 500]);
% subplot('Position', [0.024 0.3 0.1 0.6]); totErr1L = surf(x,y,z,errL(:,:,1)); % Plot surface
% title({'Total Spectral','Err. (1)','(L)'},'HorizontalAlignment', 'left', 'Position', [0.7 1 0]); 
% formatDAFxGraph(totErr1L); 
% subplot('Position', [0.126 0.3 0.1 0.6]); totErr1R = surf(x,-y,z,errR(:,:,1)); % Plot surface
% title({'Total Spectral','Err. (1)','(R)'},'HorizontalAlignment', 'right', 'Position', [0.7 -1 0]); 
% formatDAFxGraph(totErr1R);
% 
% subplot('Position', [0.226 0.3 0.05 0.3]); ITDD1 = surf(x,y,z,ITDD(:,:,1)); % Plot surface
% title({'(1)','ITDD'},'HorizontalAlignment', 'left', 'Position', [-1.2 1 0]); 
% formatDAFxGraph(ITDD1); 
% subplot('Position', [0.277 0.3 0.05 0.3]); ILDD1 = surf(x,-y,z,ILDD(:,:,1)); % Plot surface
% title({'(1)','ILDD'},'HorizontalAlignment', 'right', 'Position', [-1.2 -1 0]); 
% formatDAFxGraph(ILDD1);
% 
% 
% subplot('Position', [0.349 0.3 0.1 0.6]); totErr3L = surf(x,y,z,errL(:,:,2)); % Plot surface
% title({'Total Spectral','Err. (3)','(L)'},'HorizontalAlignment', 'left', 'Position', [0.7 1 0]);
% formatDAFxGraph(totErr3L);
% subplot('Position', [0.451 0.3 0.1 0.6]); totErr3R = surf(x,-y,z,errR(:,:,2)); % Plot surface
% title({'Total Spectral','Err. (3)','(R)'},'HorizontalAlignment', 'right', 'Position', [0.7 -1 0]); 
% formatDAFxGraph(totErr3R);
% 
% subplot('Position', [0.551 0.3 0.05 0.3]); ITDD3 = surf(x,y,z,ITDD(:,:,2)); % Plot surface
% title({'(3)','ITDD'},'HorizontalAlignment', 'left', 'Position', [-1.2 1 0]);
% formatDAFxGraph(ITDD3);
% subplot('Position', [0.602 0.3 0.05 0.3]); ILDD3 = surf(x,-y,z,ILDD(:,:,2)); % Plot surface
% title({'(3)','ILDD'},'HorizontalAlignment', 'right', 'Position', [-1.2 -1 0]);
% formatDAFxGraph(ILDD3);
% 
% 
% subplot('Position', [0.674 0.3 0.1 0.6]); totErr5L = surf(x,y,z,errL(:,:,3)); % Plot surface
% title({'Total Spectral','Err. (5)','(L)'},'HorizontalAlignment', 'left', 'Position', [0.7 1 0]);
% formatDAFxGraph(totErr5L);
% subplot('Position', [0.776 0.3 0.1 0.6]); totErr5R = surf(x,-y,z,errR(:,:,3)); % Plot surface
% title({'Total Spectral','Err. (5)','(R)'},'HorizontalAlignment', 'right', 'Position', [0.7 -1 0]);
% formatDAFxGraph(totErr5R);
% 
% subplot('Position', [0.876 0.3 0.05 0.3]); ITDD5 = surf(x,y,z,ITDD(:,:,3)); % Plot surface
% title({'(5)','ITDD'},'HorizontalAlignment', 'left', 'Position', [-1.2 1 0]);
% formatDAFxGraph(ITDD5);
% subplot('Position', [0.927 0.3 0.05 0.3]); ILDD5 = surf(x,-y,z,ILDD(:,:,3)); % Plot surface
% title({'(5)','ILDD'},'HorizontalAlignment', 'right', 'Position', [-1.2 -1 0]);
% formatDAFxGraph(ILDD5);
% 
% colormap(hot);