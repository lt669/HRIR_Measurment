HRTF = 'C:\Users\ca718\Google Drive\PhD\Large_Data\Test_Files\Binaural\HRTF_Render_Calum_(LebDev50pt)';
compar = 'C:\Users\ca718\Google Drive\PhD\Publications\!Conferences\AES 2017\AudioFiles\multifocus';

[phirad, thetarad, phiDeg, thetaDeg] = iasxLDir('50ptLebedev');

x = cos(thetarad) .* cos(phirad); 
y = cos(thetarad) .* sin(phirad);
z = sin(thetarad);
                
[rL, rR, sL, sR, ILDD, ITDD] = deal(zeros(50,5));

% Filt = designfilt('lowpassfir', 'PassbandFrequency', 516.8,...
%          'StopbandFrequency', 742.9, 'PassbandRipple', 0.5, ...
%          'StopbandAttenuation', 75, 'SampleRate', 48000,...
%          'DesignMethod', 'kaiserwin');

% Filt = designfilt('highpassfir', 'PassbandFrequency', 4000,...
%          'StopbandFrequency', 3000, 'PassbandRipple', 0.5, ...
%          'StopbandAttenuation', 75, 'SampleRate', 48000,...
%          'DesignMethod', 'kaiserwin');     

Filt = 'N/A';

for d = 1:2:5 % order
    for i = 1:50

        phi = phiDeg(i);
        theta = thetaDeg(i);

        [rL(i,d),...
         rR(i,d),...
         sL(i,d),...
         sR(i,d),...
         ITDD(i,d),...
         ILDD(i,d)] = DAFxFileComparison_func(HRTF, compar, d, phi, theta, Filt);


    end
end

norm_sL = min(sL, 1./sL);
norm_sR = min(sR, 1./sR);
ILDD = min(ILDD, 1./ILDD);

% avgrL = mean(rL, 1);
% avgrR = mean(rR, 1);
% avgsL = mean(norm_sL, 1);
% avgsR = mean(norm_sR, 1);
% avgITDD = mean(ITDD, 1);
% avgILDD = mean(ILDD, 1);

% crL = rL;
% crR = rR;
% csL = norm_sL;
% csR = norm_sR;
% cITDD = ITDD;
% cILDD = ILDD;
% cavgrL = avgrL;
% cavgrR = avgrR;
% cavgsL = avgsL;
% cavgsR = avgsR;
% cavgITDD = avgITDD;
% cavgILDD = avgILDD;

% diffrL = rL - crL;
% diffrR = rR - crR;
% diffsL = norm_sL - csL;
% diffsR = norm_sR - csR;
% diffITDD = ITDD - cITDD;
% diffILDD = ILDD - cILDD;
% diffavgrL = avgrL - cavgrL;
% diffavgrR = avgrR - cavgrR;
% diffavgsL = avgsL - cavgsL;
% diffavgsR = avgsR - cavgsR;
% diffavgITDD = avgITDD - cavgITDD;
% diffavgILDD = avgILDD - cavgILDD;