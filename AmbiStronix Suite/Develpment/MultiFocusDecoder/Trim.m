for x = 1:length(phiDeg)
    
    A = 'C:\Users\ca718\Google Drive\PhD\Large_Data\HRTF_Data_Sets\RIG\PinkCapsules\Calum\HRTFs\Calum_MF'; % Load
    B = 'C:\Users\ca718\Google Drive\PhD\Large_Data\HRTF_Data_Sets\RIG\PinkCapsules\Calum\HRTFs\Calum_MF_trimmed'; % Save
    
    filename = dir(sprintf('%s/azi_%d_ele_%d*.wav',...
                                   A, phiDeg(x), thetaDeg(x)));
                               
    [HRTFA, FSA] = audioread(sprintf('%s/%s', A, filename(1).name));

    samp = 2330;
    
    HRTFB = HRTFA(samp:samp+256, :);

    
    audiowrite(sprintf('%s/%s', B, filename(1).name), HRTFB, FSA);
                               
end