% Script seperates and combines 2 seperate HRTF recordings (each focussed
% on a single ear) into a a single Multi-Focused HRTF

for x = 1:length(phiDeg)
    
    L = 'C:\Users\ca718\Google Drive\PhD\Large_Data\HRTF_Data_Sets\RIG\PinkCapsules\WhiteKEMAR\HRTFs\KEMAR_Left' ; % Left channel
    R = 'C:\Users\ca718\Google Drive\PhD\Large_Data\HRTF_Data_Sets\RIG\PinkCapsules\WhiteKEMAR\HRTFs\KEMAR_Right'; % Rite channel
    MF = 'C:\Users\ca718\Google Drive\PhD\Large_Data\HRTF_Data_Sets\RIG\PinkCapsules\WhiteKEMAR\HRTFs\KEMAR_MF'; % Save
    
    filename = dir(sprintf('%s/azi_%d_ele_%d*.wav',...
                                   L, phiDeg(x), thetaDeg(x)));
                               
    [HRTFL, FSL] = audioread(sprintf('%s/%s', L, filename(1).name));
    [HRTFR, ~  ] = audioread(sprintf('%s/%s', R, filename(1).name));
    
    HRTFMF = [HRTFL(:, 1), HRTFR(:, 2)];
    
    audiowrite(sprintf('%s/%s', MF, filename(1).name), HRTFMF, FSL);
                               
end