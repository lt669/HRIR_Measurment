% Script adjusts leading 0s to ITD-compensate Multi-focused HRTFs based on
% sample delay values stored in sampDelL and sampDelR variables

for x = 1:length(phiDeg)
    
    A = 'C:\Users\ca718\Google Drive\PhD\Large_Data\HRTF_Data_Sets\RIG\PinkCapsules\WhiteKEMAR\HRTFs\KEMAR_MF'; % Load
    B = 'C:\Users\ca718\Google Drive\PhD\Large_Data\HRTF_Data_Sets\RIG\PinkCapsules\WhiteKEMAR\HRTFs\KEMAR_MF_ITDComp'; % Save
    
    filename = dir(sprintf('%s/azi_%d_ele_%d*.wav',...
                                   A, phiDeg(x), thetaDeg(x)));
                               
    [HRTFA, FSA] = audioread(sprintf('%s/%s', A, filename(1).name));
    
    switch sign(sampDelL(x))
        case 1
            HRTFB(:, 1) = [HRTFA(sampDelL(x)+1:end, 1); zeros(sampDelL(x), 1)];
        case -1
            HRTFB(:, 1) = [zeros(-sampDelL(x), 1); HRTFA(1:end+sampDelL(x), 1)];
    end
    
    switch sign(sampDelR(x))
        case 1
            HRTFB(:, 2) = [HRTFA(sampDelR(x)+1:end, 2); zeros(sampDelR(x), 1)];
        case -1
            HRTFB(:, 2) = [zeros(-sampDelR(x), 1); HRTFA(1:end+sampDelR(x), 2)];
    end
    
    audiowrite(sprintf('%s/%s', B, filename(1).name), HRTFB, FSA);
                               
end