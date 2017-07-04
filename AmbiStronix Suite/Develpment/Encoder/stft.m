% close all
% clear all

r = 0.0175; % Array radius
c = 343; % Speed of sound
Fs = 44100; % Sample rate
nfft = 1024;

i = 1;
for f = 0:Fs/nfft:Fs-Fs/nfft
    w = 2*pi*f;
    
    % Compute W compensation response

    Num = 1 + (j*w*r)/c - (1/3)*(((w*r)/c)*((w*r)/c));
    Den = 1 + (1/3)*(j*w*r)/c;
    FW(i) = Num/Den;
    
    % Compute XYZ compensation response
    
    Num = sqrt(6)*(1 + (1/3)*(j*w*r)/c - (1/3)*(((w*r)/c)*((w*r)/c)));
    Den = 1 + (1/3)*(j*w*r)/c;
    FXYZ(i) = Num/Den;
        
    i = i+1;
end

f = 0:Fs/nfft:Fs-Fs/nfft;


% figure; % 0-Fs mag / phase response
% subplot(2, 1, 1); plot(abs(FW));
% subplot(2, 1, 2); plot(angle(FW));

%% Mirror Complex Frequency Response
FW(1, (nfft/2)+2:nfft) = conj(fliplr(FW(1, 2:nfft/2)));
%FW(1, nfft/2+1) = abs(FW(1, nfft/2)); % Wrong!! Use line below instead...
FW(1, nfft/2+1) = 0;

FXYZ(1, (nfft/2)+2:nfft) = conj(fliplr(FXYZ(1, 2:nfft/2)));
FXYZ(1, nfft/2+1) = 0;

%%
%figure; % 0-Fs mag / phase response with 0 as Fs/2
%subplot(2, 1, 1); plot(abs(FW));
%subplot(2, 1, 2); plot(angle(FW));

%% IFFT respose and circshift
timeFW = real(ifft(FW));
shiftTimeFW = circshift(timeFW, nfft/2, 2);
timeFXYZ = real(ifft(FXYZ));
shiftTimeFXYZ = circshift(timeFXYZ, nfft/2, 2);

%%
figure; % Time domain plot
plot(shiftTimeFW);

%% FFT of impulse response for testing
testFW = fft(shiftTimeFW);
testFXYZ = fft(shiftTimeFXYZ);

%%
%figure; % 0-Fs mag / phase response with 0 as Fs/2 after conversion to / from time domain!
%subplot(2, 1, 1); plot(abs(testFW))
%subplot(2, 1, 2); plot(angle(testFW))

%% Plot
figure;
subplot(2,1,1);
semilogx(f,20*log10(abs(testFW)), 'b'); hold on;
semilogx(f,20*log10(abs(testFXYZ)), 'k--o');
legend('W-channel', 'XYZ-channel', 'location', 'northwest');
xlim([0 Fs/2]);
ylim([0 20]);
title(strcat('Magnitude compensation, radius = ', int2str(r*1000), 'mm'));
subplot(2,1,2);
semilogx(f,angle(testFW)*180/pi, 'b'); hold on;
semilogx(f,angle(testFXYZ)*180/pi, 'k--');
legend('W-channel', 'XYZ-channel', 'location', 'northwest');
title(strcat('Phase compensation, radius = ', int2str(r*1000), 'mm'));
xlabel('Frequency');
xlim([0 Fs/2]);
ylim([0 90]);


% phasedifference = angle(testFW) - angle(testFXYZ);
% figure;
% semilogx(f,phasedifference, 'b');