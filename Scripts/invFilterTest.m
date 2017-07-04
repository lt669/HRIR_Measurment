%{

Check the difference between using the invFIR() function and simply
inverting a fft spectrum

%}

[x,Fs] = audioread('/Users/Lewis 1/Documents/MATLAB/EAD_HRIR/Audio/Loudspeaker_Audio/IR(raw)/ir_azi_18_ele_-18_DFC.wav');
x = x(:,1);
[m,n] = size(x);
xArray = zeros(1,m,n);

type = 'minphase'; % Minimum Phase response
Nfft = 4096;
Noct = 1; % Octave band smoothing (0 = off, 1 = Octave, 2 = 1/2 Octave etc)
range = [400 20000]; % Range for inversion
reg = [15 20]; % In band and out of band regularisation parameters (dB)

[invFIRx]=invFIR(type,x,Nfft,Noct,Nfft,range,reg,1, Fs);
myInv = calcInvFilter(xArray);
X = abs(fft(x,Nfft));
myInvx = X(:).*myInv(:);

figure;
hold on;
GK_Freqplot(x,Fs, Nfft, 'b', 3, 16,'Diffuse Field Response', 'Frequency', 'Amplitude',[-40 30], [70 20000],  0, 0);
%GK_Freqplot(invFIRx,Fs, Nfft, 'r', 3, 16,'Diffuse Field Response', 'Frequency', 'Amplitude',[-40 30], [70 20000], 0, 0);
GK_Freqplot(X,Fs, Nfft, 'g', 3, 16,'Diffuse Field Response', 'Frequency', 'Amplitude',[-40 30], [70 20000], 0, 0);

legend('x', 'invFIRx','myInvx','location', 'southeast');


