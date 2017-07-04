%{
Plot the results of the FFEQ IRs aginst the original IRs
%}
k = 1;
xPath = '/Users/Lewis 1/Documents/MATLAB/EAD_HRIR/Audio/HRIR_Trim/test_01/';
xD = dir('/Users/Lewis 1/Documents/MATLAB/EAD_HRIR/Audio/HRIR_Trim/test_01/*.wav');
x = audioread(sprintf('%s%s',xPath,xD(k).name));

eqPath = '/Users/Lewis 1/Documents/MATLAB/EAD_HRIR/Audio/HRIR_FFEQ/48k_24bit/test_01/';
dEQ = dir('/Users/Lewis 1/Documents/MATLAB/EAD_HRIR/Audio/HRIR_FFEQ/48k_24bit/test_01/*.wav');
xEQ = audioread(sprintf('%s%s',eqPath,dEQ(k).name));

invPath = '/Users/Lewis 1/Documents/MATLAB/EAD_HRIR/Audio/HRIR_InvFilters/test_01/';
dInv = dir('/Users/Lewis 1/Documents/MATLAB/EAD_HRIR/Audio/HRIR_InvFilters/test_01/*.wav');
xInv = audioread(sprintf('%s%s',invPath,dInv(k).name));

diffInvPath = '/Users/Lewis 1/Documents/MATLAB/EAD_HRIR/Audio/HRIR_DFFEQ/48k_24bit/test_01/';
dDiff = dir('/Users/Lewis 1/Documents/MATLAB/EAD_HRIR/Audio/HRIR_DFFEQ/48k_24bit/test_01/*.wav');
xDiff = audioread(sprintf('%s%s',diffInvPath,dDiff(k).name));

Nfft = 4096;
Fs = 48000;
norm = 0;

%{
figure;
hold on;
%plot(x);
plot(xEQ);
%}

figure;
hold on;
GK_Freqplot(x(:,1),Fs, Nfft, 'b', 3, 16,'Diffuse Field Response', 'Frequency', 'Amplitude',[-40 30], [70 20000],  norm, 0);
GK_Freqplot(x(:,2),Fs, Nfft, 'b:', 3, 16,'Diffuse Field Response', 'Frequency', 'Amplitude',[-40 30], [70 20000],  norm, 0);

GK_Freqplot(xEQ(:,1),Fs, Nfft, 'g', 3, 16,'Diffuse Field Response', 'Frequency', 'Amplitude',[-40 30], [70 20000], norm, 0);
GK_Freqplot(xEQ(:,2),Fs, Nfft, 'g:', 3, 16,'Diffuse Field Response', 'Frequency', 'Amplitude',[-40 30], [70 20000], norm, 0);

GK_Freqplot(xInv(:,1),Fs, Nfft, 'r', 3, 16,'Diffuse Field Response', 'Frequency', 'Amplitude',[-40 30], [70 20000], norm, 0);
GK_Freqplot(xInv(:,2),Fs, Nfft, 'r:', 3, 16,'Diffuse Field Response', 'Frequency', 'Amplitude',[-40 30], [70 20000], norm, 0);

GK_Freqplot(xDiff(:,1),Fs, Nfft, 'k', 3, 16,'Diffuse Field Response', 'Frequency', 'Amplitude',[-40 30], [70 20000], norm, 0);
GK_Freqplot(xDiff(:,2),Fs, Nfft, 'k:', 3, 16,'Diffuse Field Response', 'Frequency', 'Amplitude',[-40 30], [70 20000], norm, 0);

legend('IR_L','IR_R', 'invFIR_L','invFIR_R','xInvFilter_L','xInvFilter_R','xDiff_L','xDiff_R','location', 'southeast');
