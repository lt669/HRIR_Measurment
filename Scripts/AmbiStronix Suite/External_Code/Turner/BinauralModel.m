function [IACC, ILD, ITD, CL, cf] = BinauralModel(IR, Fs, N, freqLow)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         BinauralModel is a function for analysing the binaural          %
%         that are used by the human brain in order to localise           %
%                   sound sources in an environment.                      %
%                                                                         %
%         Input Variables:                                                %
%                         IR = The binaural impulse response matrix.      %
%                         Fs = The sampling frequency of the Impulse      %
%                              response.                                  %
%                         N  = The numbe of Gammatone filters used        %
%                              in the auditory model.                     %
%                    freqLow = The center frequency of the lowest         %
%                              frequency band.                            %
%                                                                         %
%         Output Variables:                                               %
%                      xCorr = The cross correlation function of          %
%                              the left and right channel.                %
%                        ILD = The interaural level difference            %
%                              between the left and right ear.            %
%                         CL = The Composite Loudness spectrum            %
%                              of the left and right channel.             %
%                         cf = Center frequencies of the gammatone        %
%                              filters.                                   %
%                                                                         %
%         This code makse use of the AuditoryToolbox MatLab               %
%                    tool box by Malcolm Slaney.                          %
%                                                                         %
%        Slaney, M. (1998). Auditory Toolbox. Palo Alto, CA.              %
%        [Online] Available:                                              %
%        https://engineering.purdue.edu/~malcolm/interval/1998-010/       %
%        [Accessed: 17, November 2016].                                   %
%                                                                         %
%        The binaural model used is based on the work of Ville Pulkki,    %
%        Matti Karjalainen and Jyri Huopaniemi in their paper:            %
%                                                                         %
%        Pulkki, V., Karjalainen, M., & Huopaniemi, J. (1999).            %
%        "Analyzing Virtual Sound Source Attributes Using Binaural        %
%        Auditory Model*." J. Audio Eng. Soc, 47(4), 203–217. [Online]    %
%        Available:                                                       %
%        http://lib.tkk.fi/Diss/2001/isbn9512255324/article5.pdf          %
%        [Accessed: 17, November 2016]                                    %
%                                                                         %
%         Coded by Michael Lovedee-Turner, PhD Music Technology           %
%                       University of York, Audio Lab                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('./AuditoryToolbox')

%% Initialisation of variables to run the model.
% Create Gammatone filters.
% Create the coefficients for the N Gammeatone filters. Starting from
% freqLow and going up to Fs/2. fcoefs is a Nx10 matrix containing the 
% poles, zeros and gain coefficients of the filters.
fcoefs = MakeERBFilters(Fs,N,freqLow);

% Calculte the center frequency of the N Gammatone filters, from freqLow
%  to Fs/2 (Nyquist Frequnecy).
cf = ERBSpace(freqLow, Fs/2, N);

% Create the filter coefficients for a 800Hz low-pass filter.            
%  This filter is applied to the half rectified filtered channels.

lpfCutOff = 800; % low-pass filter cut-off frequency in Hz.

% Get the filter coefficients of the butterworth low-pass filter.
[B,A] = butter(1, lpfCutOff/(Fs/2));

%% Set the maximum lag for the cross-correlation calculation.
maxLag = round(0.0011*Fs);

%% Apply the gammatone filters to the IR.
% Split the left and right chanel out.
xL = IR(:,1);
xR = IR(:,2);

% Apply the ERB filters to the audio.
xLFiltered = ERBFilterBank(xL, fcoefs);
xRFiltered = ERBFilterBank(xR, fcoefs);

% Flip the matricies so it goes from low to high frequency instead of high
% to low.
xLFiltered = flipud(xLFiltered);
xRFiltered = flipud(xRFiltered);
cf = flipud(cf);

%% Half Wave Rectify the filtered Audio.

xLRectified = xLFiltered .* (sign(xLFiltered)+1) ./ 2;
xRRectified = xRFiltered .* (sign(xRFiltered)+1) ./ 2;

%% Apply the low pass filter to the rectified signal.

xLLowPass = filter(B,A,xLRectified')';
xRLowPass = filter(B,A,xRRectified')';

%% Calculate the loudness of the left and right audio channels.

% lengthIR = max(size(IR));
% loudnessL = sqrt(sqrt(sum(xLLowPass.^2, 2) ./ lengthIR));
% loudnessR = sqrt(sqrt(sum(xRLowPass.^2, 2) ./ lengthIR));

% Compute the loudness of the signal in dB
loudnessL = 20 * log10(abs(fft(xLFiltered',Fs)));
% Cut the mirrored out.
loudnessL(end/2:end,:) = [];
% Sum over the frequency band.
loudnessL = 1/(Fs/2) * sum(loudnessL,1);

% Compute the loudness of the signal in dB
loudnessR = 20 * log10(abs(fft(xRFiltered',Fs)));
% Cut the mirrored out.
loudnessR(end/2:end,:) = [];
% Sum over the frequency band.
loudnessR = 1/(Fs/2) * sum(loudnessR,1);

%% Initialise vectors

IACC = zeros(N, maxLag * 2 + 1);

%% Calculate the cross-correlation.
for i = 1  : N
   
    % Create temporary variables to store the Left and right signal being
    % analysed.
    xLTemp = xLLowPass(i,:);
    xRTemp = xRLowPass(i,:);
    
    % Calculate the auto correlation variable
    xAutoCorrelation = sqrt(xLTemp * xLTemp' * xRTemp * xRTemp');
    
    % Calcualte the cross correlation between the left and right ear.
    [xCorr,lags] = xcorr(xLTemp, xRTemp, maxLag);
    
    % Calculate the Interaural cross correlation by dividing the cross
    % correlation with the auto correlation variable.
    IACC(i,:) = xCorr ./ xAutoCorrelation;
    
end
%% Calculate the binaural cues.

% Calculate the interaural time difference.
[~, ITD] = max(IACC, [] , 2);
% Time difference in seconds.
ITD = lags(ITD)./Fs;

% % Calculate the Composite Loudness of the two ears.
% CL = log2(loudnessL + loudnessR)*10+40;
% Compute the composite Loudness.
CL = loudnessL + loudnessR;

% Calculate the ILD as ratio in dB.
ILD = loudnessL ./ loudnessR;
% 
% % Calculate the interaural level difference
% ILD = (log2(loudnessL) - log2(loudnessR)) * 10;

rmpath('./AuditoryToolbox')

end

