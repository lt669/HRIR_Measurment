
% Averages over 10 stereo IRs taken with binaural mics with headphones on

function [] = Headphone_Deconvolution(projectName,subjectName,fileLength,fs,bits)

    disp('--- Running Headphone Deconvolve Script ---');
    
    % Variable for FS folder name as a string
    fsFolder = int2str(round(fs/1000));
    
    % Paths to audio files
    headphoneSweep = strcat('Audio/Projects/',projectName,'/Headphones/Sweeps/',subjectName,'/', fsFolder); % Used to load files
    headphoneSweepDir = dir(sprintf('%s/*.wav',headphoneSweep)); % Used to find names of files
    
    % Check file exists
    if 0==exist(headphoneSweep,'file')
        error(sprintf('The Folder... \n\n %s \n\n ...does not exists. Have you selected the correct sampling rate?',headphoneSweep));
    end
    
    %Load inverse sweep
    if fsFolder=='48'
        inv = audioread('Audio/GlobalAudio/Sweeps/InvSweep_20to22050_48000_Pad0s.wav');
    else
        inv = audioread('Audio/GlobalAudio/Sweeps/InvSweep_20to22050_44100_Pad0s.wav');
    end

    % ---------- Average Sweeps ----------%
    %=====================================%
    
    % Load Sweep to calculate file length
    x = audioread(strcat('',headphoneSweep,'/',headphoneSweepDir(1).name));
    [r,c,p] = size(x);
    nfft = r;
    
    % Load in sweeps
    for k = 1:length(headphoneSweepDir)
       hpSweep(:,:) = audioread(strcat('',headphoneSweep,'/',headphoneSweepDir(k).name));
       fftLeft(:,k) = fft(hpSweep(:,1),nfft);
       fftRight(:,k) = fft(hpSweep(:,2),nfft);
    end
    
    % Calculate power average frequency response
    for k=1:nfft
%         % Power Average
%         fftLeftAverage = sqrt(sum(fftLeft(:,:)^2))/length(fftLeft);
%         fftRightAverage = sqrt(sum(fftRight(:,:)^2))/length(fftRight);
            
        % Average
        fftLeftAverage(k) = sum(fftLeft(k,:))/length(headphoneSweepDir);
        fftRightAverage(k) = sum(fftRight(k,:))/length(headphoneSweepDir);
    end
    
    % Invers FFT back to a sweep
    headphoneSweepAverage(:,1) = ifft(fftLeftAverage);
    headphoneSweepAverage(:,2) = ifft(fftRightAverage);
    
    
    % ---------- Deconvolve ----------%
    %=================================%
    
    % Check that the file is the correct length
    % 48k   = 144,000 Samples
    % 44.1k = 132,300 Samples
    
    % Variables for 3 second sweeps
    if fsFolder=='48'
        samples = 144000;
        trimStart = 2414;
    elseif fsFolder=='44'
        samples = 132300;
        trimStart = 2218;
    end
    
    % Deconvolve average headphone sweep with inverse sweep
    headphoneDec = deconvolve(inv,headphoneSweepAverage);
    disp(sprintf('length(headphoneDec): %d',length(headphoneDec)));
    
    % Seperate channels
    hpL = headphoneDec(:,1);
    hpR = headphoneDec(:,2);
    
    % Normalise left and right for peak detection only
    hpL = hpL/max(hpL);
    hpR = hpR/max(hpR);
    
    % Also Normalise IR for trimming and output
    headphoneDec = headphoneDec/max(max(abs(headphoneDec)));
    
    % Find peak for trimming
    [peakLeft,loc1] = findpeaks(hpL,fs,'MinPeakHeight',0.99);
    [peakRight,loc2] = findpeaks(hpR,fs,'MinPeakHeight',0.99);
    loc = [loc1,loc2];
    
    % Find shortest onset time
    minLoc = min(loc);
    disp(sprintf('minloc: %d',minLoc));
    
    preTrim = 0.01;
    trimStart = (fs*minLoc) - (fs*preTrim);

    % Trim
    headphoneDec_trim = headphoneDec(trimStart:trimStart+fileLength,:);
    
    % Create output directories
    headphoneIRDir = strcat('Audio/Projects/',projectName,'/Headphones/IR/',subjectName,'/',fsFolder);
    headphoneIRTrimDir = strcat('Audio/Projects/',projectName,'/Headphones/IR_Trim/',subjectName,'/',fsFolder);
    
    makeDIR(headphoneIRDir);
    makeDIR(headphoneIRTrimDir);
    
    % Output
    name_untrim = strcat('',headphoneIRDir,'/headphone_average_IR_untrim.wav');
    name_trim = strcat('',headphoneIRTrimDir,'/headphone_average_IR_trim.wav');
    audiowrite(name_untrim,headphoneDec,fs,'BitsPerSample',bits);
    audiowrite(name_trim,headphoneDec_trim,fs,'BitsPerSample',bits);
    
    
    % 1 == Plot response
    plotResponse(0);
    
   
    function [] = makeDIR(path)
        if 0==exist(path,'file')
            disp(strcat('Making directory: ',path));
            mkdir(path);
        end
    end

    function [] = plotResponse(isTrue)
        % --------- PLOTTING --------- %
        %==============================%
        if(isTrue == 1)
            % Plot sweeps
                ax(1) = subplot(4,1,1); plot(hpSweep(:,1,1));
                ax(2) = subplot(4,1,2); plot(headphoneSweepAverage(:,1));
                linkaxes(ax,'xy');

            % Plot Frequency Responses

                % Frequency Vector
                freq = 0:fs/length(hpSweep(:,1,1)):fs/2;

                % Take abs values of FT of the first sweep for plotting and plot on Log
                fftLeftPlotTemp(:,1) = fft(hpSweep(:,1,1),nfft);
                fftLeftPlot(:,1) = fftLeftPlotTemp(1:round(length(hpSweep(:,1,1))/2,1));
                fftLeftPlot(:,1) = 20*log10(abs(fftLeftPlot(:,1)));
                subplot(4,1,3); semilogx(freq(1:length(fftLeftPlot(:,1))),fftLeftPlot(:,1));

                fftLeftAverageTemp = fftLeftAverage;
                fftLeftAverage = fftLeftAverageTemp(1:round(length(hpSweep(:,1,1))/2,1));
                fftLeftAverage = 20*log10(abs(fftLeftAverage));
                subplot(4,1,4); semilogx(freq(1:length(fftLeftAverage)),fftLeftAverage);
        end

    end
    

end