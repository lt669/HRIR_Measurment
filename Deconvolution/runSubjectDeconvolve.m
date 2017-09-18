%{
    Deconvoles and trims the HRIRs of each of the test subjects

    Trimming
    ========
    Time for sound to travel from speaker to rig centre:
    1.5/340 = 0.004411764706
    
    Shorter value used to reduce computation of zero values
    0.001

%}

function [decStereoOut] = runSubjectDeconvolve(projectName,subjectName,fileLength,fs,bits,microphones)

    disp('--- Running Subject Deconvolve Script ---');
    
    % Variable for FS folder name as a string
    fsFolder = int2str(round(fs/1000));
    
    % Variables for 3 seconds in samples (Expected sweep length)
    if fsFolder=='48'
        samples = 144000;
    elseif fsFolder=='44'
        samples = 132300;
    end
    
    % Paths to audio files
    subjectSweepsPath = strcat('Audio/Projects/',projectName,'/Subject_Sweeps/',subjectName,'/',fsFolder); % Used to load files
    subjectDir = dir(sprintf('%s/*.wav',subjectSweepsPath)); % Used to find names of files

    % Check folder exists
    if 0==exist(subjectSweepsPath,'file')
        error(sprintf('The Folder... \n\n %s \n\n ...does not exists. Have you selected the correct sampling rate?',subjectSweepsPath));
    end

    
    %Load inverse sweep
    if fsFolder=='48'
        %inv = audioread('Audio/GlobalAudio/Sweeps/InvSweep_20to24000_48000_Pad0s.wav');
        inv = audioread('Audio/GlobalAudio/Sweeps/InvSweep_20to22050_48000_Pad0s.wav');
    else
        inv = audioread('Audio/GlobalAudio/Sweeps/InvSweep_20to22050_44100_Pad0s.wav');
    end

    % Create output array
    decStereoOut=zeros(50,fileLength,2);

    % Create output directory
    mkdir(strcat('Audio/Projects/',projectName,'/HRIR_Raw/',subjectName,'/',fsFolder));
    mkdir(strcat('Audio/Projects/',projectName,'/HRIR_Trim/',subjectName,'/',fsFolder));

    % Display Info
    disp(sprintf('Saving Raw HRIRs to:Audio/%s/HRIR_Raw/%s',projectName,subjectName));
    disp(sprintf('Saving Trimmed HRIRs to:Audio/%s/HRIR_Trim/%s',projectName,subjectName));
    disp(sprintf('length(subjectDir) = %d',length(subjectDir)));

    
    %=============--------- DECONVOLUTION
    
    for k = 1:length(subjectDir)

        file = sprintf('%s/%s',subjectSweepsPath,subjectDir(k).name);
        [pathstr,inputName,ext] = fileparts(file);
        fileName(k) = cellstr(char(inputName));

        sweep = audioread(file);

        % Check that the file is the correct length (48k = 144,000 Samples 44.1k = 132,300 Samples)
        %(48k = 144384 samples ,44.1k = 132653 samples)


        len = length(sweep);
        [r,c] = size(sweep); % Debug Only

        if(len>samples)
            disp('runSubjectDeconvolve: Trimming');
            sweep = sweep(1:samples,:);
        elseif(len<samples)
            disp('runSubjectDeconvolve: Padding');
            sweep = padarray(sweep,(samples-len),'post');
        end

        % Deconvolve
        dec(k,:,:) = deconvolve(inv,sweep,0);
    end
    
    %=============--------- Convolve with Microphone Inverse IR
    
    %Check if binaural mics were used (not the case using KU_100)
    if microphones ~= 'n'

        leftMicName = strcat('',microphones{:,1},'_Anechoic_InvFilter_1.wav');
        rightMicName = strcat('',microphones{:,2},'_Anechoic_InvFilter_1.wav');

        % Load binaural mic inv IR
        leftMic = audioread(strcat('Audio/GlobalAudio/Microphone_EQ/Microphones/InverseFilters/',leftMicName));
        rightMic = audioread(strcat('Audio/GlobalAudio/Microphone_EQ/Microphones/InverseFilters/',rightMicName));

        % Convolve with HRTFs
        for k = 1:50
            decMicEQ(k,:,1) = conv(dec(k,:,1),leftMic);
            decMicEQ(k,:,2) = conv(dec(k,:,2),rightMic);
        end
        
        % Normalise all HRIRs with respect to each other
        decNorm = normHRIR(decMicEQ);
    else
        decNorm = dec;
    end        
    
    %=============--------- TRIMMING
    
    % Find name of file with max value
    maxVal(1) = 0; 
    for k=1:50
        x(:,1) = decNorm(k,:,1);
        x(:,2) = decNorm(k,:,2);
        maxVal(2) = max(max(abs(x)));
        if(maxVal(2) > maxVal(1))
            maxVal(1) = maxVal(2);
            index = k;
            disp(k);
        end
    end
    disp(strcat('Max: ',subjectDir(index).name));
    disp(sprintf('Value(%i): %i',index,maxVal(1)));
    
    % Initial reference value
    minSort(1) = 1000;
    
    % Find shortest onset time
    for n = 1:length(subjectDir)

        left = decNorm(k,:,1);
        right= decNorm(k,:,2);

        % Normalise for easy direct sound peak detection
        left = left/max(abs(left));
        right = right/max(abs(right));

        % Find peak for trimming
        [peakLeft,loc1] = findpeaks(abs(left),fs,'MinPeakHeight',0.99);
        [peakRight,loc2] = findpeaks(abs(right),fs,'MinPeakHeight',0.99);
        
        % Store Left and Right peak locations
        loc(1) = loc1;
        loc(2) = loc2;
        
        % Find the minimum of the two
        minSort(2) = min(loc);
        
        % If the current minimum is smaller than the previous
        % replace it.
        if minSort(2) < minSort(1)
            minSort(1) = minSort(2);
        end
    end
    
    firstPeak = minSort(1);
    trimStart = round((firstPeak*fs) - (0.0005*fs));
    
    %=============--------- OUTPUT FILES
    
    for n = 1:length(subjectDir)

        % Write Raw HRIR
        outputRaw(:,:) = decNorm(n,:,:);
        name = strcat('Audio/Projects/',projectName,'/HRIR_Raw/',subjectName,'/',fsFolder,'/',char(fileName(n)),'_RawLong.wav');
        audiowrite(name,outputRaw,fs,'BitsPerSample', bits);
        
        % Write Trimmed HRIR
        decStereoOut(n,:,1) = decNorm(n,trimStart:(trimStart+fileLength-1),1); % What number should this be for 44.1k??
        decStereoOut(n,:,2) = decNorm(n,trimStart:(trimStart+fileLength-1),2);
        name = strcat('Audio/Projects/',projectName,'/HRIR_Trim/',subjectName,'/',fsFolder,'/',char(fileName(n)),'_Raw.wav');

        output(:,:) = decStereoOut(n,:,:);
        audiowrite(name,output,fs,'BitsPerSample', bits);
    end


    function [nonZero] = findNonZero(file)
       len = length(file);
       for samp = 1:len
          if(file(samp) > 0)
              nonZero = samp;
              break;
          end
       end
    end
    
end
