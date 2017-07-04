%{
    Deconvoles and trims the HRIRs of each of the test subjects
%}

function [decStereoOut] = runSubjectDeconvolve(projectName,subjectName,fileLength,fs,bits)

    disp('--- Running Subject Deconvolve Script ---');
    
    % Variable for FS folder name as a string
    fsFolder = int2str(round(fs/1000));
    
    % Variables
    if fsFolder=='48'
        samples = 144000;
        trimStart = 2414;
    elseif fsFolder=='44'
        samples = 132300;
        trimStart = 2218;
    end
    
    % Paths to audio files
    subjectSweepsPath = strcat('Audio/',projectName,'/SubjectSweeps_Raw/',subjectName,'/',fsFolder); % Used to load files
    subjectDir = dir(sprintf('%s/*.wav',subjectSweepsPath)); % Used to find names of files

    % Check file exists
    if 0==exist(subjectSweepsPath,'file')
        error(sprintf('The Folder... \n\n %s \n\n ...does not exists. Have you selected the correct sampling rate?',subjectSweepsPath));
    end

    
    %Load inverse sweep
    if fsFolder=='48'
        inv = audioread('Audio/Sweeps/InvSweep_20to22050_48000_Pad0s.wav');
    else
        inv = audioread('Audio/Sweeps/InvSweep_20to22050_44100_Pad0s.wav');
    end

    decStereoOut=zeros(50,fileLength,2);

    % Create output directory
    mkdir(strcat('Audio/',projectName,'/HRIR_Raw/',subjectName,'/',fsFolder));
    mkdir(strcat('Audio/',projectName,'/HRIR_Trim/',subjectName,'/',fsFolder));

    % Display Info
    disp(sprintf('Saving Raw HRIRs to:Audio/%s/HRIR_Raw/%s',projectName,subjectName));
    disp(sprintf('Saving Trimmed HRIRs to:Audio/%s/HRIR_Trim/%s',projectName,subjectName));
    disp(sprintf('length(subjectDir) = %d',length(subjectDir)));


    for k = 1:length(subjectDir)

        file = sprintf('%s/%s',subjectSweepsPath,subjectDir(k).name);
        [pathstr,inputName,ext] = fileparts(file);
        fileName(k) = cellstr(char(inputName));
        %disp(strcat('path: ',pathstr, 'Name: ',fileName, 'ext: ',ext));

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
        dec(k,:,:) = deconvolve(inv,sweep);
    end

    % Normalise all HRIRs with respoect to each other
    decNorm = normHRIR(dec);

    for n = 1:length(subjectDir)

        % Write Raw HRIR
        outputRaw(:,:) = decNorm(n,:,:);
        name = strcat('Audio/',projectName,'/HRIR_Raw/',subjectName,'/',fsFolder,'/',char(fileName(n)),'_RawLong.wav');
        audiowrite(name,outputRaw,fs,'BitsPerSample', bits);

        %{
            To cut 0.05s off start of the IR
        
            48k = 2414 (This was the original value used in EAD)
            44.1k = 2218 (This is approximate equivalent to the previously
            used value)
        
            This value should be better determined.
        
        ----------------------------------------------------------------------------
        
            FIND THE HRTF WITH THE SHORTEST ON SET TIME, THEN TRIM ALL
            OTHER HRTFS RELATIVE TO THAT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        ----------------------------------------------------------------------------
        
        %}
        
        
        
        % Write Trimmed HRIR
        decStereoOut(n,:,1) = decNorm(n,trimStart:(trimStart+fileLength-1),1); % What number should this be for 44.1k??
        decStereoOut(n,:,2) = decNorm(n,trimStart:(trimStart+fileLength-1),2);
        name = strcat('Audio/',projectName,'/HRIR_Trim/',subjectName,'/',fsFolder,'/',char(fileName(n)),'_Raw.wav');

        output(:,:) = decStereoOut(n,:,:);
        audiowrite(name,output,fs,'BitsPerSample', bits);
    end


end
