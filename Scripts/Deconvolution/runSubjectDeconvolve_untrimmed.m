function [decStereoOut] = runSubjectDeconvolve(subjectName,fileLength)
disp('--- Running Subject Deconvolve Script ---');
%{
    Deconvoles and trims the HRIRs of each of the test subjects
%}

subjectSweepsPath = sprintf('Audio/SubjectSweeps_Raw/%s/',subjectName); % Used to load files
subjectDir = dir(sprintf('%s/*.wav',subjectSweepsPath)); % Used to find names of files

%Load inverse sweep
inv = audioread('Audio/Sweeps/InvSweep_20to22050_48000_startPad0s_endPad0s.wav');

decStereoOut=zeros(50,fileLength,2);

% Create output directory
mkdir(sprintf('Audio/HRIR_Raw/%s',subjectName));
mkdir(sprintf('Audio/HRIR_Trim/%s',subjectName));

disp(sprintf('Saving Raw HRIRs to:Audio/HRIR_Raw/%s',subjectName));
disp(sprintf('Saving Trimmed HRIRs to:Audio/HRIR_Trim/%s',subjectName));

disp(sprintf('length(subjectDir) = %d',length(subjectDir)));

for k = 1:length(subjectDir)

    file = sprintf('%s%s',subjectSweepsPath,subjectDir(k).name);
    [pathstr,inputName,ext] = fileparts(file);
    fileName(k) = cellstr(char(inputName));
    %disp(strcat('path: ',pathstr, 'Name: ',fileName, 'ext: ',ext));
    [sweep,fs] = audioread(file);
    
    % Check that the file is the correct length (144384 samples)
    len = length(sweep);
    if(len>144384)
        disp('Trimming');
        sweep = sweep(1:144384,:);
        disp(length(sweep))
    elseif(len<144384)
        sweep = padarray(sweep,(144384-len));
    end

    % Deconvolve
    dec(k,:,:) = deconvolve(inv,sweep);
end

% Normalise all HRIRs with respect to each other
decNorm = normHRIR(dec);

for n = 1:length(subjectDir)
    
    % Write Raw HRIR
    outputRaw(:,:) = decNorm(n,:,:);
    name = strcat('Audio/HRIR_Raw/',subjectName,'/',char(fileName(n)),'_RawLong.wav');
    audiowrite(name,outputRaw,fs,'BitsPerSample', 24);

    % Write Trimmed HRIR
    decStereoOut(n,:,1) = decNorm(n,2414:(2414+fileLength-1),1);
    decStereoOut(n,:,2) = decNorm(n,2414:(2414+fileLength-1),2);
    name = strcat('Audio/HRIR_Trim/',subjectName,'/',char(fileName(n)),'_Raw.wav');
    output(:,:) = decStereoOut(n,:,:);
    audiowrite(name,output,fs,'BitsPerSample', 24);
end

    
end
