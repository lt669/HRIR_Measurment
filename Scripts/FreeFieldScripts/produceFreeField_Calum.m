% Produce Free-Field EQ filters
function [FFHRIR] = produceFreeFieldEQ(subjectName,fileLength,microphone)

disp('--- Running Free Field Filter Script ---');
% Variables for invFIR()
type = 'minphase'; % Minimum Phase response
Nfft = 4096;
Noct = 1; % Octave band smoothing (0 = off, 1 = Octave, 2 = 1/2 Octave etc)
range = [400 20000]; % Range for inversion
reg = [15 20]; % In band and out of band regularisation parameters (dB)
Fs = 48000;

%For the loudspeaker IRs only
fileLengthSpeakers = 56;


% Load File Paths

% HRIRs
hrirPath = sprintf('Audio/HRIR_Raw/%s/',subjectName);
hrirDirectory = dir(sprintf('Audio/HRIR_Raw/%s/*.wav',subjectName));



% Loudspeakers IRs
irPath_Left = strcat('Audio/Loudspeaker_Audio/IR_untrimmed/',microphone{:,1});
irDirectory_Left = dir(strcat('Audio/Loudspeaker_Audio/IR_untrimmed/',microphone{:,1},'/*.wav'));

irPath_Right = strcat('Audio/Loudspeaker_Audio/IR_untrimmed/',microphone{:,2});
irDirectory_Right = dir(strcat('Audio/Loudspeaker_Audio/IR_untrimmed/',microphone{:,2},'/*.wav'));



% Create output directory
mkdir(sprintf('Audio/HRIR_FFEQ/48K_24bit/%s',subjectName));
mkdir(sprintf('Audio/SpeakerIR_InvFilters/48K_24bit/%s',subjectName));

disp(sprintf('Creating Directory: Audio/HRIR_FFEQ/%s',subjectName));
disp(sprintf('Creating Directory: Audio/SpeakerIR_InvFilters/%s',subjectName));

% Create arrays for IR data
%irData = zeros(length(irDirectory_Left),fileLength,2);
%hrirData = zeros(length(hrirDirectory),fileLength,2);

disp('Loading Speaker IR Files...');

% Load in audio files
for i=1:length(irDirectory_Left)
    irFilePath_Left = sprintf('%s/%s',irPath_Left,irDirectory_Left(i).name);
    irFilePath_Right = sprintf('%s/%s',irPath_Right,irDirectory_Right(i).name);
    
    disp(irFilePath_Left );
    
    % Load Files
    irFile_Left = audioread(irFilePath_Left);
    [r,c] = size(irFile_Left);
    disp(sprintf('r: %i, c: %i',r,c));
    irFile_Right = audioread(irFilePath_Right);
    
    irData(i,:,1) = irFile_Left;
    irData(i,:,2) = irFile_Right;
end

disp('Loading HRIR Sweep Files...');
disp(sprintf('Length hrirDirectory: %d',length(hrirDirectory)));

for j=1:length(hrirDirectory)
    
    % Save file names
    filePath = sprintf('%s%s',hrirPath,hrirDirectory(j).name);
    disp(sprintf('j=%d filePath: %s',j,filePath));
    shortName(j) = findName(filePath);
    disp(sprintf('shortName: %s',char(shortName(j))));
    
    % Load Files
    hrirFile = audioread(filePath);
    hrirData(j,:,1) = hrirFile(:,1);
    hrirData(j,:,2) = hrirFile(:,2);
    
end

disp('Calculating Inverse Filters...');
disp(sprintf('Saving InvFilter to: EAD_HRIR/Audio/SpeakerIR_InvFilters/%s',subjectName));

% Calculate and store the frequency response of each ir file
% inverseFilters = calcInvFilter(irData); <- Old function

inverseFilters = zeros(length(irDirectory_Left),Nfft,2);

for k=1:length(irDirectory_Left)
    invInput(:,:) = irData(k,:,:);    
    inverseFilters(k,:,:) = invFIR(type,invInput,Nfft,Noct,Nfft,range,reg,1, Fs);
    saveFilter(:,:) = inverseFilters(k,:,:);
    
    inverseFiltersName = strcat('Audio/SpeakerIR_InvFilters/48K_24bit/',subjectName,'/',char(shortName(k)),'_IF.wav');
    audiowrite(inverseFiltersName,saveFilter,Fs,'BitsPerSample', 24);
end

disp('Applying Free Field Filters');

% Produce output array
FFHRIR = zeros(50,fileLength,2);

% Convolve a source (HRTFs) with the inverse filters
for k=1:50
    filter(:,:) = inverseFilters(k,:,:);
    signal(:,:) = hrirData(k,:,:);
    hrirEQ = zeros(4152:2);
    hrirEQ(:,:) = GK_Convolution(filter,signal,0);
    hrirEQ = hrirEQ(1:fileLength,:); % Truncate 
    
    %hrirEQ = 0.98* hrirEQ ./max(abs(hrirEQ (:))); % and normalise
    FFHRIR(k,:,:) = hrirEQ(:,:);
    hrirEQName = strcat('Audio/HRIR_FFEQ/48K_24bit/',subjectName,'/',char(shortName(k)),'_FFC.wav');
    audiowrite(hrirEQName,hrirEQ,Fs,'BitsPerSample', 24);
end
