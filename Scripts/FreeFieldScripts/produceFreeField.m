% Produce Free-Field EQ filters
function [FFHRIR] = produceFreeFieldEQ(projectName,subjectName,fileLength,fs,bits,microphone,type,Nfft,Noct,range,reg)

    disp('--- Running Free Field Filter Script ---');
    
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
    
    

    % ---LOAD FILE PATHS--- %

    % HRIRs
    hrirPath = strcat('Audio/',projectName,'/HRIR_Trim/',subjectName,'/',fsFolder);
    hrirDirectory = dir(sprintf('%s/*.wav',hrirPath));

    % Loudspeakers IRs
    irPath_Left = strcat('Audio/Loudspeaker_Audio/IR_untrimmed/',fsFolder,'k_',int2str(bits),'bit/',microphone{:,1});
    irDirectory_Left = dir(strcat('Audio/Loudspeaker_Audio/IR_untrimmed/',fsFolder,'k_',int2str(bits),'bit/',microphone{:,1},'/*.wav'));
    
    irPath_Right = strcat('Audio/Loudspeaker_Audio/IR_untrimmed/',fsFolder,'k_',int2str(bits),'bit/',microphone{:,2});
    irDirectory_Right = dir(strcat('Audio/Loudspeaker_Audio/IR_untrimmed/',fsFolder,'k_',int2str(bits),'bit/',microphone{:,2},'/*.wav'));

    % --------------------- %
    

    % Create output directory
    mkdir(strcat('Audio/',projectName,'/HRIR_FFEQ/',subjectName,'/',fsFolder));
    mkdir(strcat('Audio/',projectName,'/SpeakerIR_InvFilters/',subjectName,'/',fsFolder));

    disp(strcat('Creating Directory: Audio/',projectName,'/HRIR_FFEQ/',subjectName,'/',fsFolder));
    disp(strcat('Creating Directory: Audio/',projectName,'/SpeakerIR_InvFilters/',subjectName,'/',fsFolder));

    
    % Create arrays for IR data
    irData = zeros(length(irDirectory_Left),samples,2);
    hrirData = zeros(length(hrirDirectory),fileLength,2);

    disp('Loading Speaker IR Files...');

    % Load in audio files
    for i=1:length(irDirectory_Left)

        irFilePath_Left = sprintf('%s/%s',irPath_Left,irDirectory_Left(i).name);
        irFilePath_Right = sprintf('%s/%s',irPath_Right,irDirectory_Right(i).name);

        % Load Files
        irFile_Left = audioread(irFilePath_Left);
        irFile_Right = audioread(irFilePath_Right);
        
        [r,c] = size(irFile_Left);
        
        irData(i,:,1) = irFile_Left;
        irData(i,:,2) = irFile_Right;
        
    end

    disp('Loading HRIR Sweep Files...');
    disp(sprintf('Length hrirDirectory: %d',length(hrirDirectory)));

    for j=1:length(hrirDirectory)

        % Save file names
        filePath = sprintf('%s/%s',hrirPath,hrirDirectory(j).name);        
        shortName(j) = findName(filePath);

        % Load Files
        hrirFile = audioread(filePath);
        hrirData(j,:,1) = hrirFile(:,1);
        hrirData(j,:,2) = hrirFile(:,2);

    end

    disp('Calculating Inverse Filters...');
    disp(strcat('Saving InvFilter to: Audio/',projectName,'/SpeakerIR_InvFilters/',subjectName,'/',fsFolder));

    % Calculate and store the frequency response of each ir file
    % inverseFilters = calcInvFilter(irData); <- Old function

    inverseFilters = zeros(length(irDirectory_Left),Nfft,2);
    
    % ------------ Calculate & Save the inverse Filters ------------ %
    
    disp('Producing & Saving Inverse Filters');
    for k=1:length(irDirectory_Left)
        invInput(:,:) = irData(k,:,:);    
        inverseFilters(k,:,:) = invFIR(type,invInput,Nfft,Noct,Nfft,range,reg,1, fs);
        saveFilter(:,:) = inverseFilters(k,:,:);
        inverseFiltersName = strcat('Audio/',projectName,'/SpeakerIR_InvFilters/',subjectName,'/',fsFolder,'/',char(shortName(k)),'_IF.wav');
        audiowrite(inverseFiltersName,saveFilter,fs,'BitsPerSample', bits);
    end
    
    
    
    % ------------ Apply the inverse Filters ------------ %
    disp('Applying Free Field Filters');

    % Produce output array
    FFHRIR = zeros(50,fileLength,2);

    % Convolve a source (HRTFs) with the inverse filters
    for k=1:50
        filter(:,:) = inverseFilters(k,:,:);
        signal(:,:) = hrirData(k,:,:);

        hrirEQ = GK_Convolution(filter,signal,0);
        hrirEQ = hrirEQ(1:fileLength,:); % Truncate 

        %hrirEQ = 0.98* hrirEQ ./max(abs(hrirEQ (:))); % and normalise
        
        FFHRIR(k,:,:) = hrirEQ(:,:); % For retuning array of hrirEQ if necessary    
    end
    
    % ------------ Normalise the FFEQ'ed HRIRs------------ %
    FFHRIR_Norm = normHRIR(FFHRIR);
    
    % ------------ Save the Normalised FFEQ'ed HRIRs------------ %
    for k=1:50
        hrirEQ(:,:) = FFHRIR_Norm(k,:,:);
        hrirEQName = strcat('Audio/',projectName,'/HRIR_FFEQ/',subjectName,'/',fsFolder,'/',char(shortName(k)),'_FFC.wav');
        audiowrite(hrirEQName,hrirEQ,fs,'BitsPerSample', bits);
    end
end
