function [] = headphoneEqualisation(projectName,subjectName,fileLength,fs,bitDepth)

    disp('--- Running headphoneEqualisation Script ---');
    
    % Variable for FS folder name as a string
    fsFolder = int2str(round(fs/1000));
    
    % Variables
    if fsFolder=='48'
        error(sprintf('KU_100 Headphone EQ is currently only avilable at 44.1KHz \n Please change your sampling frequency'));
    end

    % ---LOAD FILE PATHS--- %
    %=======================%

    % HRIRs
    hrirPath = strcat('Audio/Projects/',projectName,'/HRIR_Trim/',subjectName,'/',fsFolder);
    %hrirDirectory = dir(sprintf('%s/*.wav',hrirPath));
    hrirDirectory = dir(strcat('',hrirPath,'/*.wav'));
    
    headphoneEQ = zeros(fileLength,2);
    headphoneEQLength = length(headphoneEQ);
    disp(sprintf('headphoneEQLength : %d',headphoneEQLength ));
    
    % Headphone EQ Inverse IR
    headphoneEQ= audioread('Audio/GlobalAudio/HeadphoneEQ/ku100_SennHD650_HpEQ_inv_filt1.wav');
    headphoneEQLength = length(headphoneEQ);
    disp(sprintf('headphoneEQLength : %d',headphoneEQLength ));  
    
    % Trim to the correct file length
    
    if(fileLength < length(headphoneEQ))
        headphoneEQ = headphoneEQ(1:fileLength,:);
    end
    
    % ----- Convolve ------ %
    %=======================%
    
    for k=1:length(hrirDirectory)

        % Seperate file name from its extension
        file = sprintf('%s/%s',hrirPath,hrirDirectory(k).name);
        [pathstr,inputName,ext] = fileparts(file);
        %fileName(k) = cellstr(char(inputName));
        
        fileName(k) = findName(file);
        
        % Load hrir
        hrir = audioread(strcat('',hrirPath,'/',hrirDirectory(k).name));
        
        % Conovlve with headphone EQ
        hrirEQ(:,1) = conv(hrir(:,1),headphoneEQ(:,1));
        hrirEQ(:,2) = conv(hrir(:,2),headphoneEQ(:,2));
        
        % Retrim to correct length
        hrirEQTrim = hrirEQ(1:fileLength,:);
        
        hrirEQTrimLength = length(hrirEQTrim );
        disp(sprintf('hrirEQTrimLength : %d',hrirEQTrimLength ));  
        
        hrirEQOut(k,:,:) = hrirEQTrim;
        
    end
    disp('--- pre-normalisation ---');
    maxVal(1) = 0; 
    for k=1:50
        x(:,1) = hrirEQOut(k,:,1);
        x(:,2) = hrirEQOut(k,:,2);
        maxVal(2) = max(max(abs(x)));
        if(maxVal(2) > maxVal(1))
            maxVal(1) = maxVal(2);
            index = k;
            disp(k);
        end
    end
    
    disp(strcat('Max: ',hrirDirectory(index).name));
    disp(sprintf('Value(%i): %i',index,maxVal(1)));
    
    % Renormalise with respect to max value
    hrirEQNorm = normHRIR(hrirEQOut);
    [r,c,p] = size(hrirEQNorm);
    
    disp('--- post-normalisation ---');
    % Find name of file with max value
    maxVal(1) = 0; 
    for k=1:r
        x(:,1) = hrirEQNorm(k,:,1);
        x(:,2) = hrirEQNorm(k,:,2);
        maxVal(2) = max(max(abs(x)));
        if(maxVal(2) > maxVal(1))
            maxVal(1) = maxVal(2);
            index = k;
        end
    end
    disp(strcat('Max: ',hrirDirectory(index).name));
    disp(sprintf('Value(%i): %i',index,maxVal(1)));
    
    % ------ Output ------- %
    %=======================%
    % Create output directory
    outputDirectory = strcat('Audio/Projects/',projectName,'/HRIR_HeadphoneEQ/',subjectName,'/',fsFolder);
    mkdir(outputDirectory);
    disp(sprintf('Saving Raw HRIRs to:Audio/%s/HRIR_HeadphoneEQ/%s',projectName,subjectName));
    
    for k=1:length(hrirDirectory)
        output(:,:) = hrirEQNorm(k,:,:);
        name = strcat('',outputDirectory,'/',char(fileName(k)),'_HEQ.wav');
        audiowrite(name,output,fs,'BitsPerSample',bitDepth);
    end
    
    
end