function [] = monoToStereoSweeps(projectName,subjectNumber,fs)

    fsFolder = int2str(round(fs/1000));
    
    % Variables
    if fsFolder=='48'
        samples = 144000;
        trimStart = 2414;
    elseif fsFolder=='44'
        samples = 132300;
        trimStart = 2218;
    end
    
    
    
    monoLeftFile = strcat('Audio/',projectName,'/SubjectSweeps_Raw/',subjectNumber,'/',fsFolder,'/Left/');
    
    % Check file exists
    if 0==exist(monoLeftFile ,'file')
        error(sprintf('The Folder... \n\n %s \n\n ...does not exists. Did you move the sweeps into seperate Left and Right folders? Or did you actually record the sweeps in stereo?',monoLeftFile));
    end
    
    monoLeftdir = dir(strcat('Audio/',projectName,'/SubjectSweeps_Raw/',subjectNumber,'/',fsFolder,'/Left/*.wav'));
    
    monoRightFile = strcat('Audio/',projectName,'/SubjectSweeps_Raw/',subjectNumber,'/',fsFolder,'/Right/');
    monoRightdir = dir(strcat('Audio/',projectName,'/SubjectSweeps_Raw/',subjectNumber,'/',fsFolder,'/Right/*.wav'));

    for k=1:50
        disp(k);
        monoLeft = audioread(sprintf('%s%s',monoLeftFile,monoLeftdir(k).name));
        monoRight = audioread(sprintf('%s%s',monoRightFile,monoRightdir(k).name));

        % Check that the file is the correct length (144384 samples)
        len = length(monoLeft);
        if(len>samples)
            disp('Trimming');
            monoLeft = monoLeft(1:samples,:);
            disp(length(monoLeft))
        elseif(len<samples)
            monoLeft = padarray(monoLeft,(samples-len));
        end
        
        len = length(monoRight);
        if(len>samples)
            disp('Trimming');
            monoRight = monoRight(1:samples,:);
            disp(length(monoRight))
        elseif(len<samples)
            monoRight = padarray(monoRight,(samples-len));
        end
        
        Stereo(:,1) = monoLeft(:,1);
        Stereo(:,2) = monoRight(:,2);

        name = strcat('Audio/',projectName,'/SubjectSweeps_Raw/',subjectNumber,'/',fsFolder,'/',monoLeftdir(k).name);
        audiowrite(name,Stereo,48000,'BitsPerSample', 24);
    end

end

   
