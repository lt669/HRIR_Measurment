% Script to put HRIR measurements into SOFA format. This version utilises
% the 'SimpleFreeFieldHRIR' SOFA convention. The data stored is the actual
% FIR filters and angles, as well as associated metadata. 
% The script reads in the correct HRIRs for a defined virtual loudspeaker
% array and stores them in a SOFA container.
% G. Kearney, 2016.


function[] = createSOFA(projectName,subjectNumber,fileLength,fs,bit,FIR_compression,ApplicationName,Organization,AuthorContact,Comment)
    
    disp('');    
    disp('---Running createSOFA---');

    % Variable for FS folder name as a string
    fsFolder = int2str(round(fs/1000));
    for h=1:2
        if h == 1
            sofaFileSource = strcat('Audio/Projects/',projectName,'/HRIR_Trim/');
            extension = '_Raw.wav';
            sofaOutputFile = strcat('Audio/Projects/',projectName,'/SOFAFiles/SOFA_RAW/');
            disp(strcat(''))
        elseif h == 2
            sofaFileSource = strcat('Audio/Projects/',projectName,'/HRIR_HeadphoneEQ/');
            extension = '_HEQ.wav';
            sofaOutputFile = strcat('Audio/Projects/',projectName,'/SOFAFiles/SOFA_HEQ/');
            disp(strcat(''))
        end

        % These are the co-ordinates for a 50ch Lebedev grid
        %NOTE: Positive azimuth angles here correspond to negative angles in max (spat)
        azimuth = [0,45,135,225,315,0,90,180,270,45,135,225,315,18,72,108,162,198,252,288,342,0,45,90,135,180,225,270,315,18,72,108,162,198,252,288,342,45,135,225,315,0,90,180,270,45,135,225,315,0];
        elevation = [90,65,65,65,65,45,45,45,45,35,35,35,35,18,18,18,18,18,18,18,18,0,0,0,0,0,0,0,0,-18,-18,-18,-18,-18,-18,-18,-18,-35,-35,-35,-35,-45,-45,-45,-45,-65,-65,-65,-65,-90];

        M = length(azimuth);

        %N=256;
        N = fileLength;
        hrirs = zeros(M,N,2); % Initialise HRIR array

        for i = 1:M
            fileloadname = strcat('',sofaFileSource,'',subjectNumber,'/',fsFolder,'/azi_', int2str(azimuth(i)), '_ele_', int2str(elevation(i)), '',extension);
            [hrirs(i,:,:),fs] = audioread(fileloadname);
        end


        %% Sofa parameters

        % Latency of the created IRs
        latency=1; % in samples, must be 1<latency<256

        % Data compression (0..uncompressed, 9..most compressed)
        compression=FIR_compression; % results in a nice compression within a reasonable processing time

        % Get an empy conventions structure
        Obj = SOFAgetConventions('SimpleFreeFieldHRIR');

        % Fill data with data
        Obj.Data.IR = NaN(length(azimuth),2,N); % data.IR must be [M R N]

        %% Sort and load data 

        % First data sort for SOFA (by azimuth)

        % sortindex gives the index of the data before it was sorted
        [g sortindex ] = sort(azimuth);

        %Sorts the azimuth angles so they are in ascending order
        for i=1:M
            j = sortindex(i);
            hrirSOFAtemp(i,:,:) = hrirs(j,:,:);
            Aztemp(i) = azimuth(j);
            Eltemp(i) = elevation(j);
        end

        % Second data sort (by elevation) now that the azimuth angles are in the
        % right order
        [g sortindex ] = sort(Eltemp);

        for i=1:M
            j = sortindex(i);
            hrirSOFA(i,:,:) = hrirSOFAtemp(j,:,:);
            Azsofa(i) = Aztemp(j);
            Elsofa(i) = Eltemp(j);
        end

        for i=1:M
            Obj.Data.IR(i,1,:)= hrirSOFA(i,:,1);
            Obj.Data.IR(i,2,:)= hrirSOFA(i,:,2);
            Obj.SourcePosition(i,:)=[Azsofa(i) Elsofa(i) 1.5];
        end

        %%
        %clc;

        % Update dimensions
        Obj=SOFAupdateDimensions(Obj);

        % Fill with attributes
        Obj.GLOBAL_ListenerShortName = subjectNumber;
        History = 'Created with a script';
        Obj.GLOBAL_DatabaseName = subjectNumber;
        Obj.GLOBAL_ApplicationName = ApplicationName;
        Obj.GLOBAL_ApplicationVersion = SOFAgetVersion('API');
        Obj.GLOBAL_Organization = Organization;
        Obj.GLOBAL_AuthorContact = AuthorContact;
        Obj.GLOBAL_Comment = Comment;
        Obj.GLOBAL_License = 'Distributed under Apache Licence.';
        Obj.Data.SamplingRate = fs;


        % save the SOFA file
        outputDirectory = strcat('',sofaOutputFile,'',subjectNumber);

        if 0==exist(outputDirectory,'file')
            mkdir(outputDirectory);
        end

        outfilename = strcat('',outputDirectory,'/',subjectNumber,'_',int2str(fileLength),'order_fir_',fsFolder,'k_',int2str(bit),'bit.sofa');

        disp(['Saving:  ' outfilename]);
        Obj=SOFAsave(outfilename, Obj, compression);

        %% Check SOFA file is created correctly

        ObjTest = SOFAload(outfilename);
    end
    disp('');    
    disp('---Exiting createSOFA---');
end


