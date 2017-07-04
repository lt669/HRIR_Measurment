        %{

    Encodes a HRTF data set into spherical harmonics using Calums stupid
    suite.

%}

clc;

projectName = 'Analysis';
subjectName = 'KU_100_copy';
Fs = 48000;
fsFolder = int2str(round(Fs/1000));
order = 5;

HRTF_Path = strcat('Audio/',projectName,'/HRIR_Trim/',subjectName,'/',fsFolder);
HRTF_Dir = dir(strcat('Audio/',projectName,'/HRIR_Trim/',subjectName,'/',fsFolder,'/*.wav'));

% Define number of Ambisonic channels
channels = (order+1)^2;

% Load in file to assess file length
file = sprintf('%s/%s',HRTF_Path,HRTF_Dir(1).name);
HRTF = audioread(file);
SH_Left = zeros(length(HRTF),channels,length(HRTF_Dir));
SH_Right = zeros(length(HRTF),channels,length(HRTF_Dir));

for k=1:length(HRTF_Dir)
    disp(sprintf('k: %i',k));
    % Load HRTF
    file = sprintf('%s/%s',HRTF_Path,HRTF_Dir(k).name);
    HRTF = audioread(file);
    
    % Extract HRTF angles
    [Azi,Ele] = findAngles(file);
    
    disp(strcat('File: ',file,...
                ' Angles: ',int2str(Azi),...
                ' ',int2str(Ele)));
    
    % Encode into Spherical Harmonics
    SH_Left(:,:,k) = asxWavEncode(HRTF(:,1),Fs,Azi,Ele,order,'bob','supOut');
    SH_Right(:,:,k) = asxWavEncode(HRTF(:,2),Fs,Azi,Ele,order,'bob','supOut');
    
end

[r,c,index] = size(SH_Left);

% Sum all values across 50 spherical harmonic matrices 
for col = 1:c
    disp(sprintf('col: %i',col));
   for row = 1:r
       disp(sprintf('row: %i',row));
        SP_Sum_Left(row,col) = sum(SH_Left(row,col,:));
        SP_Sum_Right(row,col) = sum(SH_Right(row,col,:));
    end
end

%%
% Normalise
SP_Sum_Left = normHRIR(SP_Sum_Left);
SP_Sum_Right= normHRIR(SP_Sum_Right);

%%
SP_filepath = strcat('Audio/',projectName,'/HRIR_SphericalHarmonic/',subjectName,'/',fsFolder);

if 0==exist(SP_filepath,'file')
    mkdir(SP_filepath);
end

SP_filename_left = strcat('',SP_filepath,'/SP_Left.wav');
SP_filename_right = strcat('',SP_filepath,'/SP_Right.wav');

% Output Ambisonic HRTF files
audiowrite(SP_filename_left,SP_Sum_Left,Fs);
audiowrite(SP_filename_right,SP_Sum_Right,Fs);

disp('Done');