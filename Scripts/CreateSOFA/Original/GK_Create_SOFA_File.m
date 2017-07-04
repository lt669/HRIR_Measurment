% Script to put HRIR measurements into SOFA format. This version utilises
% the 'SimpleFreeFieldHRIR' SOFA convention. The data stored is the actual
% FIR filters and angles, as well as associated metadata. 
% The script reads in the correct HRIRs for a defined virtual loudspeaker
% array and stores them in a SOFA container.
% G. Kearney, 2016.


close all;
clear all;
clc;

subjectnumber = 'KEMAR_rig_trimmed2_comp0'; % Subject number to create unique SOFA file

%% Load in HRIRs for Virtual Loudspeaker array of interest

% These are the co-ordinates for a 50ch Lebedev grid
%NOTE: Positive azimuth angles here correspond to negative angles in max (spat)
azimuth = [0,45,135,225,315,0,90,180,270,45,135,225,315,18,72,108,162,198,252,288,342,0,45,90,135,180,225,270,315,18,72,108,162,198,252,288,342,45,135,225,315,0,90,180,270,45,135,225,315,0];
elevation = [90,65,65,65,65,45,45,45,45,35,35,35,35,18,18,18,18,18,18,18,18,0,0,0,0,0,0,0,0,-18,-18,-18,-18,-18,-18,-18,-18,-35,-35,-35,-35,-45,-45,-45,-45,-65,-65,-65,-65,-90];

M = length(azimuth);

%N=256;
N = 577;
hrirs = zeros(M,N,2); % Initialise HRIR array


for i = 1:M
    fileloadname = strcat('trimmedIRs/azi_', int2str(azimuth(i)), '_ele_', int2str(elevation(i)), '_DFC.wav');
    hrirs(i,:,:) = audioread(fileloadname);
end


%% Sofa parameters

% Latency of the created IRs
latency=1; % in samples, must be 1<latency<256

% Data compression (0..uncompressed, 9..most compressed)
compression=0; % results in a nice compression within a reasonable processing time

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
clc;

% Update dimensions
Obj=SOFAupdateDimensions(Obj);

% Fill with attributes
Obj.GLOBAL_ListenerShortName = subjectnumber;
Obj.GLOBAL_History = 'Created with a script';
Obj.GLOBAL_DatabaseName = 'rigTestDataBase';
Obj.GLOBAL_ApplicationName = 'rigTest';
Obj.GLOBAL_ApplicationVersion = SOFAgetVersion('API');
Obj.GLOBAL_Organization = 'University of York';
Obj.GLOBAL_AuthorContact = 'lewis.thresh@york.ac.uk';
Obj.GLOBAL_Comment = '50 source positions. KEMAR subject. Microphones used were KEMAR built in microphones via an MOTU 24IO interface (x3).  Diffuse field compensated minimum phase HRIRs.';
Obj.GLOBAL_License = 'Distributed under Apache Licence.';
Obj.Data.SamplingRate = 48000;


% save the SOFA file
outfilename = strcat('SOFAFiles/', subjectnumber,'.sofa');

disp(['Saving:  ' outfilename]);
Obj=SOFAsave(outfilename, Obj, compression);

%% Check SOFA file is created correctly

ObjTest = SOFAload(outfilename);



