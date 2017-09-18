%{
    HRIR_MEASURMENT_V2

    New scripts for HRIR measurement and SOFA file production.

    Previous Scripts: Equalises HRTFs using the free field response of the
    binaural microhpones and Genelec speakers in the rig

    New Scripts: Equalises the HRTFs using free field responses of the
    microhpones taken in an anechoic environment with an equalised Genelec
    speaker. Headphone EQ is also incorporated.
    

%}
clc;
clear;

% IF SWEEPS WERE RECORDED WITH ONE MIC AT A TIME
% 1 = Mono 0 = Stereo
mono2stereo = 0;


% Change for output format
fs = 48000;
bitDepth = 16;

projectName = 'AES_Lewis';
subjectName = 'KU_100_AES_13230';

% Which microphones were used {'Left','Right'}
microphones = {'Yellow','Green'};
microphones = 'n';

% FFEQ Variables
type = 'minphase'; % Minimum Phase response
Nfft = 4096;
Noct = 0; % Octave band smoothing (0 = off, 1 = Octave, 2 = 1/2 Octave etc)
range = [100 20000]; % Range for inversion
reg = [20 10]; % In band and out of band regularisation parameters (dB)

% FIR Filter Options
FIR_compression = 9;
ApplicationName = 'EAD Measurements';
Organization = 'University of York';
AuthorContact = 'lewis.thresh@york.ac.uk';
Comment = '50 source positions. Human subject. Microphones used were sennheiser in-ear microphones via an MOTU 24IO interface (x3). Free Field and Diffuse field compensated minimum phase HRIRs.';

% IIR Filter Options
FilterType = 'IIR';
order = 24;
compression_IIR = 0;

fileLength = 256; % This can/should be changed accordingly
fileLength = 13230;
%fileLength = 4096;

%Convert seperatley recorded HRTFs into a stereo file
if mono2stereo == 1
    monoToStereoSweeps(projectName,subjectName,fs);
end
%% DECONVOLUTION

%Headphone_Deconvolution(projectName,subjectName,fileLength,fs,bitDepth);

% Deconvolution and averaging script based off Tom's script
%HeadphoneEQ_Edit(projectName,subjectName,fileLength,fs,bitDepth);

% Deconvolve HRIR sweeps
runSubjectDeconvolve(projectName,subjectName,fileLength,fs,bitDepth,microphones);

%% EQUALISATION

headphoneEqualisation(projectName,subjectName,fileLength,fs,bitDepth);

% Apply Free Field Equalisation
%FFHRIR = produceFreeField(projectName,subjectName,fileLength,fs,bitDepth,microphones,type,Nfft,Noct,range,reg);

% Apply Diffuse Field Equalisation
% produceDiffuseField(subjectName,fileLength);
%%
% Produce IIR Lookup Table
ITD_Lookup_Table_Generation(projectName,subjectName,fs);


createSOFA(projectName,subjectName,fileLength,fs,bitDepth,FIR_compression,ApplicationName,Organization,AuthorContact,Comment);
%FIRtoIIR(projectName,subjectName,fileLength,fs,bitDepth,order,compression_IIR);
