%{
    EAD automatic HRIR capture, equalisation and sofa file production
    script

    NOTE: Can not use DIR() to load audio files, must use direct file path!


    NOTE: Top file should no longer be Documents/MATLAB. To avoid conflict
    with other scripts, have all required scripts in your personal top
    folder.

    All audio is recorded at 48000. Change fs to output all files at the
    desired sample rate. They will be converted at the deconvole stage
    (after deconvolution).
    

%}
clc;
clear;

% IF SWEEPS WERE RECORDED WITH ONE MIC AT A TIME
% 1 = Mono 0 = Stereo
mono2stereo = 0;


% Change for output format
fs = 48000;
bitDepth = 16;

projectName = 'CombFilterTest';
subjectName = 'KU_100';

% Which microphones were used {'Left','Right'}
microphones = {'Yellow','Green'};

% FFEQ Variables
type = 'minphase'; % Minimum Phase response
Nfft = 4096;
Noct = 0; % Octave band smoothing (0 = off, 1 = Octave, 2 = 1/2 Octave etc)
range = [100 20000]; % Range for inversion
reg = [20 10]; % In band and out of band regularisation parameters (dB)



% OPTIONS -- 1 = Yes -- 0 = No

sofaFile = 1;

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


% HRIR_SCRIPT(projectName, subjectName, FS, bitDepth,)

fileLength = 256; % This can/should be changed accordingly


%Convert seperatley recorded HRTFs into a stereo file
if mono2stereo == 1
    monoToStereoSweeps(projectName,subjectName,fs);
end
%%
% Deconvolve HRIR sweeps
rawHRIR = runSubjectDeconvolve(projectName,subjectName,fileLength,fs,bitDepth);

% Apply Free Field Equalisation
%%
FFHRIR = produceFreeField(projectName,subjectName,fileLength,fs,bitDepth,microphones,type,Nfft,Noct,range,reg);

% Apply Diffuse Field Equalisation
% produceDiffuseField(subjectName,fileLength);
%%
% Produce IIR Lookup Table
ITD_Lookup_Table_Generation(projectName,subjectName,fs);


createSOFA(projectName,subjectName,fileLength,fs,bitDepth,FIR_compression,ApplicationName,Organization,AuthorContact,Comment);
FIRtoIIR(projectName,subjectName,fileLength,fs,bitDepth,order,compression_IIR);
