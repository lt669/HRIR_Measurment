% Down Sampling Script
function [] = resample_lt()
    clear;
    clc;
    filePath = ('/Users/Lewis 1/Desktop/DiNAR/Unity/Audio Samples/Trader/Horse/');
    fileDir = [dir(strcat(filePath,'48/*.wav'));dir(strcat(filePath,'48/*.mp3'))];
    %nameArray = fileDir(2).name;
    %fileDir = [dir(strcat(filePath,'48/*.wav'))];

    for k=1:length(fileDir)

        file = sprintf('%s48/%s',filePath,fileDir(k).name);%strcat('',filePath,'',fileDir.name(k));

        disp(file);
        [file_in_stereo,fs] = audioread(file);
        file_in = file_in_stereo(:,1);

        %file_44 = resample(file_in,fs,44100);

        file_44 = resample_LT(file_in);
        
        outFile = strcat(filePath,'44/',fileDir(k).name);
        disp(outFile);
        audiowrite(outFile,file_44,44100);
    end
    
    % Resampling function
    function [downsampleFile] = resample_LT(file)
        upsampleFile = interp(file,147);
        downsampleFile = downsample(upsampleFile,160);
    end
end