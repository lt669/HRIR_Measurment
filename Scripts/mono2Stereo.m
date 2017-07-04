% Make mono speaker IRs into stereo files for testing

clc;
clear;
lPath = 'Audio/Loudspeaker_Audio/Sweeps/Left/';
rPath = 'Audio/Loudspeaker_Audio/Sweeps/Right/';

l = dir(sprintf('%s/*.wav',lPath));
r = dir(sprintf('%s/*.wav',rPath));
fs = 48000;

for k=1:length(l)
    fileNameLeft = sprintf('Audio/Loudspeaker_Audio/Sweeps/Left/%s',l(k).name);
    Left(k,:) = audioread(fileNameLeft);
end

for i=1:length(r)
    fileNameRight = sprintf('Audio/Loudspeaker_Audio/Sweeps/Right/%s',r(i).name);
    Right(i,:) = audioread(fileNameRight);
end
    stereoFile=zeros(144384,2);
for j=1:50
    stereoFile(:,1) = Left(j,:);
    stereoFile(:,2) = Right(j,:);
    name = sprintf('Audio/Loudspeaker_Audio/Sweeps/Stereo/%s',l(j).name);
    audiowrite(name,stereoFile,fs);
end


