%{
    Quick plot to determine where to trim audio files

    Calculations:
    
    Rig Radius = 1.5
             C = 344m/s
            Fs = 48000

Start of impulse = 1.5/344 = 4.36*10^-3 s
                           = 209.28 samples

Files look like they can be cut at startPoint+56 samples

   
%}
subjectName = 'test_01';
subjectSweepsPath = sprintf('Audio/HRIR_Trim/%s/',subjectName); % Used to load files
subjectDir = dir(sprintf('%s/*.wav',subjectSweepsPath)); % Used to find names of files
file = sprintf('%s%s',subjectSweepsPath,subjectDir(1).name);

x = audioread(file);
plot(x);