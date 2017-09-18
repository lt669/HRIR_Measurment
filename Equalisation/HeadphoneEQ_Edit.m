%%%%%%%%%%%%%%%%%%%%%% Headphone Equalisation Filter Calculation  %%%%%%%%%
% This script loads a set of recorded sine sweeps and the corresponding
% inverse filter for deconvolution to IRs.
%
% It then calculates an average response of them and makes an inverse
% filter. It is adapted from Gavin Kearney's 2016 Diffuse HRIR script,
% which utilises Angelo Farina's inverse filter functions.
%
% (from GK 2016):
% The diffuse field equalisation filter should be created based on audition
% of the resultant HRTF dataset to ensure good sound quality. Care must be
% taken to make sure that inversion does not create unnescessary
% resonances. Consequently, smoothing can be applied to the inverse filter
% to create a more suitable response without lots of peaks and troughs.
% Note: In the plot of the equalised response - this is the average
% response and therefore a flat response here may likely result in audible
% resonances, whereas a response with some ripple is considered more
% gentle.
%
% The code utilises Angelo Farina's inverse filter functions.
%

% Adapted for Headphone Equalisation based on multiple measurements by Tom
% McKenzie, 29.3.2017

% 
%  close all
%  clear all
%  clc

function [] = HeadphoneEQ_Edit(projectName,subjectName,fileLength,fs,bits)
    %% Deconvolve Sweeps
    

    % Diffuse Field Equalisation parameters
    type = 'minphase'; % Minimum Phase response
    Nfft = 4096;
    Noct = 2; % Octave band smoothing (0 = off, 1 = Octave, 2 = 1/2 Octave etc)
    range = [200 15760]; % Range for inversion
    reg = [15 3]; % In band and out of band regularisation parameters (dB)

    disp('--- Running Headphone Deconvolve Script ---');

    % Variable for FS folder name as a string
    fsFolder = int2str(round(fs/1000));

    % Paths to audio files
    headphoneSweep = strcat('Audio/Projects/',projectName,'/Headphones/Sweeps/',subjectName,'/', fsFolder); % Used to load files
    headphoneSweepDir = dir(sprintf('%s/*.wav',headphoneSweep)); % Used to find names of files

    % Check file exists
    if 0==exist(headphoneSweep,'file')
        error(sprintf('The Folder... \n\n %s \n\n ...does not exists. Have you selected the correct sampling rate?',headphoneSweep));
    end

    %Load inverse sweep
    if fsFolder=='48'
        inv = audioread('Audio/GlobalAudio/Sweeps/InvSweep_20to22050_48000_Pad0s.wav');
    else
        inv = audioread('Audio/GlobalAudio/Sweeps/InvSweep_20to22050_44100_Pad0s.wav');
    end



    % Instantiate headphone sweep array
    hpSweep = zeros(length(audioread(strcat('',headphoneSweep,'/',headphoneSweepDir(1).name))),2);

    for  n=1:5

        hpNameL = strcat('headphoneSweep_L_',int2str((n)),'.wav');
        hpNameR = strcat('headphoneSweep_R_',int2str((n)),'.wav');
        
        % Load sweeps and create stereo array
        hpSweep_L = audioread(strcat('',headphoneSweep,'/',hpNameL));
        hpSweep_R = audioread(strcat('',headphoneSweep,'/',hpNameR));
        
        hpSweep = [hpSweep_L(:,1) hpSweep_R(:,2)];

        % Deconvolve stereo sweep
        hpIR = deconvolve(inv,hpSweep,1);
        % Add to array
        hpIR_array(n,:,:) = hpIR;
    end
    
    num_meas = length(hpIR_array(:,1,1)); % Number of measurements

    %% Compute Diffuse Field EQ

    disp('Computing Headphone EQ ...');

    hp_ir_mp = zeros(num_meas, (fileLength-1), 2); % Initialise output HRIRs
    hp_L_AVG = zeros(Nfft,1); % Initialise average left ear response
    hp_R_AVG = zeros(Nfft,1); % Initialise average right ear response

    % Get minimum phase versions of HRIRs
    for i = 1:num_meas
        [~, hp_ir_mp(i,:,1)] = rceps(hpIR_array(i,1:(fileLength-1),1));
        [~, hp_ir_mp(i,:,2)] = rceps(hpIR_array(i,1:(fileLength-1),2));
    end

    % Contribution of HRIR to average response is dependent on solid angle
    for i = 1:num_meas
        hp_IR_MP_L = fft(hp_ir_mp(i,:,1)', Nfft);
        hp_IR_MP_R = fft(hp_ir_mp(i,:,2)', Nfft);

        % TM - to use the average response based on the solid angle:
        %     L_AVG = L_AVG + s(i)*abs(HRIR_MP_L).^2;
        %     R_AVG = R_AVG + s(i)*abs(HRIR_MP_R).^2;

        % TM - or to not use the solid angle (because my spherical
        % distribution *should* be even (the 492 points on the sphere should be
        % even) then use these two lines instead of the above. This has removed
        % the 's(i)*' bit
        hp_L_AVG = hp_L_AVG + abs(hp_IR_MP_L).^2;
        hp_R_AVG = hp_R_AVG + abs(hp_IR_MP_R).^2;

    end

    hp_L_AVG = sqrt(hp_L_AVG/num_meas);
    hp_R_AVG = sqrt(hp_R_AVG/num_meas);

    hp_df_avg_L = rotate_vect(real(ifft(hp_L_AVG)),Nfft/2);
    hp_df_avg_R = rotate_vect(real(ifft(hp_R_AVG)),Nfft/2);

    hp_df_avg_LNONORM = hp_df_avg_L; % for the plotting
    hp_df_avg_RNONORM = hp_df_avg_R;

    hp_df_avg_L = hp_df_avg_L./max(abs(hp_df_avg_L(:)));
    hp_df_avg_R = hp_df_avg_R./max(abs(hp_df_avg_R(:)));



    % Compute Inverse Filters
    L = Nfft;
    window = 1;

    [hp_ih_L]=invFIR(type,hp_df_avg_L,Nfft,Noct,L,range,reg,window, fs);
    [hp_ih_R]=invFIR(type,hp_df_avg_R,Nfft,Noct,L,range,reg,window, fs);

    % Windowing
    hanlen = Nfft;
    myhan = hanning(hanlen);
    hp_ih_L(end-hanlen/2 +1:end,1) = hp_ih_L(end-hanlen/2 +1:end,1) ...
        .*myhan(end-hanlen/2 +1:end);
    hp_ih_R(end-hanlen/2 +1:end,1) = hp_ih_R(end-hanlen/2 +1:end,1) ...
        .*myhan(end-hanlen/2 +1:end);

        % Save as .wav file
    disp('Saving Inver Headphone IR');
    % Normalise
    hp_ih_stereo = ([hp_ih_L hp_ih_R] ./ max(max(abs(hp_ih_R)),max(abs(hp_ih_L))));

    % Create output directories
    headphoneIRDir = strcat('Audio/Projects/',projectName,'/Headphones/Inv_IR/',subjectName,'/',fsFolder);

    makeDIR(headphoneIRDir);

    % Output
    name = strcat('',headphoneIRDir,'/headphone_average_IR.wav');
    audiowrite(name,hp_ih_stereo,fs,'BitsPerSample',bits);

    
    
    % Check average response results in nice inversion
    hp_comp_L = conv(hp_ih_L, hp_df_avg_L);
    hp_comp_R = conv(hp_ih_R, hp_df_avg_R);

    hold off
    hd650Fig = figure();
    GK_Freqplot(hp_df_avg_L,fs, Nfft, 'r', 1, 12, ...
        '', 'Frequency (Hz)', 'Amplitude (dB)',  ...
        [-40 20], [70 20000],  0, 0);
    hold on
    GK_Freqplot(hp_df_avg_R,fs, Nfft, 'r-.', 1, 12, ...
        '', 'Frequency (Hz)', 'Amplitude (dB)',  ...
        [-40 20], [70 20000],  0, 0);
    GK_Freqplot(hp_ih_L,fs, Nfft, 'g', 1, 12, ...
        '', 'Frequency (Hz)', 'Amplitude (dB)',  ...
        [-40 20], [70 20000],  0, 0);
    GK_Freqplot(hp_ih_R,fs, Nfft, 'g-.', 1, 12, ...
        '', 'Frequency (Hz)', 'Amplitude (dB)',  ...
        [-40 20], [70 20000],  0, 0);
    GK_Freqplot(hp_comp_L,fs, Nfft, 'k', 1, 12, ...
        '', 'Frequency (Hz)', 'Amplitude (dB)',  ...
        [-40 20], [70 20000],  0, 0);
    GK_Freqplot(hp_comp_R,fs, Nfft, 'k-.', 1, 12, ...
        '', 'Frequency (Hz)', 'Amplitude (dB)',  ...
        [-40 20], [70 20000],  0, 0);

    legend('Average HpTF (L)', 'Average HpTF (R)', ...
        'Inverse Filter (L)', 'Inverse Filter (R)',  'Result (L)', ...
        'Result (R)','location', 'southwest');

    set(hd650Fig,'Units','Inches');
    pos = get(hd650Fig,'Position');
    % set(hd650Fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    % print(hd650Fig,'Sennheiser_HD650_EQ','-dpdf','-r0')
    hold off

    % Print left HpTFs and average

    fig111 = figure();
    GK_Freqplot(hp_df_avg_LNONORM,fs, Nfft, 'r', 2, 12, ...
        '', 'Frequency (Hz)', 'Amplitude (dB)',  ...
        [-50 20], [70 20000],  0, 0);
     hold on
    for  i = 1:20
        GK_Freqplot(hp_ir_mp(i,:,1),fs, Nfft, 'k:', 1, 12, ...
        '', 'Frequency (Hz)', 'Amplitude (dB)',  ...
            [-50 20], [70 20000],  0, 0);
    %     if i==1
    %         hold on
    %     end
    end
    GK_Freqplot(hp_df_avg_LNONORM,fs, Nfft, 'r', 2, 12, ...
        '', 'Frequency (Hz)', 'Amplitude (dB)',  ...
        [-50 20], [70 20000],  0, 0);
    % GK_Freqplot(hp_ih_L,fs, Nfft, 'r', 2, 12, ...
    %     'Sennheiser HD650 HpTF (L)', 'Frequency', 'Amplitude', ...
    %     [-50 20], [70 20000],  0, 0);

    legend('Average','location', 'southwest');
    set(fig111,'Units','Inches');
    pos = get(fig111,'Position');
    % set(fig111,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    % print(fig111,'Sennheiser_HD650_HpTF_L','-dpdf','-r0')
    hold off

    % Print Right HpTFs and average
    figure()
    for  i = 1:20
        GK_Freqplot(hp_ir_mp(i,:,2),fs, Nfft, 'k:', 1, 12, ...
            'Sennheiser HD650 HpTF (R)', 'Frequency', 'Amplitude', ...
            [-50 20], [70 20000],  0, 0);
        if i == 1
            hold on
        end
    end
    GK_Freqplot(hp_df_avg_RNONORM,fs, Nfft, 'b', 2, 12, ...
        'Sennheiser HD650 HpTF (R)', 'Frequency', 'Amplitude', ...
        [-50 20], [70 20000],  0, 0);

    hold off



    function [] = makeDIR(path)
        if 0==exist(path,'file')
            disp(strcat('Making directory: ',path));
            mkdir(path);
        end
    end
end