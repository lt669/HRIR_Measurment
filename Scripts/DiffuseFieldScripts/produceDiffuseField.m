%%%%%%%%%%%%%%%%%%%%%% Compute diffuse field HRIRs %%%%%%%%%%%%%%%%%%%%%%%%
% This script loads a set of RAW HRIRs and computes the diffuse field
% compensated versions. The contribution of each HRIR to the creation of 
% the diffuse field equalisation filter is dependent on each measurements' 
% solid angle. The script extracts the angles from the HRIR filenames with 
% the naming convention 'azi_X_ele_Y_RAW.wav' where X and Y are the azimuth 
% and elevation angles respectively. 
%
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
% Script Adapted by Lewis Thresh 29/12/2016
%%%%%%%%%%%%%%%%%%%%%%% G Kearney 2016 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [s] = produceDiffuseField(subjectNumber,hrirLength)
    
    disp('--- Running Diffuse Field Filter Script ---');

    % Diffuse Field Equalisation parameters
    type = 'minphase'; % Minimum Phase response
    Nfft = 4096;
    Noct = 1; % Octave band smoothing (0 = off, 1 = Octave, 2 = 1/2 Octave etc)
    range = [400 20000]; % Range for inversion
    reg = [15 20]; % In band and out of band regularisation parameters (dB)


    %% Load RAW Files

    disp('Loading RAW HRIRs ...');

    % Create a list of recorded files in directory
    %s = dir(strcat(subjectNumber, 'EAD_HRIR/Audio/HRIR_FFEQ/RAW/48K_24bit/*.wav')); % s is structure array
    s = dir(sprintf('EAD_HRIR/Audio/HRIR_FFEQ/48K_24bit/%s/*.wav',subjectNumber));

    file_list = {s.name}'; % convert the name field into cell array of strings.
    num_meas = length(s); % Number of measurements
  disp(sprintf('DFC: file_list: %s',length(file_list)));
    hrir_raw = zeros(num_meas, hrirLength, 2); % Create empty array for HRIRs
    Azout = zeros(1, num_meas);
    Elout = zeros(1, num_meas);

    for i = 1:num_meas

        % Put each pair of RAW HRIRs into a big array (Reading from 
        % 44K, 16 bit versions here)
        [hrir_raw(i,:,:), Fs] = audioread(char(strcat('EAD_HRIR/Audio/HRIR_FFEQ/48K_24bit/',subjectNumber, ...
            '/', file_list(i,:))));

        % Extract Azimuth and Elevation angle based on filename

        filenamestr = char(file_list(i,:)); % Get current filename
        
        disp(sprintf('DFC: filenamestr: %s',filenamestr));
        
        IndexAzi = strfind(file_list(i,:), 'azi_'); % Find the text 'azi_'
        Azout(i) = sscanf(filenamestr(1,cell2mat(IndexAzi) + ...
            length('azi_'):end), '%g', 1); % Get azimuth value

        IndexEle = strfind(file_list(i,:), 'ele_'); % Find the text 'ele_'
        Elout(i) = sscanf(filenamestr(1,cell2mat(IndexEle) + ...
            length('ele_'):end), '%g', 1); % Get elevation value

    end


    %% Compute Diffuse Field EQ

    disp('Computing Diffuse Field EQ ...');

    hrir_mp = zeros(num_meas, hrirLength-1, 2); % Initialise output HRIRs
    L_AVG = zeros(Nfft,1); % Initialise average left ear response
    R_AVG = zeros(Nfft,1); % Initialise average right ear response

    % Get minimum phase versions of HRIRs
    for i = 1:num_meas
        [~, hrir_mp(i,:,1)] = rceps(hrir_raw(i,1:hrirLength-1,1));
        [~, hrir_mp(i,:,2)] = rceps(hrir_raw(i,1:hrirLength-1,2));
    end

    % Get weights based on solid angle
    s = GetVeroniPlotandSolidAng(Azout, Elout, 0); 
    s = s./max(abs(s(:)));

    % Contribution of HRIR to average response is dependent on solid angle
    for i = 1:num_meas
        HRIR_MP_L = fft(hrir_mp(i,:,1)', Nfft);
        HRIR_MP_R = fft(hrir_mp(i,:,2)', Nfft);
        L_AVG = L_AVG + s(i)*abs(HRIR_MP_L).^2;
        R_AVG = R_AVG + s(i)*abs(HRIR_MP_R).^2;
    end

    L_AVG = sqrt(L_AVG/num_meas);
    R_AVG = sqrt(R_AVG/num_meas);

    df_avg_L = rotate_vect(real(ifft(L_AVG)),Nfft/2);
    df_avg_R = rotate_vect(real(ifft(R_AVG)),Nfft/2);

    df_avg_L = df_avg_L./max(abs(df_avg_L(:)));
    df_avg_R = df_avg_R./max(abs(df_avg_R(:)));


    % Compute Inverse Filters
    L = Nfft;
    window = 1;

    [ih_L]=invFIR(type,df_avg_L,Nfft,Noct,L,range,reg,window, Fs); 
    [ih_R]=invFIR(type,df_avg_R,Nfft,Noct,L,range,reg,window, Fs); 

    % Windowing
    hanlen = Nfft;
    myhan = hanning(hanlen);
    ih_L(end-hanlen/2 +1:end,1) = ih_L(end-hanlen/2 +1:end,1) ...
        .*myhan(end-hanlen/2 +1:end);
    ih_R(end-hanlen/2 +1:end,1) = ih_R(end-hanlen/2 +1:end,1) ...
        .*myhan(end-hanlen/2 +1:end);

    % Check average response results in nice inversion
    comp_L = conv(ih_L, df_avg_L);
    comp_R = conv(ih_R, df_avg_R);

    
    hold off
    figure()
    GK_Freqplot(df_avg_L,Fs, Nfft, 'b', 3, 16, ...
        'Diffuse Field Response', 'Frequency', 'Amplitude',  ...
        [-40 30], [70 20000],  0, 0);
    hold on
    GK_Freqplot(df_avg_R,Fs, Nfft, 'b:', 3, 16, ... 
        'Diffuse Field Response', 'Frequency', 'Amplitude', ...
        [-40 30], [70 20000],  0, 0);
    GK_Freqplot(ih_L,Fs, Nfft, 'r', 3, 16, ... 
        'Diffuse Field Response', 'Frequency', 'Amplitude', ...
        [-40 30], [70 20000],  0, 0);
    GK_Freqplot(ih_R,Fs, Nfft, 'r:', 3, 16, ...
        'Diffuse Field Response', 'Frequency', 'Amplitude', ...
        [-40 30], [70 20000],  0, 0);
    GK_Freqplot(comp_L,Fs, Nfft, 'k', 3, 16, ... 
        'Diffuse Field Response', 'Frequency', 'Amplitude', ...
        [-40 30], [70 20000],  0, 0);
    GK_Freqplot(comp_R,Fs, Nfft, 'k:', 3, 16, ...
        'Diffuse Field Response', 'Frequency', 'Amplitude', ...
        [-40 30], [70 20000],  0, 0);

    legend('Diffuse Field Response (L)', 'Diffuse Field Response (R)', ...
        'Inverse Filter (L)', 'Inverse Filter (R)',  'Result (L)', ...
        'Result (R)','location', 'southeast');


    %% Apply Diffuse Field EQ to dataset

    disp('Applying Diffuse Field EQ...');

    % Initialise array for diffuse field hrirs
    hrir_diff = zeros(num_meas, size(ih_L(:, 1), 1) + ...
        size(hrir_raw(1,:,1)', 1),2);

    for i = 1:num_meas
        stereo_ih(:,1) = ih_L(:,1);
        stereo_ih(:,2) = ih_R(:,1);
        hrir(:,:) = hrir_raw(i,:,:);
        hrir_diff(i,:,:) = GK_Convolution(stereo_ih(:,:), hrir(:,:), 0);
        %hrir_diff(i,:,2) = GK_Convolution(ih_R(:, :), hrir_raw(i,:,2)', 0);
    end

    hrir_diff = hrir_diff(:,1:hrirLength,:); % Truncate 
    %hrir_diff = 0.98* hrir_diff./max(abs(hrir_diff(:))); % and normalise


    %% Write out diffuse field compensated HRIRs

    disp('Writing out diffuse field HRIRs...');

    % Create folder
    mkdir(sprintf('EAD_HRIR/Audio/HRIR_DFFEQ/48K_24bit/%s',subjectNumber));

    for i = 1:num_meas
        outfilename = strcat('EAD_HRIR/Audio/HRIR_DFFEQ/48K_24bit/',subjectNumber,'/azi_', ... 
            int2str(Azout(i)),'_ele_',int2str(Elout(i)), '_DFFC.wav');
        output(:,:) = hrir_diff(i,:,:);
        audiowrite(outfilename, output, ...
            44800, 'BitsPerSample', 24);
    end

end