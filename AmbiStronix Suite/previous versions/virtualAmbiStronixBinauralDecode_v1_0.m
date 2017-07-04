function   virtualAmbiStronixBinauralDecode_v1_0(ambisonicFile,...
                                                 HRTFFolder,...
                                                 saveAs,...
                                                 varargin)
%
% Version 1.0 Dec 2016
%
% Binaural File Created From f(Ambisonic Source File,
%                              HRTF collection,
%                              Output Filename, 
%                              Options)
%
% Options:  Degree:           B-Format  (bFormat)       - default auto-
%                             First     (1)               detected
%                             Second    (2)
%                             Third     (3)
%                             Forth     (4)
%                             Fifth     (5)
%
%           Dimension:        3D        (3d)            - default auto-
%                             2D        (2d)              detected
%
%           Speaker Config.:  quad                      - default depends
%                             itu5.1                      on order /
%                             octohedron                  dimension
%                             cube
%                             birectangle
%                             hex
%                             octagon
%                             26ptLebedev
%                             9chCircular
%                             pentakisDodecahedron
%                             12chCircular
%                             pentakisIcosidodecahedron
%                             50ptLebedev
%
%           Normalisation:    N3D        (n3d)          - default
%                             SN3D       (sn3d)
%                             Orthonomal (orthonormal)
%
% Script takes an Ambisonic source file (ambisonicFile) and decodes it 
% for a virtual speaker configuration as chosen by user (options) or 
% selected by default. N3D normalisation is assumed unless otherwise 
% chosen (options). The apropriate degree and dimension of ambisonics 
% can usually be determined from the number of channels in the input 
% file, however this can also be explicitly chosen (options). 
%
% Virtual loudspeaker signals are then convolved with a set of HRTFs
% located in the folder (HRTFFolder). The SADIE database of loudspeaker
% configurations and matching HRTF measurements is used for convinience 
% and to prevent the need for relitively inaccurate HTRF interpolation. 
% The resultant signals from each convolution are finally summed
% together and output as a .wav file with filename (saveAs).
% Written by Calum Armstrong, Department of Electronics, The University 
% of York.
%
% Updates ---
% 12 / 12 / 2016: Initial coding
% 13 / 12 / 2016: general improvements


%% Handle Inputs and Set Variables

    % Get Ambisonic data
    [ambisonic, aFS] = audioread(ambisonicFile);
    nChannels = length(ambisonic(1, :)); % Get number on ambisonic 
                                         % channels
                                         
    % Auto-select degree / dimension based on number of Ambisonic 
    % channels
    degreeDimensionAmbiguity = 0; % Instantiate variable incase of 
                                  % ambiguous degree / dimension
    switch nChannels
        case 3
            degree = 1;
            dimension = '2d';
        case 4
            degree = 1;
            dimension = '3d';
        case 5
            degree = 2;
            dimension = '2d';
        case 9                            % 9 channels can either be 
            degree = 2;                   % degree 2 3D or degree 4 2D.
            dimension = '3d';             % Default to degree 2 3D but
            degreeDimensionAmbiguity = 1; % set ambiguous variable.
        case 7
            degree = 3;
            dimension = '2d';
        case 16
            degree = 3;
            dimension = '3d';
%       case 9
%           degree = 4;
%           dimension = '2d';
        case 25
            degree = 4;
            dimension = '3d';
        case 11
            degree = 5;
            dimension = '2d';
        case 36
            degree = 5;
            dimension = '3d';
        otherwise
            warning('unusual number of ambisonic channels');
    end

    % Register and handle input options
    for k = 1:length(varargin) % For each input option given...
       switch varargin{k}
           case 'bFormat' % B-Format decode
               temp = ambisonic(:, 2);
               ambisonic(:, 2) = ambisonic(:, 3);
               ambisonic(:, 3) = ambisonic(:, 4);
               ambisonic(:, 4) = temp;
               degree = 1;
           case 1 % 1st degree decode
               degree = 1;
           case 2 % 2nd degree decode
               degree = 2;
               degreeDimensionAmbiguity = 0;   % Remove ambiguity marker   
           case 3 % 3rd degree decode          % if degree is
               degree = 3;                     % explicitely stated
               degreeDimensionAmbiguity = 0;
           case 4 % 4th degree decode
               degree = 4;
           case 5 % 5th degree decode
               degree = 5;
               
           case '2d' % 2D decode               % Remove ambiguity marker   
               dimension = '2d';               % if dimension is
               degreeDimensionAmbiguity = 0;   % explicitely stated
           case '3d' % 3D decode
               dimension = '3d';
               degreeDimensionAmbiguity = 0;
               
           case 'quad'                                   % 2D 1st Degree 
                                            % quad speaker configuration
               configuration = 'quad';
           case 'itu5.1'                                 % 2D 1st Degree 
                                         % ITU 5.1 speaker configuration
               configuration = 'itu5.1';
           case 'octohedron'                             % 3D 1st Degree 
                                       % otohedron speaker configuration
               configuration = 'octohedron';
           case 'cube'                                   % 3D 1st Degree 
                                            % cube speaker configuration
               configuration = 'cube';
           case 'birectangle'                            % 3D 1st Degree 
                                     % birectangle speaker configuration
               configuration = 'birectangle';
           case 'hex'                                    % 2D 2nd Degree 
                                             % hex speaker configuration
               configuration = 'hex';
           case 'dodecahedron'                           % 3D 2nd Degree 
                                    % dodecahedron speaker configuration
               configuration = 'dodecahedron';
           case 'octagon'                                % 2D 3rd Degree 
                                         % octagon speaker configuration
               configuration = 'octagon';
           case '16chSpherePacking'                      % 3D 3rd Degree 
                       % 16 channel sphere packing speaker configuration
               configuration = '16chSpherePacking';
           case '26ptLebedev'                            % 3D 3rd Degree 
                                % 26 point Lebedev speaker configuration
               configuration = '26ptLebedev';
           case '9chCircular'                            % 2D 4th Degree 
                              % 9 channel circular speaker configuration
               configuration = '9chCircular';
           case 'pentakisDodecahedron'                   % 3D 4th Degree 
                            % pentakisDodecahedron speaker configuration
               configuration = 'pentakisDodecahedron';
           case '12chCircular'                           % 2D 5th Degree 
                             % 12 channel circular speaker configuration
               configuration = '12chCircular';
           case 'pentakisIcosidodecahedron'              % 3D 5th Degree 
                       % pentakisIcosidodecahedron speaker configuration  
               configuration = 'pentakisIcosidodecahedron';
           case '50ptLebedev'                            % 3D 5th Degree 
                                % 50 point Lebedev speaker configuration
               configuration = '50ptLebedev';
               
           case 'n3d' % N3D normalisation
               norm = 'n3d';
           case 'sn3d' % SN3D normalisation
               norm = 'sn3d';
           case 'orthonormal' % Orthonormal normalisation
               norm = 'orthonormal';
           otherwise % Unexpected input
               error('Invalid options')
       end
    end

%% Set Default Values if None Were Given

    % Set default configuration if none selected
    if not(exist('degree', 'var'))
        degree = 1;
        disp('Default degree being used (1)');
    end
    
    % Set default configuration if none selected
    if not(exist('dimension', 'var'))
        dimension = '3d';
        disp('Default dimension being used (3D)');
    end
    
    % Set default configuration if none selected
    if not(exist('norm', 'var'))
        norm = 'n3d';
        disp('Default normalisation being used (N3D)');
    end
    
    % Set default configuration if none selected
    if not(exist('configuration', 'var'))
        configuration = sprintf('default%ddegree%s', degree, dimension);
    end
    
    % Output warning if there is still degree / dimension ambiguity
    if degreeDimensionAmbiguity == 1;
       warning(['Ambiguous number of Ambisonic channels detected'...
                ' with no dimension or degree specified. Defaulting'...
                ' to 2nd degree 3D decode. Please check!'])
    end
    
%% Assign Speaker Configuration Co-Ordinates

    % Set virtual speaker coordinates acording to prefferences
    switch configuration
        case {'quad', 'default1degree2d'}
            phiDeg = [45 135 225 315];
            thetaDeg = [0 0 0 0];
        case 'itu5.1'
            phiDeg = [330 30 0 0 250 110];
            thetaDeg = [0 0 0 0 0 0];

        case {'octohedron', 'default1degree3d'}
            phiDeg = [0 45 135 225 315 0];
            thetaDeg = [90 0 0 0 0 -90];
        case 'cube'
            phiDeg = [45 135 225 315 45 135 225 315];
            thetaDeg = [35 35 35 35 -35 -35 -35 -35];
        case 'birectangle'
            phiDeg = [90 270 45 135 225 315 90 270];
            thetaDeg = [45 45 0 0 0 0 -45 -45];

        case {'hex', 'default2degree2d'}
            phiDeg = [0 60 120 180 240 300];
            thetaDeg = [0 0 0 0 0 0];

        case {'dodecahedron', 'default2degree3d'}
            phiDeg = [180 50 310 118 242 0 180 62 298 130 230 0];
            thetaDeg = [63 46 46 16 16 0 0 -16 -16 -46 -46 -63];

        case {'octagon', 'default3degree2d'}
            phiDeg = [0 45 90 135 180 225 270 315];
            thetaDeg = [0 0 0 0 0 0 0 0];

        case {'26ptLebedev', 'default3degree3d'}
            phiDeg = [0 0 90 180 270 45 135 225 315 0 45 90 135 180 ...
                            225 270 315 45 135 225 315 0 90 180 270 0];
            thetaDeg = [90	45 45 45 45 35 35 35 35	0 0 0 0	0 0	0 0	...
                                  -35 -35 -35 -35 -45 -45 -45 -45 -90];
        case '16chSpherePacking'
            phiDeg = [0 90 180 270 45 135 225 315 0 90 180 270 45 ...
                                                        135 225 315];
            thetaDeg = [51 51 51 51 14 14 14 14 -14 -14 -14 -14 -51 ...
                                                          -51 -51 -51];

        case {'9chCircular', 'default4degree2d'}
            phiDeg = [0 40 80 120 160 200 240 280 320];
            thetaDeg = [0 0 0 0 0 0 0 0 0 0];

        case {'pentakisDodecahedron', 'default4degree3d'}
            phiDeg = [90 270 0	180 45 135 225 315 90 270 0	180	32   ...
                      69 111 148 212 249 291 328 0 180 90 270 45 135 ...
                                                  225 315 0 180 90 270];
            thetaDeg = [69	69 58 58 35 35 35 35 32	32 21 21 0 0 0 0 ...
                        0 0 0 0 -21 -21 -32 -32 -35 -35 -35 -35 -58  ...
                                                           -58 -69 -69];

        case {'12chCircular', 'default5degree2d'}
            phiDeg = 	[0 30 60 90 120 150 180 210 240 270 300 330];
            thetaDeg = [0 0 0 0 0 0 0 0 0 0 0 0];

        case 'pentakinIcosidodecahedron'
            phiDeg = [0 0 180 58 122 238 302 90 270 21	159	201	339 ...
                      58 122 238 302 0 32 90 148 180 212 270 328 58 ...
                      122 238 302 21 159 201 339 90 270	58 122 238  ...
                                                          302 0 180 0];
            thetaDeg = [90	58 58 54 54 54 54 32 32	30 30 30 30	18   ...
                        18 18 18 0 0 0 0 0 0 0 0 -18 -18 -18 -18 -30 ...
                        -30 -30	-30	-32	-32 -54	-54	-54	-54 -58	-58	 ...
                                                                  -90 ];
        case {'50ptLebedev', 'default5degree3d'}
            phiDeg = 	[0 45 135 225 315 0 90 180 270 45 135 225   ...
                         315 18 72 108 162 198 252 288 342 0 45 90  ...
                         135 180 225 270 315 18 72 108 162 198 252  ...
                         288 342 45 135 225 315 0 90 180 270 45 135 ...
                                                            225 315 0];
            thetaDeg = [90 65 65 65	65 45 45 45	45 35 35 35 35 18 18 ...
                        18 18 18 18	18 18 0	0 0	0 0	0 0 0 -18 -18    ...
                        -18 -18 -18 -18 -18 -18 -35 -35 -35 -35 -45  ...
                                       -45 -45 -45 -65 -65 -65 -65 -90];

        otherwise
            error('invalid configuration');
    end

    % Convert angles to radians
    phi = phiDeg * pi / 180;
    theta = thetaDeg * pi / 180;

%% Projection Decode to Virtual Loudspeaker Array
    
    % Declare Legendre function
    syms fP_LM(u, l, m)

    % Initiate degree and index
    [L, M] = deal(0); 
    
    % Initiate variable to hold all loudspeaker signals
    speakerSignals = zeros(length(ambisonic(:, 1)), length(phi));

    for i = 1:nChannels % For each ambisonic channel

        % Define Legendre Formula
        M_L = abs(M) + L; % Pre-calculation
        fP_LM(u, l, m) = (1 / (2^l * factorial(l))) *...
                         ((1 - u^2)^(m/2)) *...
                         diff((u^2 - 1)^l,...
                         u ,...
                         M_L);

        % Calculate Legendre function for given degree, index, elevation
        logic =  abs(sin(theta)) == 1; % Exception that otherwise...
                                       % throws zero divide error
        P_LM(logic) = double(not(boolean(M)));
        P_LM(not(logic)) = fP_LM(sin(theta(not(logic))), L, abs(M)); 
                                                   % General calculation   


        % Calculate partial normalization factor (not including sqrt(2) 
        % for m != 0)
        switch norm
            case 'n3d'
                N_LM = sqrt((2*L+1)*...
                            factorial(L-abs(M))/factorial(L+abs(M)));
            case 'sn3d'
                N_LM = sqrt(factorial(L-abs(M))/factorial(L+abs(M)));
            case 'orthonomal'
                N_LM = sqrt((2*L+1)/(4*pi)*...
                            factorial(L-abs(M))/factorial(L+abs(M)));
            otherwise
                error('unknown normalisation');
        end

        % Calculate sperical harmonic function coefficients
        if M>0 % For positive indexes...
            r = sqrt(2) * N_LM * P_LM .* cos(M*phi);
        elseif M==0 % For zeroth index...
            r = N_LM * P_LM;
        else % For negative indexes...
            r = sqrt(2) * N_LM * P_LM .* sin(abs(M)*phi);
        end
        
        % Convert to double format
        r = double(r);

        % Scale channel by coefficient and sum with previous
        speakerSignals = speakerSignals + ambisonic(:, i) * r;

        % Calculate the degree and index of the next sequential ACN 
        % ordered spherical harmonic
        M = M+1; % Increase index
        if M > L % If index now exceeds degree...
            L=L+1; % then the degree requires incrementing
            M = -L; % and the index needs re-setting for the new degree
        else
            % Do nothing
        end
    end

%% Decode Virtual Loudspeaker Signals to Binaural
    
    % Read in dummy HRTF to calculate size needed for output vector
    [HTRF, hrtfFS] = audioread(sprintf('%s/azi_0_ele_0_RAW.wav',...
                                       HRTFFolder));
    
    % Throw error if sampling frequencies of ambisonic file and HRTFs
    % dont match
    if aFS ~= hrtfFS
        error(['Sampling frequency mis-match between ambisonic'...
               'file and HRTFs']);
    end
    
    % Instantiate output vector
    output = zeros(max([...
                   length(speakerSignals(:, 1))+length(HTRF(:, 1))-1,...   
                   length(speakerSignals(:, 1)),...
                   length(HTRF(:, 1))]), 2);

    % Loop through each loudspeaker signal, convolve with HRTF 
    % of matching angle, and sum outputs
    for i = 1:length(phi)
       filename = sprintf('%s/azi_%d_ele_%d_RAW.wav',...
                          HRTFFolder, phiDeg(i), thetaDeg(i));
       [HTRF, FS] = audioread(filename);
       output(:, 1) = output(:, 1) + ...
                                conv(speakerSignals(:, i),  HTRF(:, 1));   
       output(:, 2) = output(:, 2) + ...
                                conv(speakerSignals(:, i),  HTRF(:, 2));
    end

    % Normalise output
    output = output / max(max(output));
    
    % Output as .wav file
    audiowrite(sprintf('%s.wav', saveAs), output, FS);

end