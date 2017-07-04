function   D = virtualAmbiStronixBinauralDecode_v1_4( ambisonicFile,...   
                                                      HRTFFolder,...
                                                      saveAs,...
                                                      varargin )
%
% Version 1.4 Dec 2016
%
% Binaural File Created From f(Ambisonic Source File,
%                              HRTF collection,
%                              Output Filename, 
%                              * Options)
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
%           Speaker Config.:  quad                      - default   
%                             itu5.1                      dependant on 
%                             octohedron                  degree / 
%                             cube                        dimension
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
%           Decode:           Projection (proj)
%                             Pseudo-inverse (pseu-inv) - default
%
% NOTE: CODE CAN CURRENTLY ONLY HANDLE 3D AMBISONICS
%
% Script takes an Ambisonic source file (ambisonicFile) and decodes it 
% for a virtual speaker configuration as chosen by user (options) or 
% selected by default. The specific decding strategy can be selected 
% (options). N3D normalisation is assumed unless otherwise 
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
%
% Written by Calum Armstrong, Department of Electronics, The University 
% of York.
%
% Updates ---
% 12 / 12 / 2016: Initial coding
% 13 / 12 / 2016: general improvements
% 20 / 12 / 2016: Implement calculation of general decoding matrix C to
% facilitate different decoding methods
% 03 / 01 / 2017: v1.2 Furthered multi-decode strategy options 
% implementation
% 04 / 01 / 2017: v1.3 Refactored input scanning
% 05 / 01 / 2017: v1.4 Include decode matrtix output. Fix error in 
% Legendre function calculations for angle (0, -90) - it does NOT always 
% equal not(boolean(M)) as previously calculated.

%% Handle Inputs and Set Variables

    % Get Ambisonic data
    [ambisonic, aFS] = audioread(ambisonicFile);
    nChannels = length(ambisonic(1, :)); % Get number on ambisonic 
                                         % channels
                                         
    % Auto-select degree / dimension based on number of Ambisonic 
    % channels
    Ambiguity = 0; % Instantiate variable incase of 
                                  % ambiguous degree / dimension
    switch nChannels
        case 3
            degree = 1;
            dimension = '2d';
            disp('Auto-detection: degree1; 2d');
        case 4
            degree = 1;
            dimension = '3d';
            disp('Auto-detection: degree1; 3d');
        case 5
            degree = 2;
            dimension = '2d';
            disp('Auto-detection: degree2; 2d');
        case 9 % ->             % 9 channels can either be 
            degree = 2;         % degree 2 3D or degree 4 2D.
            dimension = '3d';   % Default to degree 2 3D but
            Ambiguity = 1;      % set ambiguous variable.
            disp('Auto-detection: degree2; 3d');
        case 7
            degree = 3;
            dimension = '2d';
            disp('Auto-detection: degree3; 2d');
        case 16
            degree = 3;
            dimension = '3d';
            disp('Auto-detection: degree3; 3d');
%       case 9
%           degree = 4;
%           dimension = '2d';
        case 25
            degree = 4;
            dimension = '3d';
            disp('Auto-detection: degree4; 3d');
        case 11
            degree = 5;
            dimension = '2d';
            disp('Auto-detection: degree5; 2d');
        case 36
            degree = 5;
            dimension = '3d';
            disp('Auto-detection: degree5; 3d');
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
               disp('User-selection: degree B-Format');
               
           case {1, 2, 3, 4, 5} % Degree
               degree = varargin{k};
               fprintf('User-selection: degree %d\n', varargin{k});
               % Adjust default options and Ambiguity marker if user 
               % clarifies an ambiguous number of channels
               if Ambiguity == 1 && degree == 2
                   dimension = '3d';
                   Ambiguity = 0; 
               elseif Ambiguity == 1 && degree == 4
                   dimension = '3d';
                   Ambiguity = 0;
                   disp('Re-Auto-detection: degree4; 2d');
               end
               
           case {'2d', '3d'} % Dimension  
               dimension = varargin{k};
               fprintf('User-selection: %s\n', varargin{k});
               % Adjust default options and Ambiguity marker if user 
               % clarifies an ambiguous number of channels
               if Ambiguity == 1 && strcmp(dimension, '2d')
                   degree = 4;
                   Ambiguity = 0; 
                   disp('Re-Auto-detection: degree4; 2d');
               elseif Ambiguity == 1 && strcmp(dimension, '3d')
                   degree = 2;
                   Ambiguity = 0;
               end

           case {'quad', 'itu5.1', 'octohedron', 'cube',...
                 'birectangle', 'hex', 'dodecahedron', 'octagon',...
                 '16chSpherePacking', '26ptLebedev', '9chCircular',...
                 'pentakisDodecahedron', '12chCircular',...
                 'pentakisIcosidodecahedron', '50ptLebedev', '7-design'} % Config.           
               configuration = varargin{k};
               fprintf('User-selection: %s\n', varargin{k});
               
           case {'n3d', 'sn3d', 'orthonormal'} % Normalisation
               norm = varargin{k};
               fprintf('User-selection: %s\n', varargin{k});
               
           case {'proj', 'pseu-inv'} % Decode Type
               decode = varargin{k};
               fprintf('User-selection: %s\n', varargin{k});
               
           otherwise % Unexpected input
               error('Invalid options')
       end
    end

    % Output warning if there is still degree / dimension ambiguity
    if Ambiguity == 1;
       warning(['Ambiguous number of Ambisonic channels detected'...
                ' with no dimension or degree specified. Defaulting'...
                ' to 2nd degree 3D decode. Please check!'])
    end    
    
    % Output warning if requested decode (degree / dimension) requires 
    % more channels than are in the ambisonic file
    if strcmp(dimension, '2d') && (2*degree + 1) > nChannels || ...
       strcmp(dimension, '3d') && (degree + 1)^2 > nChannels
       error(['The requested decode cxannot be forfilled as the'...
              ' number of ambisonic channels that would be'...
              ' required exceed the number of channels in the file'])
    end 
    
%% Set Default Values if None Were Given

    % Set default configuration if none selected
    if not(exist('degree', 'var'))
        degree = 1;
        disp('Default degree being used (1)');
    end
    
    % Set default dimension if none selected
    if not(exist('dimension', 'var'))
        dimension = '3d';
        disp('Default dimension being used (3D)');
    end
    
    % Set default normalisation if none selected
    if not(exist('norm', 'var'))
        norm = 'n3d';
        disp('Default normalisation being used (N3D)');
    end
    
    % Set default decode strategy if none selected
    if not(exist('decode', 'var'))
        decode = 'pseu-inv';
        disp('Default decode strategy being used (Pseudo-Inverse)');
    end
    
    % Set default configuration if none selected
    if not(exist('configuration', 'var'))
        configuration = sprintf('default%ddegree%s', degree, dimension);
        disp('Default speaker configuration used');
    end
    
%% Reset number of ambisonic channels to use based on degree choice
    if strcmp(dimension, '2d')
        nChannels = (2 * degree + 1);
    elseif strcmp(dimension, '3d')
        nChannels = (degree + 1)^2;
    end
    
%% Assign Speaker Configuration Co-Ordinates

    % Set virtual speaker coordinates acording to prefferences
    switch configuration
        case {'7-design'}
            phiDeg = [26.0011675216559,333.998832478344,17.1086452559122,342.891354744088,153.998832478344,206.001167521656,162.891354744088,197.108645255912,72.8913547440879,107.108645255912,116.001167521656,63.9988324783441,252.891354744088,287.108645255912,296.001167521656,243.998832478344,32.2544599366034,212.254459936603,302.254459936603,122.254459936603,327.745540063397,147.745540063397,57.7455400633966,237.745540063397];
            thetaDeg = [15.4641512961471,-15.4641512961471,-24.9937030546433,24.9937030546433,-15.4641512961471,15.4641512961471,24.9937030546433,-24.9937030546433,24.9937030546433,-24.9937030546433,15.4641512961471,-15.4641512961471,24.9937030546433,-24.9937030546433,15.4641512961471,-15.4641512961471,60.0253819510733,60.0253819510733,60.0253819510733,60.0253819510733,-60.0253819510733,-60.0253819510733,-60.0253819510733,-60.0253819510733];
        case {'quad', 'default1degree2d'}
            configuration = 'quad';
            phiDeg = [45 135 225 315];
            thetaDeg = [0 0 0 0];
        case 'itu5.1'
            phiDeg = [330 30 0 0 250 110];
            thetaDeg = [0 0 0 0 0 0];

        case {'octohedron', 'default1degree3d'}
            configuration = 'octohedron';
            phiDeg = [0 45 135 225 315 0];
            thetaDeg = [90 0 0 0 0 -90];
        case 'cube'
            phiDeg = [45 135 225 315 45 135 225 315];
            thetaDeg = [35 35 35 35 -35 -35 -35 -35];
        case 'birectangle'
            phiDeg = [90 270 45 135 225 315 90 270];
            thetaDeg = [45 45 0 0 0 0 -45 -45];

        case {'hex', 'default2degree2d'}
            configuration = 'hex';
            phiDeg = [0 60 120 180 240 300];
            thetaDeg = [0 0 0 0 0 0];

        case {'dodecahedron', 'default2degree3d'}
            configuration = 'dodecahedron';
            phiDeg = [180 50 310 118 242 0 180 62 298 130 230 0];
            thetaDeg = [63 46 46 16 16 0 0 -16 -16 -46 -46 -63];

        case {'octagon', 'default3degree2d'}
            configuration = 'octagon';
            phiDeg = [0 45 90 135 180 225 270 315];
            thetaDeg = [0 0 0 0 0 0 0 0];

        case {'26ptLebedev', 'default3degree3d'}
            configuration = '26ptLebedev';
            phiDeg = [0 0 90 180 270 45 135 225 315 0 45 90 135 180 ...
                            225 270 315 45 135 225 315 0 90 180 270 0];
            thetaDeg = [90	45 45 45 45 35 35 35 35	0 0 0 0	0 0	0 0	...
                                  -35 -35 -35 -35 -45 -45 -45 -45 -90];
        case '16chSpherePacking'
            phiDeg = [0 90 180 270 45 135 225 315 0 90 180 270 45   ...
                                                          135 225 315];
            thetaDeg = [51 51 51 51 14 14 14 14 -14 -14 -14 -14 -51 ...
                                                          -51 -51 -51];

        case {'9chCircular', 'default4degree2d'}
            configuration = '9chCircular';
            phiDeg = [0 40 80 120 160 200 240 280 320];
            thetaDeg = [0 0 0 0 0 0 0 0 0 0];

        case {'pentakisDodecahedron', 'default4degree3d'}
            configuration = 'pentakisDodecahedron';
            phiDeg = [90 270 0	180 45 135 225 315 90 270 0	180	32   ...
                      69 111 148 212 249 291 328 0 180 90 270 45 135 ...
                                                  225 315 0 180 90 270];
            thetaDeg = [69	69 58 58 35 35 35 35 32	32 21 21 0 0 0 0 ...
                        0 0 0 0 -21 -21 -32 -32 -35 -35 -35 -35 -58  ...
                                                           -58 -69 -69];

        case {'12chCircular', 'default5degree2d'}
            configuration = '12chCircular';
            phiDeg = 	[0 30 60 90 120 150 180 210 240 270 300 330];
            thetaDeg = [0 0 0 0 0 0 0 0 0 0 0 0];

        case 'pentakisIcosidodecahedron'
            phiDeg = [0 0 180 58 122 238 302 90 270 21	159	201	339  ...
                      58 122 238 302 0 32 90 148 180 212 270 328 58  ...
                      122 238 302 21 159 201 339 90 270	58 122 238   ...
                                                           302 0 180 0];
            thetaDeg = [90	58 58 54 54 54 54 32 32	30 30 30 30	18   ...
                        18 18 18 0 0 0 0 0 0 0 0 -18 -18 -18 -18 -30 ...
                        -30 -30	-30	-32	-32 -54	-54	-54	-54 -58	-58	 ...
                                                                  -90 ];
        case {'50ptLebedev', 'default5degree3d'}
            configuration = '50ptLebedev';
            phiDeg = 	[0 45 135 225 315 0 90 180 270 45 135 225    ...
                         315 18 72 108 162 198 252 288 342 0 45 90   ...
                         135 180 225 270 315 18 72 108 162 198 252   ...
                         288 342 45 135 225 315 0 90 180 270 45 135  ...
                                                             225 315 0];
            thetaDeg = [90 65 65 65	65 45 45 45	45 35 35 35 35 18 18 ...
                        18 18 18 18	18 18 0	0 0	0 0	0 0 0 -18 -18    ...
                        -18 -18 -18 -18 -18 -18 -35 -35 -35 -35 -45  ...
                                       -45 -45 -45 -65 -65 -65 -65 -90];

        otherwise
            error('invalid configuration');
    end

    fprintf('Speaker Configuration: %s\n', configuration);
    
    % Convert angles to radians
    phi = phiDeg * pi / 180;
    theta = thetaDeg * pi / 180;

%% Calculate Re-Encoding Matrix For Decode to Virtual Loudspeaker Array   
    
    % Declare Legendre function
    syms fP_LM(u, l_, m_)
    fP_LM(u, l_, m_) = 1;
    
    % Initiate degree and index
    [m, M] = deal(0); 
    
    % Initiate re-encoding matrix C and Legendre function output
    C = zeros(nChannels, length(phi));
    P_LM = zeros(1, length(phi));
    
    for i = 1:nChannels % For each ambisonic channel

        % Define Legendre Formula
        M_L = abs(M) + m; % Pre-calculation
        fP_LM(u, l_, m_) = (1 / (2^l_ * factorial(l_))) *...
                         ((1 - u^2)^(m_/2)) *...
                         diff((u^2 - 1)^l_,...
                         u ,...
                         M_L);

        % Calculate Legendre function for given degree and index over 
        % all loudspeaker elevations
        if M ~= 0 % Exception that otherwise throws zero divide error *
            logic = abs(sin(theta)) == 1;                                         
            P_LM(logic) = 0;                                                               
            P_LM(not(logic)) = fP_LM(sin(theta(not(logic))), m, abs(M));                                   
        else
            P_LM = fP_LM(sin(theta), m, abs(M));
        end
        
        
        % * fP_LM(1, L, 0) throws a zero divide error. In simple terms this
        % is when the legendre function should equal 0 on the Z-axis. To
        % solve this, these '0' results are entered manually whilst the
        % remaining results are calculated by MatLab.
        
        % Calculate partial normalization factor (not including sqrt(2) 
        % for m != 0)
        switch norm
            case 'n3d'
                N_LM = sqrt((2*m+1)*...
                            factorial(m-abs(M))/factorial(m+abs(M)));
            case 'sn3d'
                N_LM = sqrt(factorial(m-abs(M))/factorial(m+abs(M)));
            case 'orthonormal'
                N_LM = sqrt((2*m+1)/(4*pi)*...
                            factorial(m-abs(M))/factorial(m+abs(M)));
            otherwise
                error('unknown normalisation');
        end

        % Calculate sperical harmonic function coefficients for each
        % loudspeaker position
        if M>0 % For positive indexes...
            r = sqrt(2) * N_LM * P_LM .* cos(M*phi);
        elseif M==0 % For zeroth index...
            r = N_LM * P_LM;
        else % For negative indexes...
            r = sqrt(2) * N_LM * P_LM .* sin(abs(M)*phi);
        end

        % Add co-efficients to re-encoding matrix C
        C(i, :) =  double(r);

        % Calculate the degree and index of the next sequential ACN 
        % ordered spherical harmonic
        M = M+1; % Increase index
        if M > m % If index now exceeds degree...
            m=m+1; % then the degree requires incrementing
            M = -m; % and the index needs re-setting for the new degree
        else
            % Do nothing
        end
    end

%% Select Decode

    switch decode
        case 'proj'
            D = (1/length(phi)) * C';
            
        case 'pseu-inv'
            D = pinv(C);
    end
    
    speakerSignals = (D * ambisonic(:, 1:nChannels)')';
    
%% Decode Virtual Loudspeaker Signals to Binaural
    
    % Read in dummy HRTF to calculate size needed for output vector
    filename = dir(sprintf('%s/azi_*.wav', HRTFFolder));
    [HTRF, hrtfFS] = audioread(sprintf('%s/%s',...
                                         HRTFFolder, filename(1).name));
    
    % Throw error if sampling frequencies of ambisonic file and HRTFs
    % dont match
    if aFS ~= hrtfFS
        error(['Sampling frequency mis-match between ambisonic '...
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
       filename = dir(sprintf('%s/azi_%d_ele_%d*.wav',...
                                   HRTFFolder, phiDeg(i), thetaDeg(i)));
       if length(filename) > 1
           error('Multiple HRTF files for direction: azi %d ele %d', ...
                                                phiDeg(i), thetaDeg(i));
       end
       [HTRF, FS] = audioread(sprintf('%s/%s',...
                                         HRTFFolder, filename(1).name));
       output(:, 1) = output(:, 1) + ...
                                conv(speakerSignals(:, i),  HTRF(:, 1));   
       output(:, 2) = output(:, 2) + ...
                                conv(speakerSignals(:, i),  HTRF(:, 2));
    end

    % Normalise output
    output = output / max(max(abs(output)));
    
    % Output as .wav file
    if isa(degree, 'numeric')
        degree = num2str(degree);
    end
    audiowrite(sprintf('BiN_%s_%s_%s_degree%s_%s_%s.wav', saveAs, ...
                                                      decode, ...
                                                      configuration, ...
                                                      degree, ...
                                                      dimension, ...
                                                      norm), ...
                                                      output, FS);

end