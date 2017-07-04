function virtualAmbiStronixBinauralDecode_v2_0( ambisonicFile,...   
                                                HRTFFolder,...
                                                saveAs,...
                                                varargin )
%
% Version 2.0 Jan 2016
%
% Files Created From f(Ambisonic Source File,
%                      HRTF collection,
%                      Output Filename, 
%                      * Options)
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
    
    [degree,...
     dimension,...
     configuration,...
     norm,...
     decode] = asxDecodeOptions(nChannels, varargin);
    
%% Reset number of ambisonic channels to use based on degree choice
    if strcmp(dimension, '2d')
        nChannels = (2 * degree + 1);
    elseif strcmp(dimension, '3d')
        nChannels = (degree + 1)^2;
    end
    
%% Assign Speaker Configuration Co-Ordinates

    % Set virtual speaker coordinates acording to prefferences
    [phi, theta] = asxSpeakerDirections (configuration);

%% Calculate Virtual Speaker Signals  
    
    [D, ~] = asxD (phi, theta, decode, norm);
    
    speakerSignals = (D * ambisonic(:, 1:nChannels)')';
    
    
%% Decode Virtual Loudspeaker Signals to Binaural
    
    [output, FS] = asxA2B (speakerSignals,...
                                         HRTFFolder,...
                                         phiDeg,...
                                         thetaDeg,...
                                         aFS);

    % Output as .wav file
    if isa(degree, 'numeric')
        degree = num2str(degree);
    end
    audiowrite(sprintf('BiN_%s_%s_%s_degree%s_%s_%s.wav', saveAs,...
                                                      decode,...
                                                      configuration,...
                                                      degree,...
                                                      dimension,...
                                                      norm),...
                                                      output, FS);

end