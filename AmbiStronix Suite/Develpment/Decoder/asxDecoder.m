function output = asxDecoder( ambisonicFile,...   
                              HRTFFolder,...
                              SaveAs,...
                              varargin )

% ===================================||=================================   
%                             _______||______________  
%   ____    _  _  ||     ||  //      ||   ___ ____     ____   ||   
%  //  \\  ||\/|| ||__      ||_____  ||  //  //  \\  |//  \\     \\  //   
% ||    || ||  || ||  || ||       || || ||  ||    || ||    || ||   /\
%  \\__//\ ||  || ||__|| ||      //  || ||   \\__//  ||    || || //  \\   
%  _____________________________//  A LIBRARY FOR AMBISONIC EN/DECODING
% ======================================================================
%
% Version 2.2 Feb 2017
%
%       Files Created From f( Ambisonic Source File,
%                             HRTF collection,
%                             Output Filename, 
%                             * Options )
%
% Options:  Decode:           Pseudo-inverse    (pinv)  - default
%           (decode)          Projection        (proj)
%                             Energy Preserving (engpres)
%
%           Weighting:        Basic      (basic)        - default
%           (weighting)       Max-rE     (maxre)
%                             In-Phase   (inphase)
%
%           Degree:           B-Format  (bFormat)       - default auto-
%           (M)               First     (1)               detected
%                             Second    (2)
%                             Third     (3)
%                             Forth     (4)
%                             Fifth     (5)
%
%           Dimension:        3D        (3d)            - default auto-
%           (dim)             2D        (2d)              detected
%
%           Normalisation:    N3D        (n3d)          - default
%           (norm)            SN3D       (sn3d)
%                             Orthonomal (orthonormal)
%
%           Speaker Config.:  quad                      - default   
%           (config)          itu5.1                      dependant on 
%                             octohedron                  degree / 
%                             cube                        dimension
%                             birectangle
%                             hex
%                             dodecahedron
%                             octagon
%                             16chSpherePacking
%                             26ptLebedev
%                             9chCircular
%                             pentakisDodecahedron
%                             12chCircular
%                             pentakisIcosidodecahedron
%                             50ptLebedev
%
% NOTE: CODE CAN CURRENTLY ONLY HANDLE 3D AMBISONICS
%
% Script takes an Ambisonic source file (ambisonicFile) and decodes it 
% for a virtual speaker configuration as chosen by user (options) or 
% selected by default. Loudspeaker signals are output as separate 
% channels of a wav file (ambiSaveAs) in the order the coordinates of 
% the speakers are given in 'iasxLDir'. A specific decding strategy 
% can be selected (options). N3D normalisation is assumed unless 
% otherwise chosen (options). The apropriate degree and dimension of 
% ambisonics can usually be determined from the number of channels in 
% the input file, however this can also be explicitly chosen (options).  
%
% Virtual loudspeaker signals are then convolved with a set of HRTFs
% located in the folder (HRTFFolder). The SADIE database of loudspeaker
% configurations and matching HRTF measurements is used for convinience 
% and to prevent the need for relitively inaccurate HTRF interpolation. 
% The resultant signals from each convolution are finally summed
% together and output as a .wav file with filename (biSaveAs).
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
% 15 / 02 / 2017: v2.1 Split into functions 
% 14 / 03 / 17: Include Decode Matrix weighting opions

%% Handle Inputs and Set Variables

% Get Ambisonic data
    [ambisonic, aFS] = audioread(ambisonicFile);

% Process input options
    [decode, weighting, M, dim, norm, config, K, normTo1, supress] = ...
                             iasxDecoderOptions(ambisonic, varargin{:});

%% Assign Speaker Configuration Co-Ordinates

    [phi, theta, phiDeg, thetaDeg] = iasxLDir(config);

%% Calculate Virtual Loudspeakerpeaker Signals  

% Get Decoding Matrix, D
    [D, ~] = iasxD(phi, theta, decode, M, dim, norm);

% Weight Decoding Matrix
    D = iasxDWeighting(D, weighting, M);
    
% Calculate loudspeaker signals
    lS = (D * ambisonic(:, 1:K)')';
    
%% Decode Virtual Loudspeaker Signals to Binaural
    
    [output, FS] = iasxA2B (lS,...
                            HRTFFolder,...
                            phiDeg,...
                            thetaDeg,...
                            aFS,...
                            normTo1);

%% Output files
    if not(supress)
    % Store degree as string for forming filenames                          
        if isa(M, 'numeric')
            M = num2str(M);
        end
    
    
    % Output loudspeaker signals file
        audiowrite(sprintf('LS_%s_%s_%s_%s_%s_%s.wav',... 
                           SaveAs,...
                           decode,...
                           M,...
                           dim,...
                           norm,...
                           config),...
                   lS, aFS);                       

     % Output binaural file
        audiowrite(sprintf('Bi_%s_%s_%s_%s_%s_%s.wav',... 
                           SaveAs,...
                           decode,...
                           M,...
                           dim,...
                           norm,...
                           config),...
                   output, FS);
    end

end