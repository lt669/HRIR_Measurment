function [ decode,...
           weighting,...
           M,...
           dim,...
           norm,...
           config,...
           K,...
           normTo1,...
           supress] = iasxDecoderOptions( ambisonic, varargin )

% Version 1.0 Feb 2017
%
% [ Decoding Strategy,
%   Degree ('Ambisonic Order' of decode)
%   dimension,
%   Normalisation,
%   Loudspeaker configuration,
%   No. of Ambisonic Channels in decode ] From f( Ambisonic file,
%                                                 * Options )
%
% Script processes the optional user inputs from 'asxDecoder' and 
% combines them with auto-detection algorythms to assign values to all 
% variables needed for an Ambisonic decode. 
%
% Written by Calum Armstrong, Department of Electronics, The University 
% of York.
%
% Updates ---
% 15 / 02 / 2017: Initial coding

%% REGISTER AND HANDLE USER INPUT OPTIONS

    K = length(ambisonic(1, :)); % Get number on ambisonic channels

    for k = 1:length(varargin) % For each input option given...
        switch varargin{k} % Scan input and assign to relavent variables
            case {'proj', 'pinv', 'engpres'} % Decode Type
               decode = varargin{k};
               fprintf('User-selection: %s\n', varargin{k});
            
            case {'basic', 'maxre', 'inphase'}
                weighting = varargin{k};
                fprintf('User-selection: %s\n', varargin{k});
                
            case 'bFormat' % B-Format decode
                % Convert to ACN N3D format
                temp = ambisonic(:, 2);
                ambisonic(:, 1) = sqrt(2) * ambisonic(:, 1);
                ambisonic(:, 2) = sqrt(3) * ambisonic(:, 4);
                ambisonic(:, 4) = sqrt(3) * ambisonic(:, 3);
                ambisonic(:, 3) = sqrt(3) * temp;
                M = 1;
                norm = 'n3d';
                disp('User-selection: degree B-Format');
               
            case {1, 2, 3, 4, 5} % Degree
               M = varargin{k};
               fprintf('User-selection: degree %d\n', varargin{k});

            case {'2d', '3d'} % Dimension  
               dim = varargin{k};
               fprintf('User-selection: %s\n', varargin{k});

            case {'n3d', 'sn3d', 'orthonormal'} % Normalisation
               norm = varargin{k};
               fprintf('User-selection: %s\n', varargin{k});
            
            case {'raw'}
                normTo1 = 0;
                fprintf('User-selection: %s\n', varargin{k});
                
            case {'supOut'}
                supress = 1;
                fprintf('User-selection: %s\n', varargin{k});
                
            otherwise % Loudspeaker Configuration
                config = varargin{k};
                % Test whether configuration actually exists
                try
                    iasxLDir(config);
                catch ME % Throw error if input doesnt exist in 
                         % loudspeaker directions file
                    if (strcmp(ME.identifier,'VerifyInput:InvalidInput'))
                        error('Invalid options: %s', varargin{k});
                    else
                        rethrow(ME)
                    end
                end
                
               fprintf('User-selection: %s\n', varargin{k});
               
       end
    end

%% AUTO-SELECT DEGREE / DIMENSION IF NONE EXPLICITLY GIVEN

    if not(exist('dim', 'var')) % If there was no dimension given...
        switch K % Select most likely dimension based on number of 
                 % ambisonic channels
            case {3, 5, 7, 11}
                dim = '2d'; disp('Auto-detection: 2D');
            case {4, 16, 25, 36}
                dim = '3d'; disp('Auto-detection: 3D');
            case {9} % 9 channels are used for both degree 2 3D, and 
                     % degree 4 2D. If a degree is given, select 
                     % apropriate dimention
                if exist('M', 'var') && (M == 4 || M == 3)
                    dim = '2d'; disp('Auto-detection: 2D');
                elseif exist('M', 'var') && M == 2
                    dim = '3d'; disp('Auto-detection: 3D');
                else % If no degree was given, default to 3D but output 
                     % warning
                    dim = '3d'; disp('Auto-detection: 3D');
                    warning(['Ambiguous number (9) of Ambisonic'...
                             ' channels detected with no clarifying'...
                             ' dimension or degree specified.'...
                             ' Defaulting to 2nd degree 3D decode.'...
                             ' Please check!'])
                end
            otherwise % If dimention cannot be auto-detected, output 
                      % warning
                warning(['Unrecognised number of Ambisonic Channels'...
                        ' in file. Unable to auto-detect degree /'...
                        ' dimention']);
        end
    end

    % If there was no degree given, but a dimention has been selected...   
    if not(exist('M', 'var')) && exist('dim', 'var') 
        switch dim % Calculate degree based on dimension and number 
                   % of amsbsonic channels
            case '2d'
                M = (K-1) / 2;
            case '3d'
                M = sqrt(K) - 1;
        end
        
        if rem(M, 1) == 0
            fprintf('Auto-detection: degree %d', M);
        else % If there is an unusual number of Ambisonic Channels that 
             % has led to an ambiguous degree calculation, output 
             % warning
            warning(['Number of Ambisonic Channels in file does not'...
                     ' match the dimention specified.']);
            M = floor(M);
            fprintf('Auto-detection: degree %d', M);
        end
    end    
    
%% VALIDATE OPTIONS    
    
    % Output warning if requested decode (degree / dimension) requires 
    % more channels than are in the ambisonic file
    if strcmp(dim, '2d') && (2*M + 1) > K || ...
       strcmp(dim, '3d') && (M + 1)^2 > K
       error(['The requested decode cannot be forfilled as the'...
              ' number of ambisonic channels that would be'...
              ' required exceeds the number of channels in the file'])
    end 
    
%% SET DEFAULT VALUES FOR ALL VARIABLES NOT YET SPECIFIED

    % Set default degree if none selected
    if not(exist('M', 'var'))
        M = 1;
        disp('Default degree being used (1)');
    end
    
    % Set default dimension if none selected
    if not(exist('dim', 'var'))
        dim = '3d';
        disp('Default dimension being used (3D)');
    end
    
    % Set default normalisation if none selected
    if not(exist('norm', 'var'))
        norm = 'n3d';
        disp('Default normalisation being used (N3D)');
    end
    
    % Set default decode strategy if none selected
    if not(exist('decode', 'var'))
        decode = 'pinv';
        disp('Default decode strategy being used (Pseudo-Inverse)');
    end
    
    % Set default decode weighitng if none selected
    if not(exist('weighting', 'var'))
        weighting = 'basic';
        disp('Default decode weighitng being used (basic)');
    end
    
    % Set default configuration if none selected
    if not(exist('config', 'var'))
        config = sprintf('default%ddegree%s', M, dim);
        disp('Default speaker configuration used');
    end
    
    % Set default normTo1 if none selected
    if not(exist('normTo1', 'var'))
        normTo1 = 1;
        disp('Default Binaural amp. normalisation used (normalise)');
    end
    
    % Set default supression if none selected
    if not(exist('supress', 'var'))
        supress = 0;
        disp('Default supression used (dont supress)');
    end
    
%% RESET NUMBER OF AMBISONIC CHANNELS TO USE BASED ON DEGREE SELECTION

    if strcmp(dim, '2d')
        K = (2 * M + 1);
    elseif strcmp(dim, '3d')
        K = (M + 1)^2;
    end
    
end