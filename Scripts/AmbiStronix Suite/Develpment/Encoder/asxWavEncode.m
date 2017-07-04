function  varargout = asxWavEncode(sourceFile,FS, phi, theta,...
                                           degree, saveFolder, varargin)
%
% Version 1.0 Dec 2016
%
% Ambisonic File Created From f(audio source, azimuth, elevation, order 
%                               of ambisonics required, options)
%
% Options:       Angular unit: Radians    (rad)          - default
%                              Degrees    (deg) 
%
%               Normalisation: N3D        (n3d)          - default
%                              SN3D       (sn3d)
%                              Orthonomal (orthonormal)
%
%           Output Supression: Supress    (supOut)
%
%         Dont normalise to 1: Don't norm.(raw)
%
% Script simplifies the process of encoding a mono audio file into
% ambisonic format. In standard operation, the function automatically 
% renders a multi-channel ambisonic .wav file encoding input source
% audio (sourceFile) in the spatial location (phi, theta). Audio is
% encoded in degree (degree) with options (varargin). One of three 
% normalisation schemes can be chosen and angles can be given in either 
% radians or degrees. If a multichannel / stereo file is input, the 
% first / Left channel is selected as the input. 
%
% If required, the function can return the matrix form of the rendered
% .wav file along with the associated sampling frequency of the source.
% Further, the saving of the .wav file can be supressed via the 
% output supression option. This is useful when this function is used by
% the parent function ambiStronixMultiWavEncode.m.
%
% An ambisonic co-ordinate system is used where phi is azimuth; 0degrees 
% is straight ahead; positive is anticlockwise. Theta is elevation; 
% 0degrees is on the horizontal axis; positive is upwards.
%
% Written by Calum Armstrong, Department of Electronics, The University 
% of York.
%
% Updates ---
% 08 / 12 / 2016: Initial coding
% 13 / 12 / 2016: Finalisation of version 1.0

    % Set default normalisation
    norm = 'n3d'; 
    % Set default output supress option (used by parent function:
    % ambiStronixMultiWavEncode.m)
    supress = 0;
    normTo1 = 1;

    % Register and handle relavant input options for current function
    k = 1;
    while k <= length(varargin) % For each input option given...
       switch varargin{k}
           case {'deg', 'rad'}
               k = k + 1;
           case 'n3d' % N3D normalisation
               norm = 'n3d';
               k = k + 1;
           case 'sn3d' % SN3D normalisation
               norm = 'sn3d';
               k = k + 1;
           case 'orthonormal' % Orthnormal normalisation
               norm = 'orthonormal';
               k = k + 1;
           case 'supOut' % Supress output
               supress = 1;
               varargin(k) = [];
           case 'raw'
               normTo1 = 0;
               varargin(k) = [];
           otherwise
               error('unexpected options');
       end
    end

    % Read in audio source
    %[source, FS] = audioread(sourceFile);
    source = sourceFile;
    source = source(:, 1); % Make mono if not already
    
    %[~, sourceFile, ~] = fileparts(sourceFile);

    % Depending of degree on encode required...
    switch degree
        case 'bFormat' % B-Format encode
            % Run encode function with desired number of outputs
            [W, Y, Z, X] = asxEncode(source, phi, theta,...
                                             varargin{:});

            % Concatinate outputs into single variable
            output = [W, X, Y, Z];

            % Normalise output
            if normTo1
                outputNorm = output / max(max(abs(output)));
            else
                outputNorm = output;
            end
            
            % Output as .wav file unless supressed
            if supress ~= 1
                filename = sprintf('%s/Ambi_%s_%s_%d_%d_%s.wav',...
                                   saveFolder,...
                                   sourceFile,...
                                   degree,...
                                   phi,...
                                   theta,...
                                   norm); % Construct filename
                audiowrite(filename, outputNorm, FS); % Save file
            end
            
            % Return matrix / sampling rate data
            varargout{1} = outputNorm;
            varargout{2} = FS;

        case 1 % 1st degree encode
            [W, Y, Z, X] = asxEncode(source, phi, theta,...
                           varargin{:});

            output = [   W,...
                      Y, Z, X];

            if normTo1
                outputNorm = output / max(max(abs(output)));
            else
                outputNorm = output;
            end
            
            if supress ~= 1
                filename = sprintf('%s/Ambi_%s_%d_%d_%d_%s.wav',...
                                   saveFolder,...
                                   sourceFile,...
                                   degree,...
                                   phi,...
                                   theta,...
                                   norm);
                audiowrite(filename, outputNorm, FS);
            end
            varargout{1} = outputNorm;
            varargout{2} = FS;

        case 2 % 2nd degree encode
            [W, Y, Z, X, V, T, R, S, U] = ... 
                            asxEncode(source, phi, theta,...
                            varargin{:});

            output = [      W,... 
                         Y, Z, X,... 
                      V, T, R, S, U];

            if normTo1
                outputNorm = output / max(max(abs(output)));
            else
                outputNorm = output;
            end
            
            if supress ~= 1
                filename = sprintf('%s/Ambi_%s_%d_%d_%d_%s.wav',...
                                   saveFolder,...
                                   sourceFile,...
                                   degree,...
                                   phi,...
                                   theta,...
                                   norm);
                audiowrite(filename, outputNorm, FS);
            end
            varargout{1} = outputNorm;
            varargout{2} = FS;

        case 3 % 3rd degree encode
            [W, Y, Z, X, V, T, R, S, U, Q, O, M, K, L, N, P] = ...
                            asxEncode(source, phi, theta,...
                            varargin{:});

            output = [         W,... 
                            Y, Z, X,... 
                         V, T, R, S, U,... 
                      Q, O, M, K, L, N, P];

            if normTo1
                outputNorm = output / max(max(abs(output)));
            else
                outputNorm = output;
            end
            
            if supress ~= 1
                filename = sprintf('%s/Ambi_%s_%d_%d_%d_%s.wav',...
                                   saveFolder,...
                                   sourceFile,...
                                   degree,...
                                   phi,...
                                   theta,...
                                   norm);
                audiowrite(filename, outputNorm, FS);
            end
            varargout{1} = outputNorm;
            varargout{2} = FS;

        case 4 % 4th degree encode
            [W, Y, Z, X, V, T, R, S, U, Q, O, M, K, L, N, P, ...
                       a4, b4, c4, d4, e4, f4, g4, h4, i4] = ...
                            asxEncode(source, phi, theta,...
                            varargin{:});

            output = [                 W,... 
                                   Y,  Z,  X,... 
                               V,  T,  R,  S,  U,... 
                           Q,  O,  M,  K,  L,  N,  P,... 
                      a4, b4, c4, d4, e4, f4, g4, h4, i4];

            if normTo1
                outputNorm = output / max(max(abs(output)));
            else
                outputNorm = output;
            end
            
            if supress ~= 1
                filename = sprintf('%s/Ambi_%s_%d_%d_%d_%s.wav',...
                                   saveFolder,...
                                   sourceFile,...
                                   degree,...
                                   phi,...
                                   theta,...
                                   norm);
                audiowrite(filename, outputNorm, FS);
            end
            varargout{1} = outputNorm;
            varargout{2} = FS;

        case 5 % 5th degree encode
            [W, Y, Z, X, V, T, R, S, U, Q, O, M, K, L, N, P, ...
                         a4, b4, c4, d4, e4, f4, g4, h4, i4, ...
                 a5, b5, c5, d5, e5, f5, g5, h5, i5, j5, k5] = ...
                            asxEncode(source, phi, theta,...
                            varargin{:});

            output = [                     W,... 
                                       Y,  Z,  X,... 
                                   V,  T,  R,  S,  U,... 
                               Q,  O,  M,  K,  L,  N,  P,... 
                          a4, b4, c4, d4, e4, f4, g4, h4, i4,... 
                      a5, b5, c5, d5, e5, f5, g5, h5, i5, j5, k5];

            if normTo1
                outputNorm = output / max(max(abs(output)));
            else
                outputNorm = output;
            end
            
            if supress ~= 1
                filename = sprintf('%s/Ambi_%s_%d_%d_%d_%s.wav',...
                                   saveFolder,...
                                   sourceFile,...
                                   degree,...
                                   phi,...
                                   theta,...
                                   norm);
                audiowrite(filename, outputNorm, FS);
            end
            varargout{1} = outputNorm;
            varargout{2} = FS;
        
    end