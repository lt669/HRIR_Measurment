function varargout     = asxEncode(source, phi, theta, varargin)   
%
% Version 1.0 Dec 2016
%
% [W, Y, Z, X, V, ...] = f(audio source, azimuth, elevation, options)
%
% varargin:   Angular unit: Radians    (rad)          - default
%                           Degrees    (deg) 
%
%            Normalisation: N3D        (n3d)          - default
%                           SN3D       (sn3d)
%                           Orthonomal (orthonormal)
%
% Script encodes a 1D vector (source) into ACN ordered ambisonic 
% channels (varargout) for a given angle (phi, theta) with options 
% (varargin). One of three normalisation schemes can be chosen and 
% angles can be given in either radians or degrees.
%
% The degree of ambisonics by which the input vector is encoded into is
% defined by the number of outputs specified when the function is
% called. The function will encode the vector into as many channels as
% there are outputs in standard ACN ordering e.g.
% [W, Y, Z, X, V, T, R, S, U] = ... would encode up to 2nd degree.
% Whilst it would be unusual, it is possible to encode up to 
% half-degrees e.g. [w, Y, Z, X, V, T] = ...
%
% An ambisonic co-ordinate system is used where phi is azimuth; 0degrees 
% is straight ahead; positive is anticlockwise. Theta is elevation; 
% 0degrees is on the horizontal axis; positive is upwards.
%
% Written by Calum Armstrong, Department of Electronics, The University 
% of York.
%
% Updates ---
% 05 / 12 / 2016: Initial coding
% 13 / 12 / 2016: Finalisation of version 1.0

    % Declare output cell array the same size as number of outputs
    varargout = cell(1,nargout);
    
    % Set default normalisation
    norm = 'n3d';
    
    % Register and handle input options
    for k = 1:length(varargin) % For each input option given...
       switch varargin{k}
           case 'deg' % Angular units in degrees
               phi = phi * pi / 180;
               theta = theta * pi / 180;
           case 'rad' % Angular units in radians
               % Do Nothing
           case 'n3d' % N3D normalisation
               norm = 'n3d';
           case 'sn3d' % SN3D normalisation
               norm = 'sn3d';
           case 'orthonormal' % Orthnormal normalisation
               norm = 'orthonormal';
           otherwise % Unexpected input
               error('Invalid options')
       end
    end
    
    % Cycle through ACN ordered spherical harmonics up to specified 
    % degree and index and calculate the normalised coefficient for the 
    % user-selected angle. Sepperately multiply the origional source by 
    % each coefficient to encode onto each relevant channel.
    [L, M] = deal(0); % Initiate degree and index
    for k = 1:nargout; % For each output channel required...
        
        [~, ~, coef] = iasxShHarmCoef(L, M, phi, theta, norm);
        varargout{k} = coef * source;   
                     % Multiply source by spherical harmonic coefficient
        
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
    
end