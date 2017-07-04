function varargout = asxWavEncodeIMF(sourceFile, phi, theta, interaural, R,...
                                           degree, saveFolder, varargin)
%
% Version 1.0 Dec 2016
%
%
% An ambisonic co-ordinate system is used where phi is azimuth; 0degrees 
% is straight ahead; positive is anticlockwise. Theta is elevation; 
% 0degrees is on the horizontal axis; positive is upwards.
%
% Written by Calum Armstrong, Department of Electronics, The University 
% of York.
% 
% Updates ---
% 15 / 06 / 2017: Initial coding
    
    X = interaural / 2;
    
    test1 = abs(R) * sin(phi) * cos(theta);
    test2 = cos(phi);


%% LEFT
    
    
    
%   [matrixL, FS] = asxWavEncode(source, phi, theta,...
 %                               degree,...
 %                               savefolder,...
 %                               varargin{:});
               
%% Right
    
    A_ = sqrt((R*cos(theta))^2 + X^2 - 2*X*R*cos(theta)*cos(pi/2+phi)) ;
    phi_ = pi/2 - asin((R*cos(theta)*sin(pi/2+phi))/(A_));
    R_ = sqrt((R*sin(theta))^2 + A_^2);
    theta_ = asin((R*sin(theta))/(R_));
    
    if test1 < -0.11
        phi_ = - phi_;
    end

%    [matrixR, FS] = asxWavEncode(source, phi_, theta_,...
%                                 degree,...
%                                 savefolder,...
%                                 varargin{:});



end