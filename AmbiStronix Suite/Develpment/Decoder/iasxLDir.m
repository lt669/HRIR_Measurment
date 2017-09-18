function [ phi, theta, phiDeg, thetaDeg ] = iasxLDir( configuration )

% Version 1.0 Feb 2017
%
% [ Azimuths, 
%   Elevations
%   Azimuths (degrees), 
%   Elevations (degrees) ] From f( loudspeaker Configuration )
%
% Script provides a convinient way to 'store' sets of loudspeaker
% directions 'saved' under the name of their configuration. It provides 
% a quick and easy way to retriev this information from withtin other 
% scripts. New custom loudspeaker configurations can be easily added as
% shown below.
% 
% Angles can be stored in radians (phi, theta) or in degrees (phiDeg,
% thetaDeg) within a case statement. The name of the configuration 
% should be put as the 'case'. Angles should be stored as 2 separate 
% arrays where matching indexes represent a single angle. An example is 
% shown for how to store a '7-design' configuration in degrees.
%
% Written by Calum Armstrong, Department of Electronics, The University 
% of York.
%
% Updates ---
% 15 / 02 / 2017: Initial coding

%% LOUDSPEAKER CONFIGURATIONS

    switch configuration
    % CUSTOM CONFIGURATIONS
    
        %!%! INSERT YOUR CONFIGURATIONS HERE !%!%
    
        case {'octohedron2'}
            phiDeg = [0 180 90 270 0 180];
            thetaDeg = [45 45 0 0 -45 -45];
            
        case '7-design' % EXAMPLE
            phiDeg = [26.0011675216559,333.998832478344,17.1086452559122,342.891354744088,153.998832478344,206.001167521656,162.891354744088,197.108645255912,72.8913547440879,107.108645255912,116.001167521656,63.9988324783441,252.891354744088,287.108645255912,296.001167521656,243.998832478344,32.2544599366034,212.254459936603,302.254459936603,122.254459936603,327.745540063397,147.745540063397,57.7455400633966,237.745540063397];
            thetaDeg = [15.4641512961471,-15.4641512961471,-24.9937030546433,24.9937030546433,-15.4641512961471,15.4641512961471,24.9937030546433,-24.9937030546433,24.9937030546433,-24.9937030546433,15.4641512961471,-15.4641512961471,24.9937030546433,-24.9937030546433,15.4641512961471,-15.4641512961471,60.0253819510733,60.0253819510733,60.0253819510733,60.0253819510733,-60.0253819510733,-60.0253819510733,-60.0253819510733,-60.0253819510733];
     
        case 'tetrahedron'
            phiDeg = [0, 90, 180, 270];
            thetaDeg = [35.26438968, -35.26438968, 35.26438968, -35.26438968];
            
        case 'frontPoint'
            phiDeg = [0, 90, 270, 0];
            thetaDeg = [35.26438968, 35.26438968, 35.26438968, -35.26438968];
            
        case 'test'
            phiDeg = [0, 270, 0, 0];
            thetaDeg = [0, 0, 90, 45];
            
        case 'atmos_9.1.6'
            phiDeg = [0, 26, 60, 100, 142.5, -142.5, -100, -60, -26, 34.596, 90, 145.404, -145.404, -90, -34.596];
            thetaDeg = [0, 0, 0, 0, 0, 0, 0, 0, 0, 39.460, 64, 39.460, 39.460, 64, 39.460];
            
        case 'atmos_9.1.6_rounded'
            phiDeg =   [0, 30, 60, 100, 145, 215, 260, 300, 330, 40, 90, 140, 220, 270, 320];
            thetaDeg = [0,  0,  0,   0,   0,   0,   0,   0,   0, 40, 60,  40,  40,  60,  40];
            
        case 'auro_13.1'
            phiDeg = [0, 30, 110, 142.5, -142.5, -110, -30, 0, 30, 110, -110, -30, 0];
            thetaDeg = [0, 0, 0, 0, 0, 0, 0, 30, 30, 30, 30, 30, 90];
            
        case 'auro_13.1_rounded'
            phiDeg = [0, 30, 110, 145, 215, 250, 330, 0, 30, 110, 250, 330, 0];
            thetaDeg = [0, 0, 0, 0, 0, 0, 0, 30, 30, 30, 30, 30, 90];
            
        case '5.1'
            phiDeg = [0, 26, 100, -100 -26];
            thetaDeg = [0, 0, 0, 0, 0];
            
        case '5.1_rounded'
            phiDeg = [0, 30, 100, 260, 330];
            thetaDeg = [0, 0, 0, 0, 0];
            
        %!%! ------------------------------- !%!%
        
    % DEFAULT CONFIGURATIONS
    % FIRST ORDER
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

    % SECOND ORDER
        case {'hex', 'default2degree2d'}
            phiDeg = [0 60 120 180 240 300];
            thetaDeg = [0 0 0 0 0 0];

        case {'dodecahedron', 'default2degree3d'}
            phiDeg = [180 50 310 118 242 0 180 62 298 130 230 0];
            thetaDeg = [63 46 46 16 16 0 0 -16 -16 -46 -46 -63];

    % THIRD ORDER
        case {'octagon', 'default3degree2d'}
            phiDeg = [0 45 90 135 180 225 270 315];
            thetaDeg = [0 0 0 0 0 0 0 0];

        case '16chSpherePacking'
            phiDeg = [0 90 180 270 45 135 225 315 0 90 180 270 45   ...
                                                          135 225 315];
            thetaDeg = [51 51 51 51 14 14 14 14 -14 -14 -14 -14 -51 ...
                                                          -51 -51 -51];

        case {'26ptLebedev', 'default3degree3d'}
            phiDeg = [0 0 90 180 270 45 135 225 315 0 45 90 135 180 ...
                            225 270 315 45 135 225 315 0 90 180 270 0];
            thetaDeg = [90 45 45 45 45 35 35 35 35 0 0 0 0 0 0 0 0 ...
                                  -35 -35 -35 -35 -45 -45 -45 -45 -90];

    % FOURTH ORDER
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
                                                       
    % FIFTH ORDER
        case {'12chCircular', 'default5degree2d'}
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
            ME = MException('VerifyInput:InvalidInput', ...
                            'Invalid Loudspeaker Configuration');
            throw(ME);
    end
    
    % Convert angles to radians if given in degrees
    if exist('phiDeg', 'var')
        phi = phiDeg * pi / 180;
        theta = thetaDeg * pi / 180;
    else
        phiDeg = phi * 180 / pi;  %#ok<NODEF>
        thetaDeg = theta * 180 / pi; %#ok<NODEF>
    end

end