function [ C ] = iasxC( phi, theta, M, dim, norm )

% Version 1.0 Mar 2017
%
% [ Decoding Matrix,
%   Re-Encoding Matrix ] From f( Azimuths,
%                                Elevations,
%                                Decode strategy,
%                                Degree,
%                                Dimension,
%                                Normalisation )
%
% Script computes a re-ecoding matrix 'c'. Matrix is calculated for a 
% loudspeaker array with loudspeaker locations ('phi', 'theta'). The 
% 'Ambisonic Order' / degree and dimention (2D / 3D) of the 
% matrix are given by 'M' and 'dim' respectively. Normalisation of 
% 'norm' is used in the calculation of the Spherical Harmonics.
%
% Written by Calum Armstrong, Department of Electronics, The University 
% of York.
%
% Updates ---
% 20 / 03 / 2017: Initial coding

%% DECLARATIONS

    if strcmp(dim, '2d')
        nChannels = (2 * M + 1);
    elseif strcmp(dim, '3d')
        nChannels = (M + 1)^2;
    end

    C = zeros(nChannels, length(phi));

%% CALCULATE RE-ENCODING MATRIX

    for m = 0:M % For each degree
        for i = -m:m % For each index
            
            % Calculate Spherical Harmonics coefficients
            [~, ~, rMesh] = iasxShHarmCoef(m, i, phi, theta, norm);
            
            % Extract relavent Spherical Harmonic coefficients from mesh
            % matrix
            idx = sub2ind(size(rMesh), 1:size(rMesh, 1), ...
                                                      1:size(rMesh, 1));
            r = rMesh(idx);
            
            % Add co-efficients to re-encoding matrix C
            ACN = (m*(m+1)+i+1);
            C(ACN, :) =  double(r);
            
        end
    end
    
end