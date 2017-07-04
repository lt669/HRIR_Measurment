function [ D, C ] = iasxD( phi, theta, decode, M, dim, norm )

% Version 1.1 Mar 2017
%
% [ Decoding Matrix,
%   Re-Encoding Matrix ] From f( Azimuths,
%                                Elevations,
%                                Decode strategy,
%                                Degree,
%                                Dimension,
%                                Normalisation )
%
% Script computes a Decoding matrix 'D' and re-encoding matrix 'C'. 
% Matricies are calculated for a loudspeaker array with loudspeaker 
% locations ('phi', 'theta') and with decoding strategy 'decode'. The 
% 'Ambisonic Order' / degree and dimention (2D / 3D) of the Decode 
% matrix are given by 'M' and 'dim' respectively. Normalisation of 
% 'norm' is used in the calculation of the Spherical Harmonics.
%
% Written by Calum Armstrong, Department of Electronics, The University 
% of York.
%
% Updates ---
% 15 / 02 / 2017: Initial coding
% 20 / 03 / 2017: Seperated calculation of re-encoding matrix C to new
% function

%% CALCULATE DECODING MATRIX

    switch decode
        case 'proj' % Projection Decode
            C = iasxC(phi, theta, M, dim, norm);
            D = (1/length(phi)) * C';
            
        case 'pinv' % Pseudo-Inverse Decode
            % D = inv(C.' * C) * C.';
            C = iasxC(phi, theta, M, dim, norm);
            D = pinv(C);
            
        case 'engpres' % Energy Preserving Decode
            % http://www.cs.toronto.edu/~jepson/csc420/notes/introSVD.pdf
            C = iasxC(phi, theta, M, dim, norm);
            [U,~,V] = svd(C, 'econ');
            D = sqrt(1/length(phi)) * V * U';
            
        case 'allrad' % All-Round Ambisonic Panning
            % <https://github.com/polarch/Vector-Base-Amplitude-Panning>
            % <https://github.com/polarch/Spherical-Harmonic-Transform>
            
            % Calculate minimum t-design for regular sampling
            t = 2*M + 1; 
            [~, J_angles] = getTdesign(t);
            J_angles_deg = J_angles*180/pi;
            
            C = iasxC(J_angles(:,1), J_angles(:,2), M, dim, norm);
            
            L_angles = [phi',theta'];
            L_angles_deg = L_angles*180/pi;
            A = vbap(J_angles_deg,...
                     findLsTriplets(L_angles_deg),...
                     invertLsMtx(L_angles_deg, findLsTriplets(L_angles_deg))).';
            
            D = A * (1/size(J_angles, 1)) * C';
    end

end