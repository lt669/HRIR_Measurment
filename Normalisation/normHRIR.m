%{
    Script for normalising all of the deconvolved sweeps with each other
%}

function [normalised] = normHRIR(input)

    [n,m,p] = size(input);
    disp('--- normHRIR ---');
    disp(sprintf('Size: n=%d m=%d p=%d',n,m,p));
    
    % Find maximum value
    maximum = max(max(max(abs(input))));
    
    % Apply Normalisation
    normalised = input/maximum;
    
    [n,m,p] = size(normalised);
    disp(sprintf('Normalised Size: n=%d m=%d p=%d',n,m,p));
    
    maximumNorm = max(max(max(abs(normalised))));
    
    disp(sprintf('maximumNorm: %d',maximumNorm));
    disp('--- normHRIR Exit---');
end