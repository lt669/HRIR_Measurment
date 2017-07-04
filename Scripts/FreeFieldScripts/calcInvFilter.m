%% Script to calculate the inverse filter of an IR
function [inverseFilter] = calcInvFilter(irArray)

    %{  
        Get size of irArray
        p = number of audio files
        m = length of audio file
        n = number of channels
    %}
    [p,m,n] = size(irArray);
    % Create array of inverse filters
    inverseFilter = zeros(p,4096,n);
    FR = zeros(4096,n);
    
    % Calculate inverse filter for each IR
    for k=1:p
        FR(:,:) = fft(irArray(k,:,:),4096);
        disp(sprintf('size(FR) = %d',size(FR)));
        disp(sprintf('irArray length = %d',length(irArray)));
        inverseFilter(k,:,:) = 1./FR(:,:);
    end
end