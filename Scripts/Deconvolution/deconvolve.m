function [ IR ] = deconvolve( inverse_filter, recording)
% DECONVOLVES Calculates the impulse response of a system using its output and
% the inverse filter of the excitation signal

    % Find out how many channels there are in the recording
    channels = min(size(recording));

    % Ensure we're dealing with columnar recorded data
    if channels ~= size(recording,2)      
        recording = recording';
    end
    
    % Ensure we're dealing with columnar filter data
    if size(inverse_filter,2) > size(inverse_filter,1)
        inverse_filter = inverse_filter';
    end
    
    % Create filter data with the right number of channels
    inverse_filter = repmat(inverse_filter(:,1),1,channels);

    % Perform deconvolution
    N = length(recording) + length(inverse_filter);
    IR = ifft( (fft(recording,N)).*fft(inverse_filter,N));

    % Normalise across all channels (REMOVED FOR MASS NORMALISATION)
    IR = IR/max(max(abs(IR)));
    
    % Trim the data (remove silence equal to the sweep - or filter - length)
    IR = IR( length(inverse_filter)+1:end,:);

end

