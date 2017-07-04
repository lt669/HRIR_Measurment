function [filtW, filtXYZ] = spatialCompensation(r, Fs)

% Define Variable
nfft = 1024; 
c = 343; % Speed of Sound
i = 1;

% Calculate Complex Frequency Response of desired filters
% (note that only the first half of these values are correct! the second
% half are replaced below - see ***)
for f = 0:Fs/nfft:Fs-(Fs/nfft)
    
    w = 2*pi*f;
    
    % Compute W compensation response

    Num = 1 + (j*w*r)/c - (1/3)*(((w*r)/c)*((w*r)/c));
    Den = 1 + (1/3)*(j*w*r)/c;
    FW(i) = Num/Den;
    
    % Compute XYZ compensation response
    
    Num = sqrt(6)*(1 + (1/3)*(j*w*r)/c - (1/3)*(((w*r)/c)*((w*r)/c)));
    Den = 1 + (1/3)*(j*w*r)/c;
    FXYZ(i) = Num/Den;
        
    i = i+1;
end

% (***) Recalculate second half of complex frequency response from the
% complex conjugate of the first half.
FW(1, (nfft/2)+2:nfft) = conj(fliplr(FW(1, 2:nfft/2)));
FXYZ(1, (nfft/2)+2:nfft) = conj(fliplr(FXYZ(1, 2:nfft/2)));

% Set middle value (at nfft / 2 + 1) to 0. This makes taking the IFFT much
% easier and doesn't effect perceptual response.
FW(1, nfft/2+1) = 0;
FXYZ(1, nfft/2+1) = 0;

% Take the IFFT
timeFW = real(ifft(FW));
timeFXYZ = real(ifft(FXYZ));

% 'Rotate' / shift impulse responcse to give a single peak. This gives you
% the FIR filters for the W channel and XYZ channels.
filtW = circshift(timeFW, nfft/2, 2);
filtXYZ = circshift(timeFXYZ, nfft/2, 2);

end

