function [Sigout, nChan] = GK_Convolution(h, r, predelayflg)        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function convolves signal r with filter h. Parameters are:
% h: filter for the deconvolution
% r: Signal to be filtered (Can be multichannel)
% Sigout: The output multichannel signal
% nChan: Number of channels in output IR (same as recorded input)
% predelayflg: If this parameter is set to 1, the response is truncated so
% the start occurs at a time after the length of the inverse filter

% Adapted by Lewis Thresh 29/12/2016 to be used with multichannel
% inverse filters
%%%%%%%%%%%%%% Gavin Kearney 2015, University of York %%%%%%%%%%%%%%%%%%%%%

nChan = size(r,2); % Number of channels in recorded file

lenr = length(r);
lenfilt = length(h);
nfft = (lenr + lenfilt);

% For debugging
%{
[n,m,l] = size(r);
disp(sprintf('GK: lenr=%d',lenr));
disp(sprintf('GK: n=%d m=%d l=%d',n,m,l));
[x,y,z] = size(h);
disp(sprintf('GK: lenfilt=%d',lenfilt));
disp(sprintf('GK: x=%d y=%d z=%d',x,y,z));
disp(sprintf('GK: nfft=%d',nfft));
%}

Sigout = zeros(nfft, nChan);


    x1=[h; zeros(lenr, 2)]; %zero pad the input vectors to avoid circular convolution

    x2=[r(:,:); zeros(lenfilt, 2)];

    X1=fft(x1, nfft); % Get Fourier transform of X1

    X2=fft(x2, nfft); % Get Fourier transform of X2

    Y=X1.*X2; %frequency domain multiplication performs time domain convolution

    %Y(round(nfft1/2)) = 0;
    Sigout(:,:)=real(ifft(Y, nfft)); % Get output impulse responses

if predelayflg ==1
    Sigout =  Sigout(lenfilt+1:end,:);
    %disp('Truncating Impulse Response..');
else
    return
end




