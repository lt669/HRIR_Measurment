    function h = GK_Freqplot(x,Fs, Nfft, style, mylinewidth, myfontsize, mytitle, myxlabel, myylabel,  my_ylim, my_xlim,  normflag, dcfiltflag);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to plot the single sided frequency spectrum of an audio signal.
% Returns handle to the figure.
% Parameters:
% x: Input signal
% Fs: Sampling rate
% style: Plot style. e.g. 'r:' plots a red dotted line. See 'plot' help
% myfontsize: Text font size
% mytitle: figure title
% myxlabel: x-axis label
% myylabel: y-axis label
% my_ylim, my_xlim: axis limits in Hz.
% normflag: Normalise the spectrum prior to plotting. Makes the peak 0dBFS
% dcfiltflag: Removes any DC offset via a 30Hz highpass filter
%%%%%%%%%%%%%%% Gavin Kearney, University of York, 2015 %%%%%%%%%%%%%%%%%%%



if dcfiltflag == 1
    dcfilt = fir1(Nfft, 30/(Fs/2),'high'); % high pass filter
    x = fftfilt(dcfilt,x); 
end

X = fft(x, Nfft); % Get Fast Fourier transform

if normflag == 1
    X = X./max(abs(X(:))); % Normalize if requested
end

frlow = round(my_xlim(1)*Nfft/Fs); % Compute freq bins for x-axis limits
frhigh = round(my_xlim(2)*Nfft/Fs);


f = Fs/Nfft:Fs/Nfft:Fs; % Frequency vector for plotting

h = semilogx(f(frlow:frhigh),20*log10(abs(X(frlow:frhigh))), style, 'linewidth', mylinewidth); % plot
%h = semilogx(f(frlow:frhigh), style, 'linewidth', mylinewidth); % plot

% Formatting
xlim(my_xlim);
ylim(my_ylim);
title(mytitle, 'FontSize', myfontsize)
xlabel(myxlabel, 'FontSize', myfontsize);
ylabel(myylabel, 'FontSize', myfontsize);
grid on;
set(gca, 'fontsize', myfontsize);
