function [rL, rR, sL, sR, ITDD, ILDD] = DAFxFileComparison_func(HRTF, compar, degree, phi, theta, Filt)

    % Data
    Fs = 48000;
    Nfft = 2^nextpow2(Fs);
    %f = 0:Fs/Nfft:Fs/2; % one-sided
    
    low_freq = 0;
    hi_freq = Fs/2;
    low_freq_samp = round(Nfft / Fs * low_freq)+1;
    hi_freq_samp  = round(Nfft / Fs *  hi_freq)+1;
    
%%
    % Read in files
    ref  = audioread(sprintf('%s/Bi_HRTF_%d_%d.wav',HRTF, phi,theta));
    test = audioread(sprintf('%s/Bi_test_%d_%d_%d.wav',compar,degree,phi,theta));
    
%     ref = filter(Filt,ref);
%     test = filter(Filt,test);
    
    refL  =  ref(:,1);  refR =  ref(:,2);
    testL = test(:,1); testR = test(:,2);
%%
    % Calculate frequency response
    freqResp_refL  = abs(fft( refL, Nfft));
    freqResp_refR  = abs(fft( refR, Nfft));
    freqResp_testL = abs(fft(testL, Nfft));
    freqResp_testR = abs(fft(testR, Nfft));

%     % Calculate Power
%     freqResp_refL  = freqResp_refL.^2;
%     freqResp_refR  = freqResp_refR.^2;
%     freqResp_testL = freqResp_testL.^2;
%     freqResp_testR = freqResp_testR.^2;

    % Extract relavent of freq. response
    freqResp_refL  = freqResp_refL(low_freq_samp:hi_freq_samp);
    freqResp_refR  = freqResp_refR(low_freq_samp:hi_freq_samp);
    freqResp_testL = freqResp_testL(low_freq_samp:hi_freq_samp);
    freqResp_testR = freqResp_testR(low_freq_samp:hi_freq_samp);

    CCL = corrcoef(freqResp_refL, freqResp_testL);
    rL = CCL(2,1);
    CCR = corrcoef(freqResp_refR, freqResp_testR);
    rR = CCR(2,1);
    sL = rms(testL) / rms(refL);
    sR = rms(testR) / rms(refR);
    
    ILDD = (rms(refR) / rms(refL)) / (rms(testR) / rms(testL));
    maxLag = round(0.0011*Fs);
    ITDD = 1 - (abs(finddelay(refL,refR,maxLag) - finddelay(testL,testR,maxLag)) / maxLag);
    
end