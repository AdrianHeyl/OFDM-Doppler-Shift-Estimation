function y_CFO=add_CFO_audio(y,CFO,Nfft)
    % add CFO (carrier frequency offset) to real audio signal with conjugate
    % symmetric spectrum
    % y : Received signal
    % CFO = IFO (integral CFO) + FFO (fractional CFO)
    % Nfft = FFT size
    flag_debug = 0;

    if isrow(y)
        nn=0:length(y)-1;
    else
        nn= [0:length(y)-1]';
    end

    y_temp = y.*exp(j*2*pi*CFO*nn/Nfft); % Eq.(5.7)
    % if y is real signal with conjugate sysmetric spectrum, we need to
    % reconstruct its spectrum accordingly
    N_sc = Nfft/2;
    carriers = [2:length(y)/2]';
    NegCarriers = length(y) - carriers + 2;

    F_temp = fft(y_temp);
    F_temp_pos = F_temp(carriers);

    yCFO_Spectrums = zeros(length(y),1);
    yCFO_Spectrums(carriers) = F_temp_pos;
    yCFO_Spectrums(NegCarriers) = conj(F_temp_pos);
    y_CFO = ifft(yCFO_Spectrums);

    if flag_debug
        Fs = 44100;
        figure;
        subplot(321);
        plot(y);
        title 'symbol with CP'
        subplot(322);
        stem(linspace(-Fs/2,Fs/2,length(y)),fftshift(abs(fft(y))));
        subplot(323);
        plot(y_temp);
        title 'y temp'
        subplot(324);
        stem(linspace(-Fs/2,Fs/2,length(y)),fftshift(abs(fft(y_temp))));
        subplot(325);
        plot(y_CFO);
        title 'symbol with CFO'
        subplot(326);
        stem(linspace(-Fs/2,Fs/2,length(y_CFO)),fftshift(abs(fft(y_CFO))));
        
    end

end