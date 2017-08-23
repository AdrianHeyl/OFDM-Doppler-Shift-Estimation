clc

sig_512 = abs(fft(audioread('../res/original-OFDM-signal/512-15-57-33.wav')));
sig_1024 = abs(fft(audioread('../res/original-OFDM-signal/1024-15-57-16.wav')));
sig_2048 = abs(fft(audioread('../res/original-OFDM-signal/2048-15-57-01.wav')));

figure;
subplot(3, 1, 1)
plot(sig_512(1:length(sig_512) / 2))
title('512 subcarriers')

subplot(3, 1, 2)
plot(sig_1024(1:length(sig_1024) / 2))
title('1024 subcarriers')

subplot(3, 1, 3)
plot(sig_2048(1:length(sig_2048) / 2))
title('2048 subcarriers')