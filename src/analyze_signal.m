clc

Fs = 44100;

sig_512 = audioread('../res/original-OFDM-signal/512-15-57-33.wav');
sig_1024 = audioread('../res/original-OFDM-signal/1024-15-57-16.wav');
sig_2048 = audioread('../res/original-OFDM-signal/2048-15-57-01.wav');

db_512 = linspace(-Fs/2, Fs/2, length(sig_512));
db_1024 = linspace(-Fs/2, Fs/2, length(sig_1024));
db_2048 = linspace(-Fs/2, Fs/2, length(sig_2048));

figure;
subplot(3, 1, 1)
plot(db_512, fftshift(10 * log10(abs(fft(sig_512)))))
title('512 subcarriers')
xlabel('frequency [Hz]')
ylabel('snr [dB]')

subplot(3, 1, 2)
plot(db_1024, fftshift(10 * log10(abs(fft(sig_1024)))))
title('1024 subcarriers')
xlabel('frequency [Hz]')
ylabel('snr [dB]')

subplot(3, 1, 3)
plot(db_2048, fftshift(10 * log10(abs(fft(sig_2048)))))
title('2048 subcarriers')
xlabel('frequency [Hz]')
ylabel('snr [dB]')