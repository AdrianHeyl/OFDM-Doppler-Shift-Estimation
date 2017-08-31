clc

Fs = 44100;
filepath_1 = '../res/adrian/';
filepath_2 = '512-halfspeed.wav';

source = [filepath_1, 'unfiltered/', filepath_2];
destination = [filepath_1, 'filtered/', filepath_2];

x = audioread(source);
t = linspace(1, length(x) / Fs, length(x));
db = linspace(-Fs/2,Fs/2,length(x));

figure;
subplot(4, 1, 1)
plot(t, x)

y = fft(x);

subplot(4, 1, 2)
plot(db, fftshift(10 * log10(abs(y))))
%plot(db, fftshift(10 * log10(abs(fft(x)))))
axis([-2.5e4 2.5e4 -50 50])

filter_512_start = [1; 1081999 - 882000; 1241999 - 882000; 1391981 - 882000];
filter_512_stop = [1001999-882000; 1161999 - 882000; 1321981 - 882000; 1764000 - 882000];

for i = 1 : length(filter_512_start)
%     for j = filter_512_start(i) : filter_512_stop(i)
%         y(j) = 0;
%         y(length(y) - j + 1) = 0;
%     end
    y(filter_512_start(i):filter_512_stop(i)) = 0;
    y(length(y) - filter_512_stop(i) + 1:length(y) - filter_512_start(i) + 1) = 0;
end

subplot(4, 1, 3)
plot(db, fftshift(10 * log10(abs(y))))
axis([-2.5e4 2.5e4 -50 50])

z = ifft(y);

subplot(4, 1, 4)
plot(t, z)