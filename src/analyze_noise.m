clc

FS = 44100;

a_stp1_half = '../res/adrian/halfspeednoise-side.wav';
a_stp1_full = '../res/adrian/fullspeednoise-side.wav';
a_stp2_half = '../res/adrian/stp2-halfspeednoise17-24-02.wav';
a_stp2_full = '../res/adrian/stp2-fullspeednoise17-24-46.wav';
m_stp1_half = '../res/mehmedali/half-speed-noise-17-07-28.wav';
m_stp1_full = '../res/mehmedali/full-speed-noise-17-08-17.wav';
m_stp2_half = '../res/mehmedali/half-speed-noise-17-24-27.wav';
m_stp2_full = '../res/mehmedali/full-speed-noise-17-25-12.wav';

a_noises = [abs(fft(audioread(a_stp1_half))), ...
	abs(fft(audioread(a_stp1_full))), ...
    abs(fft(audioread(a_stp2_half))), ...
	abs(fft(audioread(a_stp2_full)))];

m_noises = [abs(fft(audioread(m_stp1_half))), ...
    abs(fft(audioread(m_stp1_full))), ...
	abs(fft(audioread(m_stp2_half))), ...
    abs(fft(audioread(m_stp2_full)))];

figure;

for i = 1 : 4
    subplot(4, 2, mod(i - 1, 4) + i)
    plot(.1 : .1 : length(a_noises) / 20, a_noises(1:length(a_noises) / 2, i))
    axis([0 FS 0 1.3e3])
    
    subplot(4, 2, mod(i - 1, 4) + i + 1)
    plot(.1 : .1 : length(m_noises) / 20, m_noises(1:length(m_noises) / 2, i))
    axis([0 FS 0 .1e3])
end