clear all, close all, clc;

Fpass = 100;
Fstop = 150;
Apass = 1;
Astop = 65;
Fs = 44100;

d = designfilt('highpassfir', 'PassbandFrequency', Fpass, ...
                'StopbandFrequency', Fstop, 'PassbandRipple', Apass, ...
                'StopbandAtenuation', Astop, 'DesignMethod', 'equiripple', ...
                'SampleRate', Fs);