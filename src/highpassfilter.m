clear all, close all, clc;

Fstop = 350;
Fpass = 400;
Astop = 65;
Apass = 0.5;
Fs = 44100;

d = designfilt('highpassfir','StopbandFrequency',Fstop, ...
  'PassbandFrequency',Fpass,'StopbandAttenuation',Astop, ...
  'PassbandRipple',Apass,'SampleRate',Fs,'DesignMethod','equiripple');

fvtool(d)
            
disp('end');