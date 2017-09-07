clear all;
close all;
clc;

%% settings

duration = 20;
destination = '../res/exp07-09-2017/setup2-sc256-No1.wav';

Fs = 44100;

%% record audio

disp('start recording...')

rec = audiorecorder(Fs, 16, 1);
recordblocking(rec, duration);

disp('end recording.')

%% write file

audiowrite(destination, getaudiodata(rec), Fs);