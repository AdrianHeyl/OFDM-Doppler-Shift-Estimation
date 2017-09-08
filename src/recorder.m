clear all;
close all;
clc;

%% settings

duration = 24;
destination = '../res/exp08-09-2017/setup2-sc512-No4.wav';

Fs = 44100;

%% record audio     

disp('start recording...')

rec = audiorecorder(Fs, 16, 1);
recordblocking(rec, duration);

disp('end recording.')

%% write file

audiowrite(destination, getaudiodata(rec), Fs);