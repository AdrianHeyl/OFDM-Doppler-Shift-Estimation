clear all;
close all;
clc;

%% settings

duration = 5;
destination = '../res/07-09-2017/file.wav';

Fs = 44100;

%% record audio

disp('start recording...')

rec = audiorecorder(Fs, 16, 1);
recordblocking(rec, duration);

disp('end recording.')

%% write file

audiowrite(destination, getaudiodata(rec), Fs);