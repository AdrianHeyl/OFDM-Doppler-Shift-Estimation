%generate real OFDM signal with 64 subcarriers, using 128 bin ifft, cp = 32
%both prefix and postfix
% y = IFFT(x), the first component is DC, the (N/2+1) th component is the highest frequency, 
% the following (N/2-1) components are negtive frequency 
clear all, close all, clc;

%% -----------generate cplx modulated ofdm signal, proof it ---------------


%% setting
t_duration = 20;  %20s
load('golaySeq.mat');

N_sc = 512;  %number of carriers
switch N_sc
    case 256
        step_sc = 8;
        golay = G32pair(:,1);
    case 512
        step_sc = 16;
        golay = G32pair(:,1);
    case 1024
        step_sc = 32;
        golay = G32pair(:,1);
    case 2048
        step_sc = 64;
        golay = G32pair(:,1);
    otherwise
        disp('Error Nsc!');
end

Fs = 44100;

bw_sc = Fs/N_sc/2;    %bandwidth of each subcarrier
ifft_size = 2*N_sc;
bidata_len = N_sc;
cp_length = 1000;  % the pre and post fix are the same length of ifft_size
symbolCP_len = ifft_size + cp_length;    % len with prefix and postfix
blank_len = 100;
N_symbol = 50;     % number of symbol in a frame
N_frame = ceil(t_duration*Fs/(N_symbol*symbolCP_len));
t = [1:symbolCP_len]'/Fs;


f_min = 10000;
f_max = 14000;
pre_len = 2048;
t_prehalf = [1:pre_len/2]/Fs;
t_lasthalf = [pre_len/2+1:pre_len]/Fs;
preamble = [chirp(t_prehalf,f_min,pre_len/2/Fs,f_max), chirp(t_lasthalf,f_max,pre_len/Fs,f_min)]';


sc_mask = zeros(N_sc,1);   % subcarrier mask

sc_active = [5:step_sc:N_sc];    % active subcarrier index
sc_mask(sc_active) = 1;

%% Modulator, BPPSK modulation, ifft modulation
hMod = comm.BPSKModulator;    % creating bpsk modulator system object
hMod.PhaseOffset = pi/16;     % phase set to pi/16

binary_data = ones(N_sc,1);
binary_data(sc_active) = golay;
BPSK_data = step(hMod,binary_data).*sc_mask;   % this is the bpsk modulated data
% BPSK_data(sc_active) = BPSK_data(sc_active) .* exp(j*2*pi*[1:length(sc_active)]/length(sc_active))';


% build up conjugate symmetry spectrum of Tx for ifft
carriers = [2:N_sc+1]';
NegCarriers = ifft_size - carriers + 2;	%find the bins for the negative frequency carriers
TxSpectrums = zeros(ifft_size,1);
TxSpectrums(carriers) = BPSK_data;
TxSpectrums(NegCarriers) = conj(BPSK_data);
% ifft modulation adding CP
symbol = ifft(TxSpectrums);
symbol_CP = [symbol(end - cp_length + 1 : end); symbol]; % add cyclic prefix

% scale the symbol to 0.8
scale = 1;
symbol_CP = symbol_CP/max(abs(symbol_CP))*scale;

% figure;
% plot(symbol_CP);
%% add preamble, make a frame

frame = [preamble;zeros(blank_len,1);repmat(symbol_CP,N_symbol,1);zeros(blank_len,1)];
seq = repmat(frame,N_frame,1);

% figure;
% plot(seq);

fileName = ['../res/probing-signal/OFDM_',num2str(pre_len),'preamble10k-14k_real',num2str(N_sc),'sc_Nsym',num2str(N_symbol),'_BPSK_cp',num2str(cp_length),'_golay32.wav'];
audiowrite(fileName,seq,Fs);

