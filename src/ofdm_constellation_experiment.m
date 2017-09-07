% this file analyses the recorded signal, shows that the change of
% constellation of different symbols
clear all;
close all;
clc;

%% setting

N_sc = 1024;  % number of carriers

Fs = 44100;   % sapmling frequency
bw_sc = Fs/N_sc;    %bandwidth of each subcarrier
ifft_size = 2*N_sc;
bidata_len = N_sc;  %binary data length

cp_length = 256;  % the pre and post fix are the same length of ifft_size, hereby I used both cyclic prefix and postfix
symbolCP_len = ifft_size + cp_length;    % symbol length with prefix and postfix
blank_len = 100;    % there is a blank interval between two frames
N_symbol = 100;     % number of symbol in a frame
N_frame = 2;        % number of frames in generated audio file
N_cluster = N_frame;
t = (1:symbolCP_len)'/Fs;
deg = zeros(N_symbol,1);    % angle of each symbol, unit degree

c= 343;
v = 0.3;
flag_MAFlt_perSC = 0;
sync_offset = -128;


sc_mask = zeros(N_sc,1);   % subcarrier mask
switch N_sc
    case 512
        sc_active = [100 200 300];    % active subcarrier index
    case 1024
        sc_active = [100 200 300 400 500];    % active subcarrier index
    case 2048
        sc_active = [100:100:1000];    % active subcarrier index
    otherwise
        disp('wrong Nsc!')     
end
sc_mask(sc_active) = 1;    % only active subcarriers are enabled to transmit data, others are blocked

hMod = comm.BPSKModulator;    % creating bpsk modulator system object
hMod.PhaseOffset = pi/16;     % phase set to pi/16
binary_data = ones(N_sc,1);
BPSK_data = step(hMod,binary_data).*sc_mask;   % this is the bpsk modulated data



filename = '../res/adrian/unfiltered/stp2-1024-fullspeed.wav'; %BH = better hardware
%% preamble

f_min = 8000;
f_max = 12000;
pre_len = 256;
t_prehalf = (1:pre_len/2)/Fs;
t_lasthalf = (pre_len/2+1:pre_len)/Fs;
preamble = [chirp(t_prehalf,f_min,pre_len/2/Fs,f_max),...
    chirp(t_lasthalf,f_max,pre_len/Fs,f_min)]';

len_frame = pre_len + symbolCP_len * N_symbol + blank_len*2;

%% synchronization, find index using passband and matched filter

% read audio file
sig_received = audioread(filename);

% passband filter 7k to 13k.
preamble_filter = firls(127, [0 6500 7000 13000 13500 Fs/2] ./ (Fs / 2), [0 0 1 1 0 0]);
sig_noist_flt = filter(preamble_filter,1,sig_received);

% matched filter
coef_MF_preamble = preamble(end:-1:1);
data_MFflted = filter(coef_MF_preamble,1,sig_noist_flt);

figure;
plot(linspace(1 / Fs, length(sig_received), length(data_MFflted)), data_MFflted);
title('data after matched filter for synchronization');

sync_threshold = 0.4;

index_arr = find(data_MFflted > sync_threshold);
index_arr_sorted = sort_index(data_MFflted,index_arr,N_cluster);
disp(index_arr_sorted);


%% sync offline
index_arr = sync_offline(sig_noist_flt,preamble,len_frame);
disp('index_arr');
disp(index_arr);
disp('diff');
disp(diff(index_arr));

N_targetFrame = length(index_arr);

index_arr = 705782;
%% design of bandpass and lowpass FIR filter

% % lowpass filter
% Fp_LP = 7000/Fs;
% Fs_LP = 8000/Fs;
% order = 127;
% coef_LP = firls(order,[0 Fp_LP Fs_LP 1],[1 1 0 0]);
% 
% %bandpass filter 16kHz ~ 18.4kHz
% Fs1_BP = 300*2/Fs;
% Fp1_BP = 600*2/Fs;
% Fp2_BP = 20750*2/Fs;
% Fs2_BP = 21050*2/Fs;
% order = 127;
% coef_BP = firls(order,[0 Fs1_BP Fp1_BP Fp2_BP Fs2_BP 1],[0 0 1 1 0 0]);

%% demodulation

v_est_diff = [];
v_est_grad = [];
deg_total = [];
deg_unwraped_total = [];
w_symbol_diff = [];
w_symbol_gra = [];

for i = 1:length(index_arr)
    deg = zeros(N_symbol,length(sc_active));
    
    for j = 1:N_symbol
        index_symbol = j;
        i_start = index_arr(i) + blank_len + symbolCP_len*(index_symbol-1) + 1 ;
        i_end = i_start + symbolCP_len - 1;
        target_sym = sig_received(i_start: i_end);

        symbol_woCP = target_sym(cp_length + 1 + sync_offset : end + sync_offset );
        frequency_data = fft(symbol_woCP);
        BPSK_demodulated = frequency_data(2:N_sc);
        % normalization
        BPSK_demodulated = BPSK_demodulated / max(abs(BPSK_demodulated(sc_active)));
        
        % calculate the angle of a symbol
        deg(j,:) = angle(BPSK_demodulated(sc_active));    

        if 0
            figure;
            xlim([-1.2 1.2]);
            ylim([-1.2 1.2]);
            hold on;
            plot(BPSK_data,'ro','MarkerSize',10,'linewidth',1.5);
            plot(BPSK_demodulated,'x','MarkerSize',10,'linewidth',1.5);

            x = real(BPSK_demodulated(sc_active));
            y = imag(BPSK_demodulated(sc_active));
            text(x,y,num2str(sc_active'));

            title (['Doppler exp ,speed=',num2str(v),'m/s, symbol #',num2str(j)]);
%             print(['symbol#',num2str(j)],'-dpng');
        end
    end

    deg_total = [deg_total;deg];
    deg_unwraped = unwrap(deg);
    deg_unwraped_total = [deg_unwraped_total;deg_unwraped];
    
    % angular speed based on phase diff
    diff_temp = diff(deg_unwraped/pi);
    w_symbol_diff = [w_symbol_diff;diff_temp(1,:);diff_temp];
    
    % angular speed based on phase gradient
    [f_temp,w_symbol_temp] = gradient(deg_unwraped/pi);
    w_symbol_gra = [w_symbol_gra;w_symbol_temp];
end

%% speed estimation based on phase information
factor = 2.85;
% factor = 2;
% symbol level speed for the whole sequence 
v_symbol_diff = bsxfun(@rdivide, w_symbol_diff.*c./factor, sc_active);
v_symbol_gra = bsxfun(@rdivide, w_symbol_gra.*c./factor, sc_active);

% apply moving average filter per subcarrier
lag_movavg = 8;
a = 1;
b = ones(1,lag_movavg)/lag_movavg;

% step2 means after the filter or not, the new thing
if flag_MAFlt_perSC
    v_symbol_diff_step2 = filter(b,a,v_symbol_diff);
    v_symbol_diff_step2(1:lag_movavg,:) = v_symbol_diff(1:lag_movavg,:);
    v_symbol_gra_step2 = filter(b,a,v_symbol_gra);
    v_symbol_gra_step2(1:lag_movavg,:) = v_symbol_gra(1:lag_movavg,:);
else
    v_symbol_diff_step2 = v_symbol_diff;
    v_symbol_gra_step2 = v_symbol_gra;
end

% select subcarriers which provide most smoothy curves
% use the variation of gradient of the curves
[fx_diff,fy_diff] = gradient(v_symbol_diff_step2);
[fx_gra,fy_gra] = gradient(v_symbol_gra_step2);
% select 10 subcarriers with least peaks to calculate the avg 
[temp1, ix_diff] = sort(var(fy_diff));
[temp2, ix_gra] = sort(var(fy_gra));


% averaged symbol speed estimation across 10 subcarriers
N_sc_avg = 3;
v_symbol_diff_avg = mean(v_symbol_diff_step2(:,ix_diff(1:N_sc_avg)),2);
v_symbol_gra_avg = mean(v_symbol_gra_step2(:,ix_gra(1:N_sc_avg)),2);
% use the subcarrier with least fluctuation
v_symbol_diff_min = v_symbol_diff_step2(:,ix_diff(1));
v_symbol_gra_min = v_symbol_gra_step2(:,ix_gra(1));


figure;
hold on;
plot(v_symbol_diff_avg);
plot(v_symbol_gra_avg);
plot(v_symbol_diff_min);
plot(v_symbol_gra_min);
legend('diff-avg','gra-avg','diff_min','gra_min');
title('speed estimation on symbol level');


%% perSC review
flag_perSC = 1;

if flag_perSC
    lag_movavg = 8;
    a = 1;
    b = ones(1,lag_movavg)/lag_movavg;

    flag_MaFlt_singleChannel = 0;

    for k = 1:length(sc_active)
        sc_index = sc_active(k);

        figure;
        % accumulated phase rotation across OFDM symbols
        subplot(221);
        plot(deg_total(:,k)/pi);
        xlabel('symbol index');
        ylabel('phase rotation/pi');
        title({['dop simu,phase info, frame#',num2str(i)];...
            ['v=',num2str(v),', N\_sc = ',num2str(N_sc),', Tsc #',num2str(sc_index)]});
        % unwraped phase
        subplot(222);
        plot(deg_unwraped_total(:,k)/pi);
        xlabel('symbol index');
        ylabel('phase rotation/pi');
        title('unwraped phase');

        % angular speed per symbol
        subplot(223);
        hold on;
        plot(w_symbol_diff(:,k));
        plot(w_symbol_gra(:,k));
        legend('diff','gradient');
        xlabel('symbol index');
        ylabel('\omega /pi');
        title('angular speed');

        if flag_MaFlt_singleChannel
            v_symbol_diff_step2 = filter(b,a,v_symbol_diff(:,k));
            v_symbol_gra_step2 = filter(b,a,v_symbol_gra(:,k));
        else
            v_symbol_diff_step2 = v_symbol_diff(:,k);
            v_symbol_gra_step2 = v_symbol_gra(:,k);
        end

        
        % plot ground truth speed and estimated speed
        subplot(224);
        hold on;
        plot(v_symbol_diff_step2);
        plot(v_symbol_gra_step2);

        ylim([-v*1.5 v*1.5]);
        legend('diff','gradient');
        xlabel('symbol index');
        ylabel('relative speed');
        title(['V Est, v=',num2str(v),'N\_sc = ',num2str(N_sc),', Tsc #',num2str(sc_index)]);
    end
end