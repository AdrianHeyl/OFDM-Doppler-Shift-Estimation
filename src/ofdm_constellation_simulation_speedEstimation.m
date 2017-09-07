
clear all, close all, clc;

%% setting

% if apply moving average filter to v_symbol on each subcarrier
flag_MAFlt_perSC = 1;

t_duration = 20;  
%golay sequence
load('golaySeq.mat');   

N_sc = 256;  %number of carriers
c = 343;  % speed of sound in the air

filename = '../res/exp07-09-2017/setup2-sc256-No1.wav';

% select step_sc accoring to N_sc, step_sc is the distance between active
% subcarriers 
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
cp_length = 256;  % the pre and post fix are the same length of ifft_size
symbolCP_len = ifft_size + cp_length;    % len with prefix and postfix
blank_len = 100;
N_symbol = 50;     % number of symbol in a frame
N_frame = ceil(t_duration*Fs/(N_symbol*symbolCP_len));
t = [1:symbolCP_len]'/Fs;

% preamble
f_min = 10000;
f_max = 14000;
pre_len = 1024;
t_prehalf = [1:pre_len/2]/Fs;
t_lasthalf = [pre_len/2+1:pre_len]/Fs;
preamble = [chirp(t_prehalf,f_min,pre_len/2/Fs,f_max), chirp(t_lasthalf,f_max,pre_len/Fs,f_min)]';


sc_mask = zeros(N_sc,1);   % subcarrier mask
sc_active = [5:step_sc:N_sc];    % active subcarrier index
sc_mask(sc_active) = 1;

hMod = comm.BPSKModulator;    % creating bpsk modulator system object
hMod.PhaseOffset = pi/16;     % phase set to pi/16

binary_data = ones(N_sc,1);
binary_data(sc_active) = golay; % the binary data transmitted is golay sequence
BPSK_data = step(hMod,binary_data).*sc_mask;   % this is the bpsk modulated data

len_frame = pre_len + blank_len*2 + symbolCP_len*N_symbol;
sync_offset = -cp_length/2;

[sig_received fs] = audioread(filename);


%% synchronization
% new function developed for synchronization
index_arr = sync_offline(sig_received,preamble,len_frame);
disp('index_arr');
disp(index_arr);
% the diff of index_arr should be around the length of frame, other wise
% synchronization is wrong
disp('diff');
disp(diff(index_arr));

% number of target frames detected by synchronization
N_targetFrame = length(index_arr);

%% phase information extraction
% following variables are all matrix, rows --> symbols, columns --> sc_active
v_est_diff = []; % speed estimated based on diff(phase)
v_est_grad = []; % speed estimated based on gradient(phase)
deg_total = [];  % original phase info 
deg_unwraped_total = []; % unwraped phase info (unwraped in frame, not between frames)
deg_unwraped_woSmooth_total = []; % without sooth
w_symbol_diff = []; % \omiga per symbol based on diff
w_symbol_gra = []; %\omiga per symbol based on gradient

% loop through multiple frames
for i_frame = 1:N_targetFrameS
    % create a matrix to record the phase info in a frame
    deg = zeros(N_symbol,length(sc_active));
    
    %loop through OFDM symbols
    for j = 1:N_symbol
        index_symbol = j;
        i_start = index_arr(i_frame) + blank_len + symbolCP_len*(index_symbol-1) + 1;
        i_end = i_start + symbolCP_len - 1;
        target_sym = sig_received(i_start: i_end);
        symbol_woCP = target_sym(cp_length + 1 + sync_offset: end + sync_offset);

        % map data to complex plane, constellation
        frequency_data = fft(symbol_woCP);
        BPSK_demodulated = frequency_data(2:N_sc+1);
        BPSK_demodulated = BPSK_demodulated / max(abs(BPSK_demodulated(sc_active)));

        deg(j,:) = angle(BPSK_demodulated(sc_active));
        % this segment is for debug, to inspect constellation symbol by symbol
        if j>40 && 0 && i_frame ==1
            figure;
            xlim([-1.2 1.2]);
            ylim([-1.2 1.2]);
            hold on;
            plot(BPSK_data,'ro','MarkerSize',10,'linewidth',1.5);
            plot(BPSK_demodulated,'x','MarkerSize',10,'linewidth',1.5);

            x = real(BPSK_demodulated(sc_active));
            y = imag(BPSK_demodulated(sc_active));
            text(x,y,num2str(sc_active'));

            title (['Doppler simu ,speed=',num2str(v),'m/s, symbol #',num2str(j)]);
%             print(['symbol#',num2str(j)],'-dpng');
        end
    end  % end of loop j

    %% angular speed extraction through phase rotation in every frame
    deg_total = [deg_total;deg];
    deg_unwraped_woSmooth = unwrap(deg);
    deg_unwraped_woSmooth_total = [deg_unwraped_woSmooth_total;deg_unwraped_woSmooth];
    % options, if smooth the unwarped phase
    for k = 1:length(sc_active)
        deg_unwraped(:,k) = smooth(unwrap(deg(:,k)),'moving');
    end
    deg_unwraped_total = [deg_unwraped_total;deg_unwraped];
    
    % angular speed based on phase diff(unit: pi)
    diff_temp = diff(deg_unwraped/pi);
    w_symbol_diff = [w_symbol_diff;diff_temp(1,:);diff_temp];
    
    % angular speed based on phase gradient(unit: pi)
    [f_temp,w_symbol_temp] = gradient(deg_unwraped/pi);
    w_symbol_gra = [w_symbol_gra;w_symbol_temp];
   
end

% this is for debug, see the effect of smooth on the phase
for k = 1:length(sc_active)
    figure;
    hold on;
    plot(deg_unwraped_woSmooth_total(:,k));
    plot(deg_unwraped_total(:,k));
    legend('wosmooth','smooth');
end

%% speed estimation based on phase information
factor = 2.85;
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
% according to the variation of gradient of the curves
[fx_diff,fy_diff] = gradient(v_symbol_diff_step2);
[fx_gra,fy_gra] = gradient(v_symbol_gra_step2);


% select N_sc_avg subcarriers with least peaks to calculate the avg 
N_sc_avg = 5;
[temp1, ix_diff] = sort(var(fy_diff));
[temp2, ix_gra] = sort(var(fy_gra));

% averaged symbol speed estimation across N_sc_avg subcarriers
v_symbol_diff_avg = mean(v_symbol_diff_step2(:,ix_diff(1:N_sc_avg)),2);
v_symbol_gra_avg = mean(v_symbol_gra_step2(:,ix_gra(1:N_sc_avg)),2);
% pick out the subcarrier with least fluctuation
v_symbol_diff_min = v_symbol_diff_step2(:,ix_diff(1));
v_symbol_gra_min = v_symbol_gra_step2(:,ix_gra(1));

dt = [1:length(v_symbol_diff_avg)]*symbolCP_len/Fs;

% show symbol speed estimation across all the target frames 
figure;
hold;
% plot(dt,v_symbol_diff_avg);
% plot(dt,v_symbol_gra_avg);
plot(dt,v_symbol_diff_min);
plot(dt,v_symbol_gra_min);
% legend('diff\_avg','gra\_avg','diff\_min','gra\_min');
legend('diff\_min','gra\_min');


% transform symbol speed to sample speed, to compare with groud truth, not
% usefull right now since we do not have a ground truth
flag_symbol2sample = 0 ;
if flag_symbol2sample
    % transform speed est from symbol level to sample level
    v_sample_gra_avg = symbol2sample(v_symbol_gra_avg,pre_len,blank_len,N_symbol,symbolCP_len);
    v_sample_diff_avg = symbol2sample(v_symbol_diff_avg,pre_len,blank_len,N_symbol,symbolCP_len);
    v_sample_gra_min = symbol2sample(v_symbol_gra_min,pre_len,blank_len,N_symbol,symbolCP_len);
    v_sample_diff_min = symbol2sample(v_symbol_diff_min,pre_len,blank_len,N_symbol,symbolCP_len);

    %% ground truth of relative speed
    % line case: target start at (-5,1), towards positive direction of x axis,
    % with different speed
    vr = [];
    t_vr = [1:len_frame*N_targetFrame]/Fs;
    switch typ
        case 'l'
            theta = atan(1./(v*t_vr-4.825));
            theta(theta<0) = theta(theta<0) + pi;
            vr = -1*v*cos(theta);
        case 'r'
            disp('TBD!');
        otherwise
            disp('typ error!');
    end

    win_len = 1024;
    win_step = 1024;
    rmse_fft = 1000;
    rmse_fft_final = 1001;
    for win_len = 1024:1024:Fs*2

        if flag_noise
            v_sinefft = sineFFTSpeed(v,win_len,snr);
        else
            v_sinefft = sineFFTSpeed(v,win_len);
        end

        if length(v_sinefft) < length(t_vr)
            v_sinefft = [v_sinefft;ones(length(t_vr)-length(v_sinefft),1)*v_sinefft(end)];
        else
            v_sinefft = v_sinefft(1:length(t_vr));
        end
        rmse_fft = sqrt(immse(vr',v_sinefft));
    %     disp([num2str(win_len),',',num2str(rmse_fft)]);
        if rmse_fft < rmse_fft_final
            rmse_fft_final = rmse_fft;
            v_sinefft_final = v_sinefft;
            win_len_final = win_len;
        end
    end


    figure;
    hold on;
    % plot(t_vr,v_sample_gra_avg,'r');
    % plot(t_vr,v_sample_gra_min,'g');
    plot(t_vr,v_sinefft_final,'m');
    plot(t_vr,vr,'b');
    plot(t_vr,ones(len_frame*N_targetFrame,1)*v,'k--');
    plot(t_vr,ones(len_frame*N_targetFrame,1)*v*(-1),'k--');
    ylim([-v*1.5 v*1.5]);
    legend('gra-avg','gra-min','fft','relative speed');
    title(['speed estimation, v=',num2str(v),' N\_sc = ',num2str(N_sc)]);

    % consider delay when applying MA filter
    if flag_MAFlt_perSC
        len_delay = symbolCP_len*lag_movavg/2;
        plot(t_vr(1:end-len_delay),v_sample_gra_avg(len_delay+1:end),'r');
        plot(t_vr(1:end-len_delay),v_sample_gra_min(len_delay+1:end),'g');
        rmse_gra_avg = sqrt(immse(vr(1:end-len_delay)',v_sample_gra_avg(len_delay+1:end)));
        rmse_gra_min = sqrt(immse(vr(1:end-len_delay)',v_sample_gra_min(len_delay+1:end)));
    else
        plot(t_vr,v_sample_gra_avg,'r');
        plot(t_vr,v_sample_gra_min,'g');
        rmse_gra_avg = sqrt(immse(vr',v_sample_gra_avg));
        rmse_gra_min = sqrt(immse(vr',v_sample_gra_min));
    end

    disp('RMSE of gra_avg, gra_min,fft:');
    disp(num2str(rmse_gra_avg));
    disp(num2str(rmse_gra_min));
    disp(num2str(rmse_fft_final));
end

%% perSC review
flag_perSC = 1;

if flag_perSC

    for k = 1:length(sc_active)
        sc_index = sc_active(k);

        figure;
        % accumulated phase rotation across OFDM symbols
        subplot(221);
        plot(deg_total(:,k)/pi);
        xlabel('symbol index');
        ylabel('phase rotation/pi');
        title(['phase info, N\_sc = ',num2str(N_sc),', Tsc #',num2str(sc_index)]);
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

        % plot estimated speed per symbol
        subplot(224);
        hold on;
        plot(v_symbol_gra(:,k));
        plot(v_symbol_gra_step2);
        legend('gra-ori','gra-MA');
        title(['sc #',num2str(sc_active(k))]);
    end
end
