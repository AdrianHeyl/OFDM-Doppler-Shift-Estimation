%{ 
this file simulated the effect of snr, STO, CFO
- snr is Signal to Noise  Ratio
- STO is Symbol Time Offset, which is time domain sync error,
- CFO is Carrier Frequency Offset, it is a normalized value: 
   (frequency shift)/(bandwidth of each subcarrier)
   and it is identical across all the subcarriers, doppler effect actually
   introduces different CFO on different subcarriers
%}
clear all, close all, clc;

%% -----------generate cplx modulated ofdm signal, proof it ---------------

%% setting
Fs = 44100;   % sapmling frequency
N_sc = 2048;  % number of carriers
bw_sc = Fs/N_sc;    %bandwidth of each subcarrier
ifft_size = 2*N_sc;
bidata_len = N_sc;  %binary data length
cp_length = 256;  % the pre and post fix are the same length of ifft_size, hereby I used both cyclic prefix and postfix
symbolCP_len = ifft_size + 2*cp_length;    % symbol length with prefix and postfix
blank_len = 100;    % there is a blank interval between two frames
N_symbol = 100;     % number of symbol in a frame
N_frame = 2;        % number of frames in generated audio file
t = [1:symbolCP_len]'/Fs;
deg = zeros(N_symbol,1);    % angle of each symbol, unit degree

sc_mask = zeros(N_sc,1);   % subcarrier mask
sc_active = [100 1000];    % active subcarrier index
sc_mask(sc_active) = 1;    % only active subcarriers are enabled to transmit data, others are blocked

snr = 50;   % in unit of dB
sync_offset = 0;    % STO value, in unit of number of samples, e.g. 5 means STO = 5 samples, -5 means STO = -5 samples

flag_figure = 0;    % to display a lot of figures or not
flag_singlesymbol = 0;   % if we want only to have a look for single symbol
flag_CFO = 0;       % if we want to introduce CFO 

if flag_singlesymbol 
    N_j = 1;        % we only want to see single symbol
else
    N_j = 100;      % if we want to see multiple symbols, choose the number of symbols we want to display
end

if flag_CFO    % introduced CFOs, could be an array; 
    CFO = [0.1:0.1:0.9];      % if multiple CFO values, it will loop through all the CFO values
    N_i = length(CFO);        % round of iteration of CFO values
%     symbol_CP = add_CFO_audio(symbol_CP,CFO,ifft_size);
else        % if you do not want to introduce CFO
    N_i = 1;        % round of iteration of CFO values
    CFO = 0;        % CFO is 0, means there is no frequency shift between Tx and Rx
end

%% modulation
% BPSKModulator
hMod = comm.BPSKModulator;    % creating bpsk modulator system object
hMod.PhaseOffset = pi/16;     % initial phase set to pi/16

binary_data = ones(N_sc,1);   % binary data are all 1
BPSK_data = step(hMod,binary_data);   % bpsk modulated data
BPSK_data = BPSK_data.*sc_mask;    % use mask to disable unwanted subcarriers

% build up conjugate symmetry spectrum of Tx for ifft
% negelect this part, what you need to know is that TxSpectrums is
% frequency data which will be used in ifft modulation of OFDM
carriers = [2:N_sc+1]';
NegCarriers = ifft_size - carriers + 2;	%find the bins for the negative frequency carriers
TxSpectrums = zeros(ifft_size,1);
TxSpectrums(carriers) = BPSK_data;
TxSpectrums(NegCarriers) = conj(BPSK_data);
% ifft modulation 
symbol = ifft(TxSpectrums);
phase_ori = angle(symbol)/pi;
phase_ori_unwraped = unwrap(angle(symbol))/pi;
% adding CP
symbol_CP = [symbol(end - cp_length + 1 : end);...
             symbol;...
             symbol(1:cp_length)]; % add cyclic prefix and postfix

% scale the amplitude of symbol to 0.8
scale = 0.8;
symbol_CP = symbol_CP/max(abs(symbol_CP))*scale;

%show OFDM symbol
figure;
hold on;
plot(symbol_CP);
plot([cp_length+1:cp_length + ifft_size],symbol/max(abs(symbol))*scale);
legend('symbol with CP','symbol');
xlabel('sample index');
ylabel('amplitude');
title('OFDM symbol in time domain');


%% add preamble and make frame
% preamble is a segment of samples known to both Tx and Rx, used for time
% sync, here we use chirp as preamble, refer to wikipedia for further info
f_min = 8000;   % min frequency
f_max = 12000;  % max frequency
pre_len = 256;  % length of preamble
t_prehalf = [1:pre_len/2]/Fs;
t_lasthalf = [pre_len/2+1:pre_len]/Fs;
preamble = [chirp(t_prehalf,f_min,pre_len/2/Fs,f_max), chirp(t_lasthalf,f_max,pre_len/Fs,f_min)]';

frame = [preamble;zeros(blank_len,1);repmat(symbol_CP,N_symbol,1)];
seq = [zeros(Fs*0.5,1);frame; zeros(Fs*1,1); frame];

% show audio signal with 2 frames
figure;
plot([1:length(seq)]/Fs,seq);
xlabel('t/s');
ylabel('amplitude');
title('audio OFDM signal with 2 frames');

%% AWGN channel, introduce noise
%sig_awgn = awgn(seq,snr);
sig_awgn = seq;

%% ---------------demodulation-------------------------------
coef_MF_preamble = preamble(end:-1:1);  % coeffient of matched filter
% sync_threshold = [105 105 83 ];
for i = 1:N_i   % loop of CFO values
    %% introduce CFO in received signal
    if flag_CFO  % if CFO exist
        sig_received = add_CFO_audio(sig_awgn,CFO(i),ifft_size);
    else    % if CFO does not exist
        sig_received = sig_awgn;
    end
    
    %% synchronization
    % filt received signal with matched filter, the peak in
    % data_MFflted indicates the start of a frame
    data_MFflted = filter(coef_MF_preamble,1,sig_received);

    % this figure is used to oberserve the peak value manully
%     figure;
%     plot(data_MFflted);

    % this threshold should be adjusted according to different snr level
    sync_threshold = 50;
    index_temp = find(data_MFflted > sync_threshold);
    % a simple function used to search for index of starting sample of a frame
    index_arr = sort_index(data_MFflted,index_temp,N_frame);  
    disp('index_arr');
    disp(index_arr);
    
%     deg_pre = 0;
    for j = 1: N_j  % loop of symbols in a frame
        index_symbol = j;
        % index of start and end of a symbol
        i_start = index_arr(1) + blank_len + symbolCP_len*(index_symbol-1) + 1 + sync_offset ;
        i_end = i_start + symbolCP_len - 1;
        target_sym = sig_received(i_start: i_end);  % target symbol

        symbol_woCP = target_sym(cp_length + 1 : end - cp_length); % remove CP
        %compute phase information
        phase_demodulated = angle(symbol_woCP)/pi;
        phase_demodulated_unwraped = unwrap(angle(symbol_woCP))/pi;
        % demodulation, map data to complex plane, constellation
        frequency_data = fft(symbol_woCP);
        BPSK_demodulated = frequency_data(1:N_sc);
        BPSK_demodulated = BPSK_demodulated / (max(abs(BPSK_demodulated)));

        q = sc_active(2)+1;
        deg(j) = angle(BPSK_demodulated(q));    % calculate the angle of a symbol

        % plot constellation for each OFDM symbol
        if flag_figure
            figure;
            xlim([-1 1]);
            ylim([-1 1]);
            hold on;
            plot(BPSK_data,'ro','MarkerSize',10,'linewidth',1.5);
            plot(BPSK_demodulated,'x','MarkerSize',10,'linewidth',1.5);
            for p = 1:N_sc
                if sc_mask(p)
                    x = real(BPSK_demodulated(p));
                    y = imag(BPSK_demodulated(p));
                    text(x,y,num2str(p-1));
                    x = real(BPSK_demodulated(p+1));
                    y = imag(BPSK_demodulated(p+1));
                    text(x,y,num2str(p));
                    x = real(BPSK_demodulated(p+2));
                    y = imag(BPSK_demodulated(p+2));
                    text(x,y,num2str(p+1));
                end
            end
            title ({['AWGN, snr = ',num2str(snr),'dB, N\_sc = ',num2str(N_sc)];...
                ['STO = ',num2str(sync_offset),', CFO=',num2str(CFO(i)),', sc:'...
                ,num2str(sc_active),', symbol:',num2str(j)]});
%             print(['symbol#',num2str(j)],'-dpng');
        end

    end  % end of loop j
    
    % show accumulate phase rotation across all the symbols in one frame
    figure;
    plot(unwrap(deg)/pi);
    xlabel('symbol index');
    ylabel('accumulated phase rotation/pi');
    title(['AWGN, snr = ',num2str(snr),'dB, N\_sc = ',num2str(N_sc),', Asc=',num2str(q-1),' CFO=',num2str(CFO(i))]);
    v_angular = gradient(unwrap(deg)/pi)';
    disp(['CFO = ',num2str(CFO(i)), ' angular v = ',num2str(mean(v_angular)),'pi, std=',num2str(std(v_angular))]);
    
    close all
    figure;
    plot(v_angular)
end  % end of loop i




