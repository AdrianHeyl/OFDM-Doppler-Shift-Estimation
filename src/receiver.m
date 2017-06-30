clear all, close all, clc;

%% setting
Fs = 44100;   % sapmling frequency
N_sc = 1024;  % number of carriers
bw_sc = Fs/N_sc;    %bandwidth of each subcarrier
ifft_size = 2*N_sc;
bidata_len = N_sc;  %binary data length
cp_length = 256;  % the pre and post fix are the same length of ifft_size, hereby I used both cyclic prefix and postfix
symbolCP_len = ifft_size + 2*cp_length;    % symbol length with prefix and postfix
blank_len = 100;    % there is a blank interval between two frames
N_symbol = 100;     % number of symbol in a frame
N_frame = 2;        % number of frames in generated audio file
t = [1:symbolCP_len]'/Fs;
R_temp = 0.2;

sc_mask = zeros(N_sc,1);   % subcarrier mask
sc_active = [100 200 300 400 500];    % active subcarrier index
% sc_active = [100 200 300 400 500 600 700 800 900 1000];    % active subcarrier index
%sc_active = [1100 1200 1300 1400 1500 1600 1700 1800 1900 2000];    % active subcarrier index
sc_mask(sc_active) = 1;    % only active subcarriers are enabled to transmit data, others are blocked
N_sca = length(sc_active);  % number of active subcarriers
deg = zeros(N_symbol,N_sca);    % angle of each symbol, unit degree

sync_offset = 0;
flag_figure = 0;    % to display a lot of figures or not

N_j = 100;      % if we want to see multiple symbols, choose the number of symbols we want to display

%% preamble
f_min = 8000;   % min frequency
f_max = 12000;  % max frequency
pre_len = 256;  % length of preamble
t_prehalf = [1:pre_len/2]/Fs;
t_lasthalf = [pre_len/2+1:pre_len]/Fs;
preamble = [chirp(t_prehalf,f_min,pre_len/2/Fs,f_max), chirp(t_lasthalf,f_max,pre_len/Fs,f_min)]';


%% loading audio file from hard drive
%13-03-46
fileName = 'C:\Users\Erdo\Desktop\Designing studies\export\17-18-53.wav';
[sig_received,fs_temp] = audioread(fileName);
load('bias.mat');
%% ---------------demodulation-------------------------------
coef_MF_preamble = preamble(end:-1:1);  % coeffient of matched filter
% sync_threshold = [105 105 83 ];
  
%% synchronization
% filt received signal with matched filter, the peak in
% data_MFflted indicates the start of a frame
data_MFflted = filter(coef_MF_preamble,1,sig_received);

% this figure is used to oberserve the peak value manully
figure;
plot(data_MFflted);

% this threshold should be adjusted according to different snr level
sync_threshold = 0.4;
disp(size(find(data_MFflted > sync_threshold)));
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
    
    for k = 1:N_sca
        q = sc_active(k)+1;
        deg(j,k) = angle(BPSK_demodulated(q));    % calculate the angle of a symbol
    end
%         v_angular(j) = deg - deg_pre;
%         deg_pre = deg;

    % plot constellation for each OFDM symbol
    if flag_figure
        figure;
        xlim([-R_temp R_temp]);
        ylim([-R_temp R_temp]);
        hold on;
%         plot(BPSK_data,'ro','MarkerSize',10,'linewidth',1.5);
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
        title ({['Exp, N\_sc = ',num2str(N_sc),', symbol:',num2str(j)]});
        print(['symbol#',num2str(j)],'-dpng');
    end

end  % end of loop j

%loop through all the active subcarriers
for k =1:N_sca
    % show accumulate phase rotation across all the symbols in one frame
    figure;
%     subplot(1,2,1);
    %plot(unwrap(deg(:,k))/pi-bias-bias(1));
    plot(unwrap(deg(:,k))/pi);
    xlabel('symbol index');
    ylabel('accumulated phase rotation/pi');
    title(['Exp, N\_sc = ',num2str(N_sc),', Asc=',num2str(sc_active(k))]);
%     subplot(1,2,2);
%     plot(unwrap(deg)/pi-bias);
%     title('subtracting bias');
    v_angular = gradient(unwrap(deg)/pi)';
    disp(['angular v = ',num2str(mean(v_angular)),'pi, std=',num2str(std(v_angular))]);

end
