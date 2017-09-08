% this function find the index of frames through offline synchronization
%   input:
%          -sig: signal
%          -preamble: preamble
%          -len_frame: length of a frame
%   output:
%          -index_arr: an array of index of starting point of frames
function index_arr = sync_offline(sig,preamble,len_frame)
    
    threshold = 1.8;
    
    flag_debug = 1;
    index_arr = [];
    % matched filter, which works like correlation function
    coef_MF_preamble = preamble(end:-1:1);
    data_MFflted = filter(coef_MF_preamble,1,sig);
    
    if flag_debug
        figure;
        plot(data_MFflted);
    end
    % get upper and lower envolop of matched filter output
    [up,lo] = envelope(data_MFflted,1000,'peak');
    if flag_debug
        figure;
        hold on;
        plot(data_MFflted);
        plot(up);
        title('index detection')
    end
    len_sig = length(sig);
    len_win = 10*length(preamble);
    i = 1;
    while i < len_sig
        if i ~=1
            % for each search segment, there is a overlap of 0.1 search windown length
            i = i - round(0.1*len_win);
        end
        % set the target search window, search a small windown every time
        % if there is a peak which is obviously higher than envolop in the
        % window, then set it as a candidate
        i_tail = i + len_win;
        if i_tail > len_sig
            i_tail = len_sig;
        end
        idx = find(data_MFflted(i:i_tail) > threshold*mean(up(i:i_tail)));

        if ~isempty(idx)
            idx = idx + i - 1;
            [y_temp,i_max] = max(data_MFflted(idx));
            index_arr = [index_arr;idx(i_max)];
        end
%         %if k-means need to be introduced, since only 1.2 frame length, so maximal 2 clusters 
%         if std(idx) > 1000
%             cluster = kmeans(idx,2);
%             for k =1:2
%                 idx_candidate = idx(cluster == k);
%                 [y_temp,i_max] = max(data_MFflted(idx_candidate));
%                 index_arr = [index_arr;idx_candidate(i_max)];
%             end
%         else %if not
%             [y_temp,i_max] = max(data_MFflted(idx));
%             index_arr = [index_arr;idx(i_max)];
%         end
        
        i = i_tail +1 ;
    end
    
    % sort index_arr and remove duplication
    index_arr = sort(index_arr);
    idx_dup = find(diff(index_arr)<len_frame*0.5);
    index_arr(idx_dup) = [];
end