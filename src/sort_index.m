function index = sort_index(seq,index_arr,k)
% this function first cluster the index_arr into k clusters, and then select 
% largest peak among each cluster as index
% seq: corr sequency
% index_arr: arr of index
% frame_len: frame length
% k: number of clusters
    if isrow(index_arr)
        idx = kmeans(index_arr',k);
    else
        idx = kmeans(index_arr,k);
    end
    index = zeros(k,1);
    idx_unique = unique(idx,'stable');
    for i = 1:k
        temp = idx_unique(i);
        cluster = index_arr(find(idx == temp));
        [y temp_i]= max(seq(cluster));
        index(i) = cluster(temp_i);
        if i == 1
            continue;
        end
%         disp(index(i)-index(i-1));
    end
end