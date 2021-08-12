function [idx0, idx1, idx2, idx3] = convert_dist_into_index(s, Nbins)
    ret_idx = zeros(size(s));
    for k=1:4
        ret_idx(k) = floor((Nbins-1)*s(k));
    end
    
    sum_a = Nbins - 1 - sum(ret_idx);
    
    [~, a_order] = sort((Nbins-1)*s - ret_idx);
    if sum_a == 1
        ret_idx(a_order(4)) = ret_idx(a_order(4)) + 1;
    elseif sum_a == 2
        ret_idx(a_order(3)) = ret_idx(a_order(3)) + 1;
        ret_idx(a_order(4)) = ret_idx(a_order(4)) + 1;
    elseif sum_a == 3
        ret_idx(a_order(2)) = ret_idx(a_order(2)) + 1;
        ret_idx(a_order(3)) = ret_idx(a_order(3)) + 1;
        ret_idx(a_order(4)) = ret_idx(a_order(4)) + 1;
    elseif sum_a >= 4
        error('Unexpected error!');
    end
    
    idx0 = ret_idx(1) + 1;
    idx1 = ret_idx(2) + 1;
    idx2 = ret_idx(3) + 1;
    idx3 = ret_idx(4) + 1;
end

