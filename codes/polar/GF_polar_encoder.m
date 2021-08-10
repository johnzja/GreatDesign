function x = GF_polar_encoder(u, n, GF_info)
% Assume: N = 2^n. u is a column symbol vector or row vector of size N.
% u(i) in GF(2^m).

    alpha = GF_info.alpha;
    add_table = GF_info.add_table;
    mult_table = GF_info.mult_table;
    
    N = 2^n;
    assert(length(u)==N, 'Incorrect length u');
    
    for k = 1:n
        % for each polarization layer
        stride = 2^k;
        stride_2 = 2^(k-1); % stride/2, half of size of kernel input.
        for kern_iter = 1:(N/stride)
            index_1 = (kern_iter-1)*stride+1;
            index_2 = index_1+stride_2;
            for inner_iter = 1:stride_2
                u(index_1) = add_table(u(index_1)+1, mult_table(u(index_2)+1, alpha+1)+1);
                index_1 = index_1 + 1;
                index_2 = index_2 + 1;
            end
        end
    end
    x = u;
end

