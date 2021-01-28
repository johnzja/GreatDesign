function x = my_polar_encoder(u, n)
% Assume: N = 2^n. u is a column vector or row vector of size N.
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
                u(index_1) = u(index_1) + u(index_2);
                index_1 = index_1 + 1;
                index_2 = index_2 + 1;
            end
        end
    end
    x = mod(u,2);   % perform mod-2 sum only once, reducing MATLAB complexity.
end