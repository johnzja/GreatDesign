function x = my_polar_encoder(u, n)
% Assume: N = 2^n. u is a column vector of size N.
    N = 2^n;
    for k = 1:n
        stride = 2^k;
        for kern_iter = 1:(N/stride)
            index_1 = (kern_iter-1)*stride+1;
            index_2 = index_1+stride/2;
            for inner_iter = 1:(stride/2)
                u(index_1) = u(index_1) + u(index_2);
                index_1 = index_1 + 1;
                index_2 = index_2 + 1;
            end
        end
    end
    x = mod(u,2);
end