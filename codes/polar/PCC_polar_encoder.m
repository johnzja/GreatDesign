function x = PCC_polar_encoder(info_bits, PCC_conf)
% PCC_POLAR_ENCODER
%   Assume: input info_bits is a row-vector.

    assert(length(info_bits) == PCC_conf.info_bits_cnt, 'input length error!');
    %% Step1: Extract all the locations of parity-check bits.
    K = PCC_conf.parity_bits_cnt;
    M = PCC_conf.info_bits_cnt;
    N = PCC_conf.N;
    nonfrozen_bits_logical = PCC_conf.nonfrozen_bits_logical;
    
    parity_bits_index = PCC_conf.parity_bits_index;
    parity_bit_locs = zeros(K,1);
    info_bit_locs_logical = true(1,M+K);
    
    for k = 1:K
        p = parity_bits_index{k}(end);
        info_bit_locs_logical(p) = false;
        parity_bit_locs(k) = p;
    end
    
    non_frozen_bits = false(1,M+K);
    non_frozen_bits(info_bit_locs_logical) = info_bits;
    
    %% Step2: Insert parity bits.
    for k = 1:K
        p = parity_bit_locs(k);
        check_srcs = parity_bits_index{k}(1:end-1);
        c = mod(sum(non_frozen_bits(check_srcs)),2);
        non_frozen_bits(p) = c;
    end
    
    u = zeros(1, N);    % Strange fact: HERE logical type doesn't work.
    u(nonfrozen_bits_logical) = non_frozen_bits;
    x = my_polar_encoder(u,log2(N));
end

