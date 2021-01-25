function [PCC_structs, ch_reliability] = get_standard_PCC(N, M, design_Ebn0)
% Construct different PCC codes for different K. K: additional parity bits.
% B segments
% standard GA construction.
% M info bits.
% 0.5*K(K-1) <= M

% "Standard" means to design the "burst-error segs" uniformly and
% parity-check the most unreliable bits.

    addpath('codes/');
    addpath('codes/polar/');
    addpath('codes/polar/GA/');

    K_max = floor(0.5*(1+sqrt(1+8*M)));
    % Construct K_max different PCC concatenation schemes.

    for K = 1:K_max
        R=(M+K)/N;
        sigma_cc = 1/sqrt(2 * R) * 10^(-design_Ebn0/20);
        [channels, ~] = GA(sigma_cc, N);  

        [~, channel_ordered] = sort(channels, 'descend');
        info_bits = sort(channel_ordered(1:K+M), 'ascend');  
        frozen_bits = ones(N , 1);                             
        frozen_bits(info_bits) = 0;
        info_bits_logical = logical(mod(frozen_bits + 1, 2));

        %% Step2: Construct valid burst-error segs.
        % subject to constraint:  S1 + S2 + ... + SK == M+K &&  Si >= K-(i-1).
        % Using uniformly-distributed method.

        must_use = (K+1)*K/2;
        remain = M+K-must_use;
        distribute = floor(remain/K);

        burst_err_seg_index = cell(K,1); 
        parity_bits_index = cell(K,1);
        s = 1;
        for k = 1:K-1
           burst_err_seg_index{k} = (s:s+K-k+distribute);
           s = s+K-k+1+distribute;
        end
        burst_err_seg_index{K} = (s:M+K);

        % Construct parity-check eqns. Protect the most unreliable bits.
        GA_unfrozen_LLRs = channels(info_bits_logical);
        parity_bits_used = false(M+K,1);
        for k = 1:K
            sel_cnt = k;
            for sel_iter = 1:sel_cnt
                burst_err_seg = burst_err_seg_index{sel_iter};        % Fetch seg. index w.r.t non-frozen bits.
                GA_LLRs_within_seg = GA_unfrozen_LLRs(burst_err_seg);   % Fetch GA-LLRs within a seg.
                used = parity_bits_used(burst_err_seg);
                if all(used)
                    error('All bits in a certain burst_err_seg are used in parity check!');
                end

                % Find the least reliable bits in "unused".
                t = find_most_unreliable(used, GA_LLRs_within_seg);
                protected_info_bit_index = burst_err_seg(t);

                tmp = parity_bits_index{k};
                tmp(sel_iter) = protected_info_bit_index;
                parity_bits_index{k} = tmp;

                parity_bits_used(protected_info_bit_index) = true;
            end
        end

        % Setup PCC-structure.
        PCC_structs(K).unfrozen_bits_cnt = M+K;
        PCC_structs(K).info_bits_cnt = M;
        PCC_structs(K).parity_bits_cnt = K;
        PCC_structs(K).N = N;
        PCC_structs(K).parity_bits_index = parity_bits_index;      % contain parity equations.
        PCC_structs(K).nonfrozen_bits_logical = info_bits_logical; 
    end
    ch_reliability = channels;
end