% Construct different PCC codes for different K.
% B segments
% K additional parity bits.
% GA construction.
% M info bits.
% 0.5*K(K-1) <= M


% load data/ch_ber.mat;
addpath('codes/');
addpath('codes/polar/');
addpath('codes/polar/GA/');
M = 256;
N = 512;
n = log2(N);

K_max = floor(0.5*(1+sqrt(1+8*M)));
% Construct K_max different PCCs.


for K = 1:K_max
    % burst_err_seg is a vector of length K.
    %% Step1: Display channel realibility by graph.
    design_Ebn0 = 1.5;
    R=(M+K)/N;
    sigma_cc = 1/sqrt(2 * R) * 10^(-design_Ebn0/20);
    [channels, ~] = GA(sigma_cc, N);  

    % Analysis based on BREV order (if necessary). Reconstruct the code.
    % channels = bitrevorder(channels);

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
    PCC(K).unfrozen_bits_cnt = M+K;
    PCC(K).info_bits_cnt = M;
    PCC(K).parity_bits_cnt = K;
    PCC(K).N = N;
    PCC(K).parity_bits_index = parity_bits_index;      % contain parity equations.
    PCC(K).nonfrozen_bits_logical = info_bits_logical; 
end

disp('PCC code construction complete.');
Ebn0 = 1.25;

%% Simulation Loop.
max_errors = 1600;
BLERs = zeros(K_max, 1);
parfor K = 1:K_max
    BLERs(K) = sim_PCC(PCC(K), Ebn0, max_errors, 16);
    fprintf('Simulation complete for K = %d\n', K);
end

disp('Simulation complete');

save(['data/PCC_K_Ebn0=', num2str(Ebn0), '.mat'], 'BLERs', 'PCC');

%% Useful Functions.
function ret_index = find_most_unreliable(is_used, LLRs)
    LLRs_unused = LLRs(~is_used);
    [~,ind] = min(LLRs_unused);
    tmp = find(~is_used, ind);
    ret_index = tmp(end);
end
