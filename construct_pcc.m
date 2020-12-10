% load data/ch_ber.mat;
addpath('codes/');
addpath('codes/polar/');
addpath('codes/polar/GA/');

%% Setup Monte-Carlo simulation parameters.
n = 9;
N = 2^n;    % N=512 polar code.
M = 256;    % info bits.
K = 8;     % parity bits,

% Using heuristic construction.
%% Step1: Display channel realibility by graph.
design_Ebn0 = 1.5;
R=(M+K)/N;
sigma_cc = 1/sqrt(2 * R) * 10^(-design_Ebn0/20);
[channels, ~] = GA(sigma_cc, N);  

% Analysis based on BREV order (if necessary). Reconstruct the code.
% channels = bitrevorder(channels);
[~, channel_ordered] = sort(channels, 'descend');
info_bits = sort(channel_ordered(1 : K+M), 'ascend');  
frozen_bits = ones(N , 1);                             
frozen_bits(info_bits) = 0;
info_bits_logical = logical(mod(frozen_bits + 1, 2));

figure;
hold on;
unfrozen_length = M+K;
GA_unfrozen_LLRs = channels(info_bits_logical);
% MonteCarlo_unfrozen_ErrProbs = cnt_err_unfrozen/cnt_total_decoding_attempts;
stem(1:unfrozen_length, 1./GA_unfrozen_LLRs, '*');
%stem(1:unfrozen_length, MonteCarlo_unfrozen_ErrProbs, 'o');
set(gca, 'yscale', 'log');
ylabel('some measure of unreliability');

%% Step2: Extract K burst-err-segments.
err_prob_log10 = zeros(1, M+K);
cutoff_thres = 40;
for k = 1:M+K
    err_prob_log10(k) = convert_to_log10(GA_unfrozen_LLRs(k));
    if err_prob_log10(k) < -cutoff_thres
        err_prob_log10(k) = -cutoff_thres;
    end
end

reliability_diff = diff(err_prob_log10);
[~, i_burst_err_seg_boundaries] = sort(reliability_diff, 'descend');
i_burst_err_seg_boundaries = sort(i_burst_err_seg_boundaries(1:K-1));   % this vector sequentially stores the end index of each burst-err-seg, excepting the last one.
logical_burst_err_seg_boundary = false(1,unfrozen_length);
logical_burst_err_seg_boundary(i_burst_err_seg_boundaries) = true;

figure;
hold on;
scatter(1:unfrozen_length, err_prob_log10, 'b*');
scatter(i_burst_err_seg_boundaries, err_prob_log10(logical_burst_err_seg_boundary), 'ro');

%% Step3: Insert parity bits to construct PCC code.
% Use these boundaries.
burst_err_seg_index = cell(K,1);        % There are K burst-error-segs. seg.index w.r.t non-frozen bits.
parity_bits_index = cell(K,1);          % [1,2,3] means the beginning 3 non-frozen bits form a parity check eqn.
parity_bits_used = false(unfrozen_length,1);

s = 1;
for k = 1:K-1
    burst_err_seg_index{k} = s:1:i_burst_err_seg_boundaries(k);
    parity_bits_index{k} = zeros(k,1);
    s = i_burst_err_seg_boundaries(k)+1;
end
burst_err_seg_index{K} = s:1:unfrozen_length;
parity_bits_index{K} = zeros(K,1);

% Construct parity-check eqns. Protect the most unreliable bits.
for k = 1:K
    sel_cnt = k;
    for sel_iter = 1:sel_cnt
        burst_err_seg = burst_err_seg_index{sel_iter};        % Fetch seg. index w.r.t non-frozen bits.
        GA_LLRs_within_seg = GA_unfrozen_LLRs(burst_err_seg);   % Fetch GA-LLRs within a seg.
        used = parity_bits_used(burst_err_seg);
        if all(used)
            error('All bits in burst_err_seg are used in parity check!');
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
PCC.unfrozen_bits_cnt = unfrozen_length;
PCC.info_bits_cnt = M;
PCC.parity_bits_cnt = K;
PCC.N = N;
PCC.parity_bits_index = parity_bits_index;
PCC.nonfrozen_bits_logical = info_bits_logical;

%% SAVE config e.mat files.
save('data/PCC_config.mat', 'PCC');


%% useful functions.
function ret=convert_to_log10(x)
    if x>5
        ret = -x*log10(exp(1));
    else
        ret = -log10(1+exp(x));
    end
end

function ret_index = find_most_unreliable(is_used, LLRs)
    LLRs_unused = LLRs(~is_used);
    [~,ind] = min(LLRs_unused);
    tmp = find(~is_used, ind);
    ret_index = tmp(end);
end