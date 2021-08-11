clear all;clc;

%% Calculate the blackwell-distribution for AWGN channels.
L = 9;  % 2^L bins, equally-spaced.
N_bins = 2^L;
bin_centers = linspace(0,1,N_bins);
bin_range_half = 0.5/(N_bins-1);

Ebn0 = 1.5;
n = 8;
N = 2^n;
K = 2^(n-1);
K_CRC = 8;
R = K/N;
sigma = 1/sqrt(2*R)*10^(-Ebn0/20);

%% using BPSK modulation, AWGN noise with variance sigma^2.
bm_dist = zeros(1, N_bins);
for iter = 1:N_bins
    bin_left = bin_centers(iter) - bin_range_half;
    if bin_left < 0
        bin_left = 0;
    end
    bin_right = bin_centers(iter) + bin_range_half;
    if bin_right > 1
        bin_right = 1;
    end
    
    % Calculate Pr(bin_left < S <= bin_right).
    if bin_left>0
        y_bin_left = sigma^2/2*log(bin_left/(1-bin_left));
    else
        y_bin_left = -Inf;
    end
    
    if bin_right<1
        y_bin_right = sigma^2/2*log(bin_right/(1-bin_right));
    else
        y_bin_right = Inf;
    end
    
    p = 0.5*(normcdf((y_bin_right-1)/sigma)-normcdf((y_bin_left-1)/sigma)) + ...
        0.5*(normcdf((y_bin_right+1)/sigma)-normcdf((y_bin_left+1)/sigma));
    
    bm_dist(iter) = p;
end

%% Transforms
figure(1);
plot(bin_centers, bm_dist);
hold on;

Wn = up_transform(bm_dist, bm_dist, bin_centers);
Wp = down_transform(bm_dist, bm_dist, bin_centers);
plot(bin_centers, Wp);
plot(bin_centers, Wn);
set(gca, 'yscale', 'log');
legend('W', 'W+', 'W-');

IW = get_I(bm_dist, bin_centers);
IWp = get_I(Wp, bin_centers);
IWn = get_I(Wn, bin_centers);

PeW = get_pe_ML(bm_dist, bin_centers);
PeWp = get_pe_ML(Wp, bin_centers);
PeWn = get_pe_ML(Wn, bin_centers);

fprintf('Channel: AWGN Eb/n0=%f\n IW=%f bits\tIW+=%f bits\tIW-=%f bits\n', Ebn0, IW, IWp, IWn);
fprintf('Pe=%f\tPe+=%f\tPe-=%f\n', PeW, PeWp, PeWn);

%%
I_synth = zeros(1, N);      % fill in this vector last.
Pe_synth = zeros(1, N);
W_mat = zeros(N, N_bins);   % in bit-rev order.
W_mat(1,:) = bm_dist;

for k = 1:n
    % just polarize n times.
    N_target_rows = 2^k;
    for iter = 1:N_target_rows/2
        W = W_mat(iter,:);
        
        Wn = up_transform(W,W,bin_centers);
        Wp = down_transform(W,W,bin_centers);
        
        W_mat(iter, :) = Wn;
        W_mat(iter+N_target_rows/2, :) = Wp;
    end
    fprintf('polarize %d complete\n', k);
end

idx = bitrevorder(1:N);
W_mat = W_mat(idx, :);

%% Calculate Mutual Information for each of the polarized sub-channels.
for k = 1:N
    I_synth(k) = get_I(W_mat(k,:), bin_centers);
    Pe_synth(k) = get_pe_ML(W_mat(k,:), bin_centers);
end

%% Get the ordered sub-channels.
[~, order_bm] = sort(I_synth.', 'descend');
info_bits = sort(order_bm(1 : K+K_CRC), 'ascend');  
frozen_bits = ones(N , 1);                             
frozen_bits(info_bits) = 0;
info_bits_logical_bm = logical(mod(frozen_bits + 1, 2));

% compare with GA.
addpath('codes/polar/GA/');
[channels, ~] = GA(sigma, N);  

% Analysis based on BREV order (if necessary). Reconstruct the code.
% channels = bitrevorder(channels);
[~, order_GA] = sort(channels, 'descend');
info_bits = sort(order_GA(1 : K+K_CRC), 'ascend');  
frozen_bits = ones(N , 1);                             
frozen_bits(info_bits) = 0;
info_bits_logical_GA = logical(mod(frozen_bits + 1, 2));


[~, order_ML] = sort(Pe_synth, 'ascend');
info_bits = sort(order_ML(1 : K+K_CRC), 'ascend');  
frozen_bits = ones(N , 1);                             
frozen_bits(info_bits) = 0;
info_bits_logical_ML = logical(mod(frozen_bits + 1, 2));

diff_1 = sum(xor(info_bits_logical_bm, info_bits_logical_ML));
diff_2 = sum(xor(info_bits_logical_bm, info_bits_logical_GA));

info_bits_logical_bm = info_bits_logical_bm.';
info_bits_logical_ML = info_bits_logical_ML.';
info_bits_logical_GA = info_bits_logical_GA.';

save('blackwell_measure/GA_bm.mat', 'info_bits_logical_bm', 'info_bits_logical_ML', 'info_bits_logical_GA');

%% Realize two measure-transform functions.
function ret_dist = up_transform(dist_1, dist_2, bin_centers)
    N_bins = length(bin_centers);
    ret_dist = zeros(1, N_bins);
    for s1 = 1:N_bins
        for s2 = 1:N_bins
            p = dist_1(s1) * dist_2(s2);    % Pr{S1 \in s1 && S2 \in s2}
            v = bin_centers(s1)*bin_centers(s2)+(1-bin_centers(s1))*(1-bin_centers(s2));
            [~,ind] = min(abs(bin_centers-v));
            
            ret_dist(ind) = ret_dist(ind) + p;
        end
    end
end

function ret_dist = down_transform(dist_1, dist_2, bin_centers)
    % Step1: Decompose symmetric W1 and W2 into sum of BSC channels.
    % Assume that the number N_bins is even.
    
    N_BSC = length(bin_centers)/2;
    BSC_1 = zeros(2, N_BSC);
    BSC_2 = zeros(2, N_BSC);
    
    N_bins = length(bin_centers);
    ret_dist = zeros(1,N_bins);
    
    for k = 1:N_BSC
        BSC_1(1,k) = dist_1(k)*2;
        BSC_2(1,k) = dist_2(k)*2;
        BSC_1(2,k) = bin_centers(k);    % give cross-over probabilities to the second row.
        BSC_2(2,k) = bin_centers(k);
    end
    
    for k1 = 1:N_BSC
        for k2 = 1:N_BSC
            p = BSC_1(2,k1);            % Cross-over probabilities.
            q = BSC_2(2,k2);
            alpha = p*q/(p*q+(1-p)*(1-q));
            beta = (1-p)*q/((1-p)*q+p*(1-q));
            beta = min([beta, 1-beta]); % Ensure beta <= 0.5.
            
            lambda  = BSC_1(1,k1);   
            mu      = BSC_2(1,k2);
            
            c1 = p*q+(1-p)*(1-q);
            c2 = 1-c1;
            [~, ind_c1] = min(abs(bin_centers - alpha));
            [~, ind_c2] = min(abs(bin_centers - beta));
            
            ret_dist(ind_c1) = ret_dist(ind_c1) + 0.5*lambda*mu*c1;
            ret_dist(N_bins+1-ind_c1) = ret_dist(N_bins+1-ind_c1) + 0.5*lambda*mu*c1;
            
            ret_dist(ind_c2) = ret_dist(ind_c2) + 0.5*lambda*mu*c2;
            ret_dist(N_bins+1-ind_c2) = ret_dist(N_bins+1-ind_c2) + 0.5*lambda*mu*c2;
        end
    end
end

function I = get_I(bm_dist, bin_centers)
    N_bins = length(bm_dist);
    % capacity function = 1+xlogx + (1-x)log(1-x).
    I = 0;
    for k = 1:N_bins/2
        p = bin_centers(k); % cross-over probability.
        if p == 0
            Ip = 1;
        else
            Ip = 1+p*log2(p) + (1-p)*log2(1-p);
        end
        I = I + Ip * (2*bm_dist(k));
    end
end

function pe = get_pe_ML(bm_dist, bin_centers)
    N_bins = length(bm_dist);
    pe = 0;
    for k = 1:N_bins/2
        p = bin_centers(k);
        pe_ML = 0.5*(1-abs(1-2*p));
        pe = pe + pe_ML * (2*bm_dist(k));
    end
end

