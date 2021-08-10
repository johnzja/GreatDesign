load('C:/Users/John Zhu/Desktop/instant_files/bm/bm_p128_N256.mat');
fprintf('Load BM data complete.\n');

%% Get capacities.
n = 8;
N = 2^n;

if ~exist('bin_centers', 'var')
    N_bins_each_dim = 128;
    bin_centers = linspace(0, 1, N_bins_each_dim);
end

polarized_channels = Synth_channels{n+1};
polarized_channels = {polarized_channels{bitrevorder(1:N)}};

C = zeros(N,1);
for idx = 1:N
    C(idx) = get_I_4D(polarized_channels{idx}, bin_centers);
end

I_total = sum(C);

%% 
total_capacity_error = (I_total - N*I)/(N*I);
fprintf('Total capacity error = %.2f %%\n', total_capacity_error*100);

%% Code construction, R=1/2.

[~, order] = sort(C, 'descend');
M = 128;
R = M/N;

info_syms = sort(order(1:M), 'ascend').';

% construct code by binary Gaussian-Approximation (GA).
% compare with GA.
addpath('codes/polar/GA/');
[channels, ~] = GA(sigma, N);  

[~, order_GA] = sort(channels, 'descend');
info_bits = sort(order_GA(1 : M), 'ascend'); 
info_bits_logical = false(1,N);
info_bits_logical(info_bits) = true;
frozen_bits = ~info_bits_logical;

%% Required simulation.
GF_info.m = 2;          % GF(2^2).
GF_info.alpha = 2;
GF_info.mult_table  = [0, 0, 0, 0; 0, 1, 2, 3; 0, 2, 3, 1; 0, 3, 1, 2];
GF_info.add_table   = [0, 1, 2, 3; 1, 0, 3, 2; 2, 3, 0, 1; 3, 2, 1, 0];
GF_info.kernel_index_vec0 = kernel_index_vec0;          % boolean vector.
GF_info.kernel_index_vec1 = int32(kernel_index_vec1);   % int32 vector.
GF_info.kernel_index_mat0 = int32(kernel_index_mat0);   % int32 matrix.
GF_info.kernel_index_mat1 = int32(kernel_index_mat1);   % int32 matrix.


Ebn0 = 1.5;
L = 16;

addpath('sim/');
% if code rate=1/2, in QPSK channel Es/n0 = Eb/n0.

% simulation.
Ebn0_arr = [1.25, 1.5, 1.75, 2.00, 2.25, 2.5, 2.75, 3.00];
min_errors = 1200;
N_Ebn0 = length(Ebn0_arr);
BLERs = zeros(1, N_Ebn0);

for idx = 1:N_Ebn0
    BLERs(idx) = sim_Qary_SCL(N, M, false, info_bits_logical, Ebn0, min_errors, L, GF_info);
end


save('data/bm_construction_q4_L16.mat', 'Ebn0_arr', 'BLERs', 'L', 'info_bits_logical', 'N', 'M', 'm');
