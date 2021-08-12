fprintf('Start loading BM data...\n');
addpath('C:/Users/dell/Desktop/instant_files/bm/exp3');
load('bm_synth_channels_Nbins150_layer8.mat');
fprintf('Load BM data complete.\n');

%% Get capacities.
n = 8;
N = 2^n;

if ~exist('bin_centers', 'var')
    N_bins_each_dim = 128;
    bin_centers = linspace(0, 1, N_bins_each_dim);
end

polarized_channels = synth_ch(bitrevorder(1:N));

C = zeros(N,1);
Pe = zeros(N,1);
for idx = 1:N
    C(idx) = get_I_4D(polarized_channels{idx}, bin_centers);
    %Pe(idx) = get_pe_4D(polarized_channels{idx}, bin_centers);
end

I_total = sum(C);

total_capacity_error = (I_total - N*I)/(N*I);
fprintf('Total capacity error = %.2f %%\n', total_capacity_error*100);

%% Code construction, R=1/2, using BM method.
[~, order] = sort(C, 'descend');
M = 128;
R = M/N;

info_syms = sort(order(1:M), 'ascend').';
info_bits_logical_BM = false(1,N);
info_bits_logical_BM(info_syms) = true;
fprintf('BM Code construction: Complete.\n');

%% BM-PE construction
[~, order] = sort(Pe, 'ascend');
info_syms = sort(order(1:M), 'ascend').';
info_bits_logical_BM_Pe = false(1,N);
info_bits_logical_BM_Pe(info_syms) = true;
fprintf('BM Pe Code construction: Complete.\n');

%% Required simulation. Blackwell-measure construction.
GF_info.m = 2;          % GF(2^2).
GF_info.alpha = 2;
GF_info.mult_table  = [0, 0, 0, 0; 0, 1, 2, 3; 0, 2, 3, 1; 0, 3, 1, 2];
GF_info.add_table   = [0, 1, 2, 3; 1, 0, 3, 2; 2, 3, 0, 1; 3, 2, 1, 0];
% GF_info.kernel_index_vec0 = kernel_index_vec0;          % boolean vector.
% GF_info.kernel_index_vec1 = int32(kernel_index_vec1);   % int32 vector.
% GF_info.kernel_index_mat0 = int32(kernel_index_mat0);   % int32 matrix.
% GF_info.kernel_index_mat1 = int32(kernel_index_mat1);   % int32 matrix.

L = 16;
addpath('sim/');
% if code rate=1/2, in QPSK channel Es/n0 = Eb/n0.
% simulation.

Ebn0_arr = [1.25, 1.5, 1.75, 2.00, 2.25, 2.5, 2.75, 3.00];
min_errors = 1600;
N_Ebn0 = length(Ebn0_arr);
BLERs = zeros(1, N_Ebn0);

for idx = 1:N_Ebn0
    Ebn0 = Ebn0_arr(idx);
    BLERs(idx) = sim_Qary_SCL(N, M, false, info_bits_logical_BM, Ebn0, min_errors, L, GF_info);
end

save('data/bm_construction_q4_L16.mat', 'Ebn0_arr', 'BLERs', 'L', 'info_bits_logical_BM', 'N', 'M');
fprintf('Files saved.\n');

%% N512 binary code.
N = 512;
M = 256;
L = 16;

addpath('codes/polar/GA/');
[channels, ~] = GA(sigma, N);
[~, order_GA] = sort(channels, 'descend');
info_bits = sort(order_GA(1 : M), 'ascend');
info_bits_logical = false(1,N);
info_bits_logical(info_bits) = true;
frozen_bits = ~info_bits_logical;

Ebn0_arr = [1.25, 1.5, 1.75, 2.00, 2.25, 2.5, 2.75, 3.00];
min_errors = 1602;
N_Ebn0 = length(Ebn0_arr);
BLERs = zeros(1, N_Ebn0);

for idx = 1:N_Ebn0
    Ebn0 = Ebn0_arr(idx);
    BLERs(idx) = sim_SCL(N,M,info_bits_logical, Ebn0, min_errors, L);
end

disp('Binary code sim complete.');

%% Use Gaussian-Approximation to construct 4-ary code.
N = 256;
M = 128;
L = 16;

Esn0 = 1.5;             % unit: dB
sigma_GA = 1 * (10^(-Esn0/20));
addpath('codes/polar/GA/');
[channels, ~] = GA(sigma_GA, N);
[~, order_GA] = sort(channels, 'descend');
info_bits = sort(order_GA(1 : M), 'ascend');
info_bits_logical_GA = false(1,N);
info_bits_logical_GA(info_bits) = true;
frozen_bits = ~info_bits_logical_GA;


% if code rate=1/2, in QPSK channel Es/n0 = Eb/n0.
% simulation.
Ebn0_arr = [1.25, 1.5, 1.75, 2.00, 2.25, 2.5, 2.75, 3.00];
min_errors = 1600;
N_Ebn0 = length(Ebn0_arr);
BLERs = zeros(1, N_Ebn0);

for idx = 1:N_Ebn0
    Ebn0 = Ebn0_arr(idx);
    BLERs(idx) = sim_Qary_SCL(N, M, false, info_bits_logical_GA, Ebn0, min_errors, L, GF_info);
end

save('data/ga_construction_q4_L16.mat', 'Ebn0_arr', 'BLERs', 'L', 'info_bits_logical_GA', 'N', 'M');

%% Plot
figure;
load data/bm_construction_q4_L16.mat;
plot(Ebn0_arr, BLERs, 'rs-'); hold on;

load data/ga_construction_q2_L16.mat;
plot(Ebn0_arr, BLERs, 'bs-');

load data/ga_construction_q4_L16.mat;
plot(Ebn0_arr, BLERs, 'gs-');

set(gca, 'yscale', 'log');
xlabel('Eb/n0(dB)');ylabel('BLER'); grid on;
legend('Q4 SCL16 N256 BM', 'Binary SCL16 N512 GA', 'Q4 SCL16 N256 GA');



