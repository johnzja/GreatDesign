%% This file generates the performance curves for CA-SCL FPGA decoders.
addpath('fpga/');
addpath('fpga/fpga_sim/');

addpath('codes/polar/GA/');
addpath('sim/');
N = 32;
M = 16;     % This M contains both info bits and CRC bits.
K_CRC = 4;
L = 8;

R = M/N;
n = log2(N);

%% Call Gaussian-Approximation code construction.
design_Ebn0 = 2.5;
sigma_cc = 1/sqrt(2 * R) * 10^(-design_Ebn0/20);
[channels, ~] = GA(sigma_cc, N);  

[~, channel_ordered] = sort(channels, 'descend');
info_bits = sort(channel_ordered(1 : M+K_CRC), 'ascend');  
frozen_bits = ones(N , 1);                             
frozen_bits(info_bits) = 0;
info_bits_logical = logical(mod(frozen_bits + 1, 2));

%% Run simulations!
Ebn0_arr = [1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0];
N_Ebn0 = length(Ebn0_arr);
min_errors = 800;
amp_factor = 3.5;

simulated_BLERs = zeros(2, N_Ebn0);

% Using parallel-computing between FPGA and computer.
for sel = 1:2
    switch sel
        case 1
            for ebn0_iter = 1:N_Ebn0
                Ebn0 = Ebn0_arr(ebn0_iter);
                simulated_BLERs(sel, ebn0_iter) = sim_CASCL_given_construction(K_CRC, N, M, info_bits_logical, Ebn0, min_errors, L);
            end
        case 2
            for ebn0_iter = 1:N_Ebn0
                Ebn0 = Ebn0_arr(ebn0_iter);
                simulated_BLERs(sel, ebn0_iter) = sim_CA_SCL_fpga(Ebn0, min_errors, amp_factor);
            end
    end
end

%% Save files.
save('data/FPGA/fpga_N32_first_test.mat');

