%% This file generates the performance curves for CA-SCL FPGA decoders.
addpath('fpga/');
addpath('fpga/fpga_sim/');

addpath('codes/polar/GA/');
addpath('sim/');
N = 32;
M = 16;     % This M contains both info bits and CRC bits.
K_CRC = 4;
K_PCC = 2;

L = 8;

R = M/N;
n = log2(N);

%% Setup code construction.
% Fixed info_bits.
info_bits_logical = logical([0   0   0   0   0   0   0   1   0   0   1   1   1   1   1   1   0   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1]).';
assert(sum(info_bits_logical) == M + K_CRC + K_PCC, 'info_bits_logical mismatch!');

% Setup PCC config structure, used in PCC-polar encoder.
PCC_conf.info_bits_cnt = M + K_CRC;
PCC_conf.parity_bits_cnt = K_PCC;
PCC_conf.N = N;
PCC_conf.nonfrozen_bits_logical = info_bits_logical;
PCC_conf.parity_bits_index = {[0,2,5]+1, [1,6,10]+1};

%% Run simulations!
Ebn0_arr = [1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0];
N_Ebn0 = length(Ebn0_arr);
min_errors = 800;
amp_factor = 3.5;
T = 8;                              % Flipping attempts.
simulated_BLERs = zeros(2, N_Ebn0);

% Using parallel-computing between FPGA and computer.
for sel = 1:2
    switch sel
        case 1
            for ebn0_iter = 1:N_Ebn0
                Ebn0 = Ebn0_arr(ebn0_iter);
                simulated_BLERs(sel, ebn0_iter) = sim_PCC_SCLF_given_construction(K_CRC, PCC_conf, N, M, Ebn0, min_errors, L, T);
            end
        case 2
            for ebn0_iter = 1:N_Ebn0
                Ebn0 = Ebn0_arr(ebn0_iter);
                simulated_BLERs(sel, ebn0_iter) = sim_CA_PCC_SCL_fpga(Ebn0, min_errors, amp_factor);
            end
    end
end

%% Save files.
if T == 0
    save('data/FPGA/fpga_N32_CA_PC_SCL.mat', 'Ebn0_arr', 'info_bits_logical', 'min_errors', 'PCC_conf', 'simulated_BLERs');
else
    save('data/FPGA/fpga_N32_CA_PC_SCLF.mat', 'Ebn0_arr', 'info_bits_logical', 'min_errors', 'PCC_conf', 'simulated_BLERs', 'T');
end

