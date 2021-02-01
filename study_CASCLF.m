clear;clc;
addpath('sim/');

%% Simulate CRC-PCC-SCLF code.
% Compare: K_CRC=4, K_PCC=8 code+SCLF v.s. K_PCC=12 code.
N = 512;
M = 256;
K_CRC = 8;


% Ebn0_arr = 1.0:0.2:2.4;
Ebn0_arr = [2.25, 2.5, 2.75, 3.00];           % All under this Eb/n0.
min_errors = 1200;
N_ebn0 = length(Ebn0_arr);
L = 16;
T = 32;

BLERs = zeros(2, N_ebn0);

%% Simulation Loop.
for ebn0_iter = 1:N_ebn0
    Ebn0 = Ebn0_arr(ebn0_iter);
    BLERs(1, ebn0_iter) = sim_CA_SCLF(K_CRC,N, M, Ebn0, min_errors, L, 0);
    BLERs(2, ebn0_iter) = sim_CA_SCLF(K_CRC,N, M, Ebn0, min_errors, L, T);
end

save('data/bler_crc8_cascl_t0_t32_ext.mat');

