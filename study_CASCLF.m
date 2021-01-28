
clear;clc;
addpath('sim/');

%% Simulate CRC-PCC-SCLF code.
% Compare: K_CRC=4, K_PCC=8 code+SCLF v.s. K_PCC=12 code.
N = 512;
M = 256;
K_CRC = 8;


% Ebn0_arr = 1.0:0.2:2.4;
Ebn0_arr = [1, 1.25, 1.5, 1.75, 2.0];           % All under this Eb/n0.
min_errors = 1200;
N_ebn0 = length(Ebn0_arr);
L = 16;
T = 0;

BLERs = zeros(1, N_ebn0);



%% Simulation Loop.
for ebn0_iter = 1:N_ebn0
    Ebn0 = Ebn0_arr(ebn0_iter);
    BLERs(ebn0_iter) = sim_CA_SCLF(K_CRC,N, M, Ebn0, min_errors, L, T);
end

save('data/bler_crc8_cascl.mat');

