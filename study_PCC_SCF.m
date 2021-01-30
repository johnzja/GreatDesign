% Study this problem: is SC-Flip effective in decoding PCC-Polar
% codes? We should first approach the answer by Monte-Carlo simulations on
% the error distribution within one PC equation. If a partially decoded
% sequence is obtained, which ends with an erroneous parity bit, is the
% uncorrectly decoded bit identical to the bit with the most "uncertain"
% LLR?

% Three possibilities: identical, non-identical, >=3 (odd number) errors.
% construct SCFlip decoder.
clear;clc;
addpath('sim/');

%% Simulate CRC-PCC-SCLF code.
% Compare: K_CRC=4, K_PCC=8 code+SCLF v.s. K_PCC=12 code.
N = 512;
M = 256;
K_SCLF_CRC = 8;
K_SCLF_PCC = 4;
K_PCC = 12;

design_Ebn0 = 1.5;
[PCC_structs, ~] = get_standard_PCC(N,M,design_Ebn0);

PCC = PCC_structs(K_PCC);

% Ebn0_arr = 1.0:0.2:2.4;
Ebn0_arr = [2.25, 2.5, 2.75];           % All under this Eb/n0.
min_errors = 1200;
N_ebn0 = length(Ebn0_arr);
L = 16;
T = 32;

BLERs = zeros(1, N_ebn0);


%% Simulation Loop.
for ebn0_iter = 1:N_ebn0
    Ebn0 = Ebn0_arr(ebn0_iter);
    BLERs(ebn0_iter) = sim_PCC_SCLF(K_SCLF_CRC, K_SCLF_PCC, N,M,Ebn0, min_errors, L, T);
end

save('data/test.mat');
