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
K_SCLF_PCC = 8;
K_PCC = 12;

design_Ebn0 = 1.5;
[PCC_structs, ~] = get_standard_PCC(N,M,design_Ebn0);

PCC = PCC_structs(K_PCC);

% Ebn0_arr = 1.0:0.2:2.4;
Ebn0_arr = [1.5];           % All under this Eb/n0.
min_errors = 1200;
N_ebn0 = length(Ebn0_arr);
L = 16;
T = 32;

BLERs = zeros(3, N_ebn0);
K_SCLF_CRC_arr = [4, 6, 8, 10, 12];
K_SCLF_PCC_arr = (4:1:16);

%% Simulation Loop.
% parfor sim_setup_iter = 1:3
%     for ebn0_iter = 1:N_ebn0
%         Ebn0 = Ebn0_arr(ebn0_iter);
%         switch sim_setup_iter
%             case 1
%                 BLERs(sim_setup_iter,ebn0_iter) = sim_PCC(PCC, Ebn0, min_errors, L);
%             case 2
%                 BLERs(sim_setup_iter,ebn0_iter) = sim_PCC_SCLF(K_SCLF_CRC, K_SCLF_PCC, N,M,Ebn0, min_errors, L, 0);
%             case 3
%                 BLERs(sim_setup_iter,ebn0_iter) = sim_PCC_SCLF(K_SCLF_CRC, K_SCLF_PCC, N,M,Ebn0, min_errors, L, T);
%         end
%     end
% end

% simulate without bit-flipping.
% using CRC-PCC-Polar concatenation scheme.
N_CRC_arr = length(K_SCLF_CRC_arr);
N_PCC_arr = length(K_SCLF_PCC_arr);

BLERs_mat = zeros(N_CRC_arr, N_PCC_arr);
Ebn0 = Ebn0_arr(1);

for k_crc_iter = 1:N_CRC_arr
    for k_pcc_iter = 1:N_PCC_arr
        K_SCLF_CRC = K_SCLF_CRC_arr(k_crc_iter);
        K_SCLF_PCC = K_SCLF_PCC_arr(k_pcc_iter);
        BLERs(k_crc_iter, k_pcc_iter) = sim_PCC_SCLF(K_SCLF_CRC, K_SCLF_PCC, N, M, Ebn0, min_errors, L, 0);
        % do not invoke bit-flip now.
    end
end




