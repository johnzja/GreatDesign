% * **************************************************************************
% * ********************                                  ********************
% * ********************      Tsinghua University         ********************
% * ********************                                  ********************
% * **************************************************************************
% *                                                                          *
% *                                   _oo8oo_                                *
% *                                  o8888888o                               *
% *                                  88" . "88                               *
% *                                  (| -_- |)                               *
% *                                  0\  =  /0                               *
% *                                ___/'==='\___                             *
% *                              .' \\|     |// '.                           *
% *                             / \\|||  :  |||// \                          *
% *                            / _||||| -:- |||||_ \                         *
% *                           |   | \\\  -  /// |   |                        *
% *                           | \_|  ''\---/''  |_/ |                        *
% *                           \  .-\__  '-'  __/-.  /                        *
% *                         ___'. .'  /--.--\  '. .'___                      *
% *                      ."" '<  '.___\_<|>_/___.'  >' "".                   *
% *                     | | :  `- \`.:`\ _ /`:.`/ -`  : | |                  *
% *                     \  \ `-.   \_ __\ /__ _/   .-` /  /                  *
% *                 =====`-.____`.___ \_____/ ___.`____.-`=====              *
% *                                   `=---=`                                *
% * **************************************************************************
% * **************************************************************************
% * *******************      				                ******************
% * ********************         Wish to Graduate          *******************
% * *********************                                 ********************
% * **************************************************************************

clear all; close all; clc;

set(0,'DefaultLineMarkerSize',4);
set(0,'DefaultTextFontSize',14);

set(0,'DefaultAxesFontSize',12);
set(0,'DefaultLineLineWidth',1.4);

%% Step1: Plot curves for SCLFlip decoders.
set(0,'defaultfigurecolor','w');
figure('color',[1 1 1]);
hold on;
% Plot the performance of pure Parity-check-concatenated Polar codes.
load data/compare_flip_with_others/bler_pcc12_t0.mat;
plot(Ebn0_arr, BLERs, 'color', [0, 0.5, 1], 'LineStyle', '-', 'marker', 's', 'markersize', 6);

% Without PCC: CASCL v.s. CASCL+Flip(T=32), all using CRC8.
load data/compare_flip_with_others/bler_crc8_casclf_t0_t32.mat;
plot(Ebn0_arr, blers_crc8_cascl, 'color', [0, 0.5, 0], 'LineStyle', '-.', 'marker', 's', 'markersize', 6);
plot(Ebn0_arr, blers_crc8_casclf_t32, 'color', [0, 0.5, 0], 'LineStyle', '-', 'marker', 's', 'markersize', 6);

% Plot CRC8-PC8 code under SCL without flipping.
load data/compare_flip_with_others/bler_crc8_pcc8_t0.mat;
plot(Ebn0_arr, BLERs, 'color', [1, 0, 1], 'LineStyle', '-.', 'marker', 'd', 'markersize', 6);

% CRC8-PC4-polar code with flip(T=32).
load data/compare_flip_with_others/bler_sclf_crc8_pcc4.mat;
plot(Ebn0_arr, BLERs, 'color', [0, 0, 1], 'LineStyle', '-', 'marker', 'x', 'markersize', 6);

% Plot CRC8-PC8 code under SCLFlip.
load data/compare_flip_with_others/bler_crc8_pcc8_t32.mat
plot(Ebn0_arr, BLERs, 'rx-');

grid on; box on;
set(gca, 'yscale', 'log');
xlabel('E_b/n_0(dB)');
ylabel('BLER');


legend('PCC12-SCL', 'CRC8-SCL', 'CRC8-SCLF(T=32)', 'CRC8-PC8-SCL', 'CRC8-PC4-SCLF(T=32)', 'CRC8-PC8-SCLF(T=32)');
set(gca,'FontName','Times New Roman');

%% Plot2: plot the code performance with K_CRC and K_PCC.
% 2-Dim plot.

load data/compare_flip_with_others/CRC-PCC_performance_test.mat;
BLERs_nf = BLERs;
figure(2);hold on;

for idx = 1:5
    plot(K_SCLF_PCC_arr, BLERs_nf(idx, :), 'x-');
end
legend('CRC4', 'CRC6', 'CRC8', 'CRC10', 'CRC12');
xlabel('K\_PCC');ylabel('BLER');
set(gca, 'yscale', 'log');
set(gca,'FontName','Times New Roman');
grid on; box on;

load data/compare_flip_with_others/CRC-PCC-SCLF_performance_test.mat
BLERs_f = BLERs;
figure(3);hold on;
for idx = 1:5
    plot(K_SCLF_PCC_arr, BLERs_f(idx, :), 'x-');
end
legend('CRC4', 'CRC6', 'CRC8', 'CRC10', 'CRC12');
xlabel('K\_PCC');ylabel('BLER');
set(gca, 'yscale', 'log');
set(gca,'FontName','Times New Roman');
grid on; box on;


%% Plot3: FPGA code performance
figure(4); hold on; grid on;
set(gca, 'yscale', 'log');
load data/FPGA/fpga_N32_first_test.mat;
plot(Ebn0_arr, simulated_BLERs(1,:),  'color', [0, 0.75, 0.9], 'LineStyle', '-', 'marker', 'x', 'markersize', 8); % CA-SCL
plot(Ebn0_arr, simulated_BLERs(2,:), 'color', [0, 0.75, 0.9], 'LineStyle', '-.', 'marker', 's', 'markersize', 8); % CA-SCL

load data/FPGA/fpga_N32_CA_PC_SCL.mat;
plot(Ebn0_arr, simulated_BLERs(1,:),'color', [0.6, 0, 1], 'LineStyle', '-', 'marker', 'x', 'markersize', 8);    % CA-PC-SCL N=32 R=1/2 PC:hand-constructed simulated.
plot(Ebn0_arr, simulated_BLERs(2,:),'color', [0.6, 0, 1], 'LineStyle', '-.', 'marker', 's', 'markersize', 8);   % CA-PC-SCL N=32 R=1/2 PC:hand-constructed FPGA.

load data/FPGA/fpga_N32_CA_PC_SCLF.mat;
plot(Ebn0_arr, simulated_BLERs(1,:), 'color', [1, 0, 0], 'LineStyle', '-', 'marker', 'x', 'markersize', 8);     % CA-PC-SCLF sim
plot(Ebn0_arr, simulated_BLERs(2,:), 'color', [1, 0, 0], 'LineStyle', '-.', 'marker', 's', 'markersize', 8);    % CA-PC-SCLF FPGA
legend('CA-SCL sim', 'CA-SCL FPGA', 'CA-PC-SCL sim', 'CA-PC-SCL FPGA', 'CA-PC-SCLF sim', ...
    'CA-PC-SCLF FPGA');
xlabel('E_b/n_0(dB)'); ylabel('BLER');
set(gca,'FontName','Times New Roman');
box on;

%% Plot4: Burst-error segment and standard PCC construction.
N = 512; M = 256;
addpath('codes/');
addpath('codes/polar/');
addpath('codes/polar/GA/');

R=M/N;
sigma_cc = 1/sqrt(2 * R) * 10^(-0/20);
[channels, ~] = GA(sigma_cc, N);  

[~, channel_ordered] = sort(channels, 'descend');
info_bits = sort(channel_ordered(1:M), 'ascend');  
frozen_bits = ones(N , 1);                             
frozen_bits(info_bits) = 0;
info_bits_logical = logical(mod(frozen_bits + 1, 2));

unfrozen_channels = channels(info_bits_logical);

% Draw scatters.
figure('color', [1,1,1]);
h = scatter(1:M, 1./(1+exp(unfrozen_channels)), 'b*');
set(gca, 'yscale', 'log');
set(gca,'FontName','Times New Roman');
ylabel('Error probability'); xlabel('Unfrozen bit channel index');
box on;

% Calculate burst-error segments.
K = 8;
must_use = (K+1)*K/2;
remain = M+K-must_use;
distribute = floor(remain/K);

burst_err_seg_index = cell(K,1); 
parity_bits_index = cell(K,1);
s = 1;
for k = 1:K-1
   burst_err_seg_index{k} = (s:s+K-k+distribute);
   s = s+K-k+1+distribute;
end
burst_err_seg_index{K} = (s:M+K);

% Draw vertical dashed lines.
% figure(6);
set(gca, 'xlim', [0, M+1]);
ylim=get(gca,'Ylim');
hold on;
for k = 1:K-1
    s = burst_err_seg_index{k}(end);
    plot([s,s],ylim,'m--'); 
end

%% Plot5: Complexity of SCLF decoders.
load data/test_Ebn0_3.00.mat;
bler_3 = BLERs; tr_3 = Trial_Rates;
load data/flip_mean_complexity.mat;
BLERs(end) = bler_3;
Trial_Rates(end) = tr_3;

% Plot this curve.
N_Ebn0 = length(BLERs);
assert(N_Ebn0 == length(Trial_Rates));
figure('color', [1,1,1]); hold on;

plot(1:0.25:3, Trial_Rates, 'LineStyle', '-', 'marker', 'x', 'markersize', 8);
plot(1:0.25:3, ones(1, N_Ebn0), 'LineStyle', '-', 'marker', 's', 'markersize', 8);

set(gca, 'ylim', [0, max(Trial_Rates)]);
grid on; xlabel('E_b/n_0 (dB)');ylabel('# Trials of SCL-decoding');
legend('SCLF decoder', 'SCL decoder' );


