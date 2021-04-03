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
plot(Ebn0_arr, BLERs, 'kx-');

% Plot CRC8-PC8 code under SCL without flipping.
load data/compare_flip_with_others/bler_crc8_pcc8_t0.mat;
plot(Ebn0_arr, BLERs, 'mx-');

% Plot CRC8-PC8 code under SCLFlip.
load data/compare_flip_with_others/bler_crc8_pcc8_t32.mat
plot(Ebn0_arr, BLERs, 'rx-');

grid on;
set(gca, 'yscale', 'log');
xlabel('Eb/n_0(dB)');
ylabel('BLER');

% CRC8-PC4-polar code with flip(T=32).
load data/compare_flip_with_others/bler_sclf_crc8_pcc4.mat;
plot(Ebn0_arr, BLERs, 'bx-');


% Without PCC: CASCL v.s. CASCL+Flip(T=32), all using CRC8.
load data/compare_flip_with_others/bler_crc8_casclf_t0_t32.mat;
plot(Ebn0_arr, blers_crc8_casclf_t32, 'gx-');
plot(Ebn0_arr, blers_crc8_cascl, 'gs-');

legend('PCC12-SCL', 'CRC8-PC8-SCL', 'CRC8-PC8-SCLF(T=32)', 'CRC8-PC4-SCLF(T=32)', 'CRC8-SCLF(T=32)', 'CRC8-SCL');
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
grid on;

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
grid on;


%% Plot3: FPGA code performance
figure(4); hold on; grid on;
set(gca, 'yscale', 'log');
load data/FPGA/fpga_N32_first_test.mat;
plot(Ebn0_arr, simulated_BLERs(1,:), 'rx-');    % CA-SCL N=32 R=1/2 simulated.
plot(Ebn0_arr, simulated_BLERs(2,:), 'ro-');    % CA-SCL N=32 R=1/2 FPGA.

load data/FPGA/fpga_N32_CA_PC_SCL.mat;
plot(Ebn0_arr, simulated_BLERs(1,:), 'bx-');    % CA-PC-SCL N=32 R=1/2 PC:hand-constructed simulated.
plot(Ebn0_arr, simulated_BLERs(2,:), 'bo-');    % CA-PC-SCL N=32 R=1/2 PC:hand-constructed FPGA.

load data/FPGA/fpga_N32_CA_PC_SCLF.mat;
plot(Ebn0_arr, simulated_BLERs(1,:), 'gx-');
plot(Ebn0_arr, simulated_BLERs(2,:), 'go-');
legend('CA-SCL sim', 'CA-SCL FPGA', 'CA-PC-SCL sim', 'CA-PC-SCL FPGA', 'CA-PC-SCLF sim', ...
    'CA-PC-SCLF FPGA');
xlabel('Eb/n0(dB)'); ylabel('BLER');
set(gca,'FontName','Times New Roman');
