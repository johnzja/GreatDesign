clear all;
%% step1: Plot curves for sclf.
figure(1);
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
xlabel('Eb/n_0');
ylabel('BLER');

% CRC8-PC4-polar code with flip(T=32).
load data/compare_flip_with_others/bler_sclf_crc8_pcc4.mat;
plot(Ebn0_arr, BLERs, 'bx-');


% Without PCC: CASCL v.s. CASCL+Flip(T=32), all using CRC8.
load data/compare_flip_with_others/bler_crc8_casclf_t0_t32.mat;
plot(Ebn0_arr, blers_crc8_casclf_t32, 'gx-');
plot(Ebn0_arr, blers_crc8_cascl, 'gs-');


legend('PCC12-SCL', 'CRC8-PCC8-SCL', 'CRC8-PCC8-SCLF(T=32)', 'CRC8-PCC4-SCLF(T=32)', 'CRC8-SCLF(T=32)', 'CRC8-SCL');

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
grid on;