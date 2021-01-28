%% step1: Plot curves for sclf.
load data/bler_pccscl_sclf_crc8_pcc8.mat;

figure(1);
hold on;
plot(Ebn0_arr, BLERs(1,:), 'bx-');
plot(Ebn0_arr, BLERs(2,:), 'gx-');
plot(Ebn0_arr, BLERs(3,:), 'rx-');
grid on;
set(gca, 'yscale', 'log');
xlabel('Eb/n_0');
ylabel('BLER');

load data/bler_sclf_crc8_pcc4.mat;
plot(Ebn0_arr, BLERs, 'cx-');
legend('PCC12-SCL', 'CRC8-PCC8-SCL', 'CRC8-PCC8-SCLF(T=32)', 'CRC8-PCC4-SCLF(T=32)');

%% Plot2: plot the code performance with K_CRC and K_PCC.
% 2-Dim plot.

load data/CRC-PCC_performance_test.mat;
BLERs_nf = BLERs;
figure(2);hold on;
for idx = 1:5
    plot(K_SCLF_PCC_arr, BLERs_nf(idx, :), 'x-');
end
legend('CRC4', 'CRC6', 'CRC8', 'CRC10', 'CRC12');
xlabel('K\_PCC');ylabel('BLER');
set(gca, 'yscale', 'log');
grid on;

load data/CRC-PCC-SCLF_performance_test.mat
BLERs_f = BLERs;
figure(3);hold on;
for idx = 1:5
    plot(K_SCLF_PCC_arr, BLERs_f(idx, :), 'x-');
end
legend('CRC4', 'CRC6', 'CRC8', 'CRC10', 'CRC12');
xlabel('K\_PCC');ylabel('BLER');
set(gca, 'yscale', 'log');
grid on;