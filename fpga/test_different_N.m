clear;clc;
addpath('codes/');
addpath('codes/polar/');
addpath('codes/polar/GA/');
addpath('sim/');
N = 32;
M = 22;

R = M/N;
n = log2(N);

%% Call Gaussian-Approximation code construction.
design_Ebn0 = 1.5;
sigma_cc = 1/sqrt(2 * R) * 10^(-design_Ebn0/20);
[channels, ~] = GA(sigma_cc, N);  

[~, channel_ordered] = sort(channels, 'descend');
info_bits = sort(channel_ordered(1 : M), 'ascend');  
frozen_bits = ones(N , 1);                             
frozen_bits(info_bits) = 0;
info_bits_logical = logical(mod(frozen_bits + 1, 2));

% Generate the frozen_bits sequence used in Verilog.
frozen_str = sprintf("%d'b", N);
for i=N:-1:1
    frozen_str = strcat(frozen_str, char('0'+frozen_bits(i)));
end

disp(frozen_str);