clear;clc;
load data/PCC_config.mat;
addpath('codes/');
addpath('codes/polar/');
addpath('codes/polar/GA/');

%% Sim parameters.
N = 512;    % transmitted bits.
n = log2(N);
M = 256;    % info bits.
L = 32;

N_sim = 1000;

R = M/N;        % overall code rate.
Ebn0 = 1.5;     % in dB
sigma_sim = 1/(sqrt(2*R)) * 10^(-Ebn0/20);

%% Calculate some handy variables.
lambda_offset = 2.^(0 : n);
llr_layer_vec = get_llr_layer(N);
bit_layer_vec = get_bit_layer(N);

%% Start Simulation

for sim_iter = 1:N_sim
    random_bits = (rand([1,M])>0.5);
    encoded_bits = PCC_polar_encoder(random_bits, PCC);
    x_BPSK = 1-2*encoded_bits;
    y = x_BPSK + sigma_sim * randn(size(x_BPSK));   % add AWGN.
    llr = 2*y/(sigma_sim^2);
    
    % PCC-SCL decoder.
    polar_info_esti = PCC_SCL_decoder(llr,L,PCC,lambda_offset,llr_layer_vec,bit_layer_vec);
    % Check correctness of PCC-SCL decoding result.
    err_cnt = sum(xor(random_bits, polar_info_esti));
    disp(err_cnt);
end


% SC decoding.
% nonfrozen_bits_logical = PCC.nonfrozen_bits_logical;




