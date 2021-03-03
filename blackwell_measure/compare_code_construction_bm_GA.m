clear all;
load blackwell_measure/GA_bm.mat;

addpath('sim/');
addpath('codes/polar/');
addpath('D:/iTsinghua/Major/github_syncs/Encoding/PolarCpp/PolarCpp/x64/Release');

%% Using cpp to simulate.
N = 256;
M = 128;
n = log2(N);

Ebn0 = 1.5;
R = M/N;    % coding rate.
sigma_sim = 1/(sqrt(2*R))*10^(-Ebn0/20);

err_cnt = 0;
N_runs = 0;
L = 32;

decoder_config.partially_frozen = false;
decoder_config.is_qary = false;
decoder_config.is_LLR = true;
decoder_config.is_Genie = false;
decoder_config.update_decoder = true;
decoder_config.is_list = true;
decoder_config.L = L;

frozen_bits = ~info_bits_logical_GA;
min_errors = 1600;

%% Sim with SCL decoder.
tic;
while err_cnt < min_errors
    random_bits = rand([1,M])>0.5;
    u = zeros(1,N);
    u(~frozen_bits) = random_bits;
    x = my_polar_encoder(u,n);
    x_bpsk = 1-2*x;
    
    noise_vec = sigma_sim * randn([1,N]);
    y = x_bpsk + noise_vec;   % add AWGN.
    llr = 2*y/(sigma_sim^2);
    
    polar_info_esti = Qary_SC_Decoder(llr, N, 1, frozen_bits, 0, decoder_config);
    if any(polar_info_esti~=random_bits)
        err_cnt = err_cnt + 1;
    end
    N_runs = N_runs + 1;
end

toc;
disp(['Using C++: BLER = ', num2str(err_cnt/N_runs)]);

%% Simulate using CA-SCL decoders.
addpath('sim/');
Ebn0_arr = [1.0, 1.25, 1.5, 1.75, 2.0];

N_Ebn0 = length(Ebn0_arr);
GA_bler = zeros(1, N_Ebn0);
BM_bler = zeros(1, N_Ebn0);

for ebn0_iter = 1:N_Ebn0
    Ebn0 = Ebn0_arr(ebn0_iter);
    GA_bler(ebn0_iter) = sim_CASCL_given_construction(8, N, M, info_bits_logical_GA,Ebn0, min_errors, 16);
    BM_bler(ebn0_iter) = sim_CASCL_given_construction(8, N, M, info_bits_logical_ML,Ebn0, min_errors, 16);
end




