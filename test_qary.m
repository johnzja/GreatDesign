addpath('codes/');
addpath('codes/polar/');
addpath('codes/polar/GA/');
addpath('sim/');
addpath('D:/iTsinghua/Major/github_syncs/Encoding/PolarCpp/PolarCpp/x64/Debug');

clear all;

%% Test binary polar code.
N = 8;
M = 4;
m = 1;
n = log2(N);
R = M/N;
frozen_bits = [1, 1, 1, 0, 1, 0, 0, 0];
Ebn0 = 5;

sigma_sim = 1/sqrt(2*R) * (10^(-Ebn0/20));

N_sim = 100;

decoder_config.partially_frozen = false;
decoder_config.is_qary = false;
decoder_config.is_LLR = true;
decoder_config.is_Genie = false;
decoder_config.update_decoder = true;   % it can be automatically modified into false.
decoder_config.is_list = true;
%decoder_config.L = 4;


for sim_iter = 1:N_sim
    random_bits = (rand([1,M])>0.5);
    u = zeros(1,N);
    u(~frozen_bits) = random_bits;
    x = my_polar_encoder(u,n);
    x_bpsk = 1-2*x;
    
    y = x_bpsk + sigma_sim * randn(size(x_bpsk));
    llr = 2*y/(sigma_sim^2);
    
    % Call C++ decoder.
    polar_info_esti = Qary_SC_Decoder(llr, N, m, frozen_bits, 0, decoder_config);
end
