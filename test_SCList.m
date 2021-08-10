%% Setup paths.
addpath('codes/polar/');
addpath('D:/iTsinghua/Major/github_syncs/Encoding/PolarCpp/PolarCpp/x64/Release');

%% TEST
N = 256;
M = 128;
frozen_bits = [1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0];


n = log2(N);
lambda_offset = 2.^(0:n);
llr_layer_vec = get_llr_layer(N);
bit_layer_vec = get_bit_layer(N);

Ebn0 = 2.5;
R = M/N;    % coding rate.
sigma_sim = 1/(sqrt(2*R))*10^(-Ebn0/20);

min_errors = 1600;
err_cnt = 0;
N_runs = 0;
L = 16;

tic;
while err_cnt < min_errors
    random_bits = rand([1,M])>0.5;
    u = zeros(1,N);
    u(info_bits_logical_bm) = random_bits;
    x = my_polar_encoder(u,n);
    x_bpsk = 1-2*x;
    
    noise_vec = sigma_sim * randn([1,N]);
    y = x_bpsk + noise_vec;   % add AWGN.
    llr = 2*y/(sigma_sim^2);
    
    polar_info_esti = SCL_decoder(llr,L,M,frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec);
    if any(polar_info_esti~=random_bits.')
        err_cnt = err_cnt + 1;
    end
    N_runs = N_runs + 1;
end

toc;
disp(['Using MATLAB: BLER = ', num2str(err_cnt/N_runs)]);

%% using cpp.
N = 512;
M = 256;

addpath('codes/polar/GA/');
[channels, ~] = GA(sigma, N);  

[~, order_GA] = sort(channels, 'descend');
info_bits = sort(order_GA(1 : M), 'ascend'); 
info_bits_logical = false(1,N);
info_bits_logical(info_bits) = true;
frozen_bits = ~info_bits_logical;

err_cnt = 0;
N_runs = 0;

decoder_config.partially_frozen = false;
decoder_config.is_qary = false;
decoder_config.is_LLR = true;
decoder_config.is_Genie = false;
decoder_config.update_decoder = true;
decoder_config.is_list = true;
decoder_config.L = L;

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