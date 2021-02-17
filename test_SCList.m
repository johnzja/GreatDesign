addpath('codes/polar/');

N = 8;
M = 4;
frozen_bits = [1, 1, 1, 0, 1, 0, 0, 0];


    
n = log2(N);
lambda_offset = 2.^(0:n);
llr_layer_vec = get_llr_layer(N);
bit_layer_vec = get_bit_layer(N);

Ebn0 = 2.5;
R = M/N;    % code rate.
sigma_sim = 1/(sqrt(2*R))*10^(-Ebn0/20);

min_errors = 800;
err_cnt = 0;
N_runs = 0;
L = 4;

while err_cnt < min_errors
    random_bits = rand([1,M])>0.5;
    u = zeros(1,N);
    u(~frozen_bits) = random_bits;
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

disp(['BLER = ', num2str(err_cnt/N_runs)]);