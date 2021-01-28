clear;
addpath('codes/');
addpath('codes/polar/');
addpath('codes/polar/GA/');



%% Setup Monte-Carlo simulation parameters.
n = 9;
N = 2^n;    % N=512 polar code.
M = 256;    % info bits.
K = 20;     % parity bits,


lambda_offset = 2.^(0 : n);
llr_layer_vec = get_llr_layer(N);
bit_layer_vec = get_bit_layer(N);

%% Construct Polar code with GA algorithm.
design_Ebn0 = 2.5;
R=(M+K)/N;
sigma_cc = 1/sqrt(2 * R) * 10^(-design_Ebn0/20);
[channels, ~] = GA(sigma_cc, N);                    
[~, channel_ordered] = sort(channels, 'descend');
info_bits = sort(channel_ordered(1 : K+M), 'ascend');   % Both real "info" bits and parity bits are all treated as info bits during construction.
frozen_bits = ones(N , 1);                              % column vector(?)
frozen_bits(info_bits) = 0;
info_bits_logical = logical(mod(frozen_bits + 1, 2));

%% Perform Genie-aided SC decoding to estimate channel reliability.
cnt_err_unfrozen = zeros(1, M+K);
cnt_min_errors = 10;
cnt_total_decoding_attempts = 0;
Ebn0 = 1.5;     % in dB
sigma_sim = 1/(sqrt(2*R)) * 10^(-Ebn0/20);

last_min_errors = 0;
min_errors = 0;
while(min_errors < cnt_min_errors && cnt_total_decoding_attempts < 10000000)
    random_bits = (rand(1, M+K) < 0.5);
    u = zeros(1, N);    
    u(info_bits_logical) = random_bits;
    x = my_polar_encoder(u, n);
    x_BPSK = 1-2*x;
    y = x_BPSK + sigma_sim * randn(size(x_BPSK));
    llr = 2*y/(sigma_sim^2);
    error_pattern = Genie_aided_SC_decoder(llr, K+M, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec, random_bits);
    cnt_err_unfrozen = cnt_err_unfrozen + error_pattern;
    cnt_total_decoding_attempts = cnt_total_decoding_attempts + 1;
    min_errors = min(cnt_err_unfrozen);
    if min_errors > last_min_errors
        last_min_errors = min_errors;
        disp(min_errors);
    end
end

disp('finished');
