% Study this problem: is SC-Flip effective in decoding PCC-Polar
% codes? We should first approach the answer by Monte-Carlo simulations on
% the error distribution within one PC equation. If a partially decoded
% sequence is obtained, which ends with an erroneous parity bit, is the
% uncorrectly decoded bit identical to the bit with the most "uncertain"
% LLR?

% Three possibilities: identical, non-identical, >=3 (odd number) errors.
% construct SCFlip decoder.
clear;clc;
addpath('codes/');
addpath('codes/polar/');
addpath('codes/polar/GA/');

%% Step1: Construct PCC-CRC-Polar codes.
K_CRC = 4;
K_PCC = 8;
N = 512;        % mid-short polar code.
M = 256;
L = 16;         % decode list size.

% Get CRC object.
[gen, det, g] = get_crc_objective(K_CRC);
n = log2(N);

[G_crc, H_crc] = crc_generator_matrix(g, M);
CRC_Qmatrix = G_crc(:, M+1 : end);        % Generator matrix for CRC bits.; G=[I,Q].

% Calculate some handy variables used in general SCL decoder.
lambda_offset = 2.^(0:n);
llr_layer_vec = get_llr_layer(N);
bit_layer_vec = get_bit_layer(N);
R = M/N;        % code rate.

% Concatenation scheme: info bits => [info bits, CRC] => [info bits, CRC, PCC] => x_1^N, polar codeword.
design_Ebn0 = 1.5;
[PCC_structs, GA_ch_reliability] = get_standard_PCC(N, M+K_CRC, design_Ebn0);

PCC_conf = PCC_structs(K_PCC);                   % Select the object with K_PCC additional parity bits.
PCC_CRC_polar_config.PCC_conf = PCC_conf;
PCC_CRC_polar_config.N = N;
PCC_CRC_polar_config.M = M;                 % cnt. info bits.
PCC_CRC_polar_config.H_CRC = H_crc;         % CRC-Check matrix.
PCC_CRC_polar_config.Q_CRC = CRC_Qmatrix;   % Redundancy Generator matrix.
PCC_CRC_polar_config.K_CRC = K_CRC;         % cnt. Additional CRC bits.


%% Construct SRCS set according to the constructed CRC-PCC-Polar code.
nbl = PCC_conf.nonfrozen_bits_logical;
CS = nbl;   % Critical set.
for k = 1:n
    step = 2^k;
    cnt_iter = N/step;
    for iter = 1:cnt_iter
        start_index = (iter-1)*step+1;
        stop_index = start_index+step-1;
        if all(nbl(start_index:stop_index))
            CS(start_index:stop_index) = zeros(step,1);
            CS(start_index) = 1;
        end
    end
end

% Sort CS according to GA-estimated channel reliability.
tmp = (1:N).';
CS = tmp(CS);
CS_reliability = GA_ch_reliability(CS);
[~, CS_order] = sort(CS_reliability, 'ascend');
CS = CS(CS_order);

% Construct the information structure needed for the PCC-CRC-Polar decoder.
PCC_CRC_polar_decoder_info.L = L;
PCC_CRC_polar_decoder_info.lambda_offset = lambda_offset;
PCC_CRC_polar_decoder_info.llr_layer_vec = llr_layer_vec;
PCC_CRC_polar_decoder_info.bit_layer_vec = bit_layer_vec;
PCC_CRC_polar_decoder_info.CS = CS;
PCC_CRC_polar_decoder_info.T = 32;  % Trials = 32.

%% CRC-PCC-Polar Encoding.
loops = 0;
total_error_cnt = 0;
first_error_pattern = zeros(1,M);

while total_error_cnt < 10000
    random_bits = (rand([1,M])>0.5);
    CRC_aided_bits = [random_bits, logical(mod(random_bits*PCC_CRC_polar_config.Q_CRC, 2))];  % row vector.
    PCC_CRC_encoded_bits = PCC_polar_encoder(CRC_aided_bits, PCC_CRC_polar_config.PCC_conf);

    %% Channel.
    Ebn0 = 1.4;
    x_PCC_CRC = 1-2*PCC_CRC_encoded_bits;

    sigma = 1/(sqrt(2*R))*10^(-Ebn0/20);
    noise_vec = sigma * randn([1,N]);
    y_PCC_CRC = x_PCC_CRC + noise_vec;

    llr_PCC_CRC = 2*y_PCC_CRC/(sigma^2);


    %% CRC-PCC-Polar Decoding.
    polar_info_esti_PCC_CRC = PCC_CRC_SCLF_decoder(llr_PCC_CRC, PCC_CRC_polar_config, PCC_CRC_polar_decoder_info);
    loops = loops + 1;
    errs = xor(polar_info_esti_PCC_CRC, random_bits);
    if any(errs)
        ind = find(errs,1);
        first_error_pattern(ind) = first_error_pattern(ind) + 1;
        total_error_cnt = total_error_cnt + 1;
    end
end



disp('Complete');


