clear;clc;
load data/PCC_config_K18.mat;
addpath('codes/');
addpath('codes/polar/');
addpath('codes/polar/GA/');

%% Get CRC8 objects.
crc_length = 8;
[gen, det, g] = get_crc_objective(crc_length);

%% Sim parameters & til variables.
N = 512;        % cnt. of transmitted bits.
n = log2(N);
M = 256;        % info bits.
L = 32;         % SCL decoder list length.

CASCL_K = M+crc_length;

[G_crc, H_crc] = crc_generator_matrix(g, CASCL_K - crc_length);
crc_parity_check = G_crc(:, CASCL_K - crc_length + 1 : end);

N_sim = 1000;
R = M/N;        % overall code rate.
Ebn0_arr = [-0.5, 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 2];
N_Ebn0 = length(Ebn0_arr);
BLER_CRC = zeros(1,N_Ebn0);
BLER_PCC = zeros(1,N_Ebn0);

min_error = 800;            % 1600 observations of error needed.
% min_errors is determined by a Bernoulli porbabilistic model with a
% confidence of 95%. This number can be adjusted smaller to ensure quick 
% simulation.
% Formula: min_error >= (z_0.975/required_maximum_relative_error)^2.

%% Calculate some handy variables.
lambda_offset = 2.^(0:n);
llr_layer_vec = get_llr_layer(N);
bit_layer_vec = get_bit_layer(N);
% Re-construct CRC-Polar code.
design_Ebn0 = 1.5;
R=CASCL_K/N;
sigma_cc = 1/sqrt(2 * R) * 10^(-design_Ebn0/20);
[channels, ~] = GA(sigma_cc, N);  

% Analysis based on BREV order (if necessary). Reconstruct the code.
% channels = bitrevorder(channels);
[~, channel_ordered] = sort(channels, 'descend');
info_bits = sort(channel_ordered(1 : CASCL_K), 'ascend');  
frozen_bits = ones(N , 1);                             
frozen_bits(info_bits) = 0;
info_bits_logical = logical(mod(frozen_bits + 1, 2));
nonfrozen_bits_logical = info_bits_logical;             % both info & CRC bits are all "non-frozen".


%% Start Simulation.
parfor i_Ebn0 = 1:N_Ebn0
    % estimate BLER at this Eb/n0.
    % model: Bin(n,p), To estimate p.
    Ebn0 = Ebn0_arr(i_Ebn0);
    sigma_sim = 1/(sqrt(2*R))*10^(-Ebn0/20);
    fprintf('Estimating BLER @ Eb/n0=%.2f dB\n', Ebn0);
    N_runs = 0;
    N_PCC_errs = 0;
    N_CRC_errs = 0;
    
    while ((N_PCC_errs < min_error) || (N_CRC_errs < min_error))
        random_bits = (rand([1,M])>0.5);
        PCC_encoded_bits = PCC_polar_encoder(random_bits, PCC);
        CRC_aided_bits = [random_bits, logical(mod(random_bits*crc_parity_check, 2))];  % row vector.
        
        u = zeros(1,N); % type: double array, row vector.
        u(nonfrozen_bits_logical) = CRC_aided_bits;
        x_CRC_BPSK = 1-2*my_polar_encoder(u,n);
        x_PCC_BPSK = 1-2*PCC_encoded_bits;
        
        noise_vec = sigma_sim * randn([1,N]);
        y_PCC = x_PCC_BPSK + noise_vec;   % add AWGN.
        llr_PCC = 2*y_PCC/(sigma_sim^2);

        y_CRC = x_CRC_BPSK + noise_vec;
        llr_CRC = 2*y_CRC/(sigma_sim^2);

        % PCC-SCL decoder.
        polar_info_esti_PCC = PCC_SCL_decoder(llr_PCC,L,PCC,lambda_offset,llr_layer_vec,bit_layer_vec);
        polar_info_esti_CRC = CASCL_decoder(llr_CRC,L,CASCL_K,~nonfrozen_bits_logical,H_crc,...
                                            lambda_offset, llr_layer_vec, bit_layer_vec);
        polar_info_esti_CRC = logical(polar_info_esti_CRC(1:M).');                          
        % Check correctness of PCC-SCL decoding result.
        % err_cnt_PCC = sum(xor(random_bits, polar_info_esti_PCC));
        % err_cnt_CRC = sum(xor(random_bits, polar_info_esti_CRC));
        if any(random_bits ~= polar_info_esti_PCC)
            N_PCC_errs = N_PCC_errs+1;
        end
        if any(random_bits ~= polar_info_esti_CRC)
            N_CRC_errs = N_CRC_errs + 1;
        end
        N_runs = N_runs + 1;
        if mod(N_runs, min_error/10) == 0
            fprintf('Estimating BLER @ Eb/n0=%.2f dB, Complete: %.2f%%\n', ...
            Ebn0, 100*(N_PCC_errs / min_error));
        end
    end
    BLER_CRC(i_Ebn0) = N_CRC_errs / N_runs;
    BLER_PCC(i_Ebn0) = N_PCC_errs / N_runs;
end

%% Plot BLER-Eb/n0 curves.
figure(1);
plot(Ebn0_arr, BLER_CRC, '-*');
hold on;
plot(Ebn0_arr, BLER_PCC, '-+');
set(gca, 'yscale', 'log');
legend('CRC-8', 'PCC-18, heuristic');
xlabel('E_b/n_0(dB)');
ylabel('BLER');
title('Comparison: CRC-Aided SCL and PCC-SCL');
grid on;

