function BLER = sim_CA_SCLF(K_CRC, N, M, Ebn0, min_errors, L, T)
    addpath('codes/');
    addpath('codes/polar/');
    addpath('codes/polar/GA/');
    addpath('sim/');
    
    %% Get CRC object.
    [~, ~, g] = get_crc_objective(K_CRC);
    n = log2(N);
    [G_crc, H_crc] = crc_generator_matrix(g, M);
    CRC_Qmatrix = G_crc(:, M+1 : end);        % Generator matrix for CRC bits.; G=[I,Q].
    
    %% Calculate some handy variables used in general SCL decoder.
    lambda_offset = 2.^(0:n);
    llr_layer_vec = get_llr_layer(N);
    bit_layer_vec = get_bit_layer(N);
    R = M/N;        % code rate.
    
    N_nonfrozen = M + K_CRC;
    %% Concatenation scheme: info bits => [info bits, CRC] => [info bits, CRC, PCC] => x_1^N, polar codeword.
    design_Ebn0 = 1.5;
    sigma_cc = 1/sqrt(2 * R) * 10^(-design_Ebn0/20);
    [channels, ~] = GA(sigma_cc, N);
    [~, channel_ordered] = sort(channels, 'descend');
    info_bits = sort(channel_ordered(1 : N_nonfrozen), 'ascend');
    frozen_bits = ones(N , 1);
    frozen_bits(info_bits) = 0;
    info_bits_logical = logical(mod(frozen_bits + 1, 2));
    
    CASCLF_config.M = M;                 % cnt. info bits.
    CASCLF_config.N = N;                 % cnt. code length.
    CASCLF_config.H_CRC = H_crc;         % CRC-Check matrix.
    CASCLF_config.Q_CRC = CRC_Qmatrix;   % Redundancy Generator matrix.
    CASCLF_config.K_CRC = K_CRC;         % cnt. additional CRC bits.
    CASCLF_config.nonfrozen_bits_logical = info_bits_logical;
    
    %% Construct SRCS set.
    nbl = info_bits_logical;
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
    
    %% Sort CS according to GA-estimated channel reliability.
    tmp = (1:N).';
    CS = tmp(CS);
    CS_reliability = channels(CS);
    [~, CS_order] = sort(CS_reliability, 'ascend');
    CS = CS(CS_order);

    %% Construct the information structure needed for the decoder.
    CASCLF_config.L = L;
    CASCLF_config.T = T;
    CASCLF_config.CS = CS;
    CASCLF_config.lambda_offset = lambda_offset;
    CASCLF_config.llr_layer_vec = llr_layer_vec;
    CASCLF_config.bit_layer_vec = bit_layer_vec;
    
    %% Start Simulation.
    % model: Binomial(n,p), To estimate p as BLER.
    sigma_sim = 1/(sqrt(2*R))*10^(-Ebn0/20);
    fprintf('Estimating BLER @ Eb/n0=%.2f dB for CRC-Polar SCL decoder.\n', Ebn0);
    
    % Add parallel support,
    N_parallel = 4;
    min_errors_each = min_errors / N_parallel;
    N_runs_each = zeros(1, N_parallel);
    
    parfor p_iter = 1:N_parallel
        N_CASCLF_errs = 0;
        N_runs = 0;
        while (N_CASCLF_errs < min_errors_each)
            random_bits = (rand([1,M])>0.5);
            CRC_aided_bits = [random_bits, logical(mod(random_bits*CASCLF_config.Q_CRC, 2))];  % row vector.
            u = zeros(1,N);
            u(info_bits_logical) = CRC_aided_bits;
            x = my_polar_encoder(u,n);

            x_bpsk = 1-2*x;
            noise_vec = sigma_sim * randn([1,N]);
            y = x_bpsk + noise_vec;   % add AWGN.

            llr = 2*y/(sigma_sim^2);

            % CRC-Polar SCLF decoding.
            polar_info_esti_CASCLF = CRC_SCLF_decoder(llr, CASCLF_config);

            % Check correctness of PCC-SCL decoding result.
            % err_cnt_PCC = sum(xor(random_bits, polar_info_esti_PCC));
            if any(random_bits ~= polar_info_esti_CASCLF)
                N_CASCLF_errs = N_CASCLF_errs+1;  % BLER.
            end

            N_runs = N_runs + 1;
            if mod(N_runs, min_errors/100) == 0
                fprintf('worker %d: Estimating BLER @ Eb/n0=%.2f dB, Complete: %.2f%%\n', ...
                p_iter, Ebn0, 100*(N_CASCLF_errs / min_errors_each));
            end
        end 
        N_runs_each(p_iter) = N_runs;
    end

    BLER = min_errors / sum(N_runs_each);
end

