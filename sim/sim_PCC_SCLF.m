function varargout = sim_PCC_SCLF(K_CRC, K_PCC, N, M, Ebn0, min_errors, L, T)
    %% Setup path references.
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
    
    %% Concatenation scheme: info bits => [info bits, CRC] => [info bits, CRC, PCC] => x_1^N, polar codeword.
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
    
    %% Sort CS according to GA-estimated channel reliability.
    tmp = (1:N).';
    CS = tmp(CS);
    CS_reliability = GA_ch_reliability(CS);
    [~, CS_order] = sort(CS_reliability, 'ascend');
    CS = CS(CS_order);

    %% Construct the information structure needed for the PCC-CRC-Polar decoder.
    PCC_CRC_polar_decoder_info.L = L;
    PCC_CRC_polar_decoder_info.lambda_offset = lambda_offset;
    PCC_CRC_polar_decoder_info.llr_layer_vec = llr_layer_vec;
    PCC_CRC_polar_decoder_info.bit_layer_vec = bit_layer_vec;
    PCC_CRC_polar_decoder_info.CS = CS;
    PCC_CRC_polar_decoder_info.T = T;
    
    %% Start Simulation.
    % model: Binomial(n,p), To estimate p as BLER.
    sigma_sim = 1/(sqrt(2*R))*10^(-Ebn0/20);
    fprintf('Estimating BLER @ Eb/n0=%.2f dB for CRC-PCC-Polar SCLF decoder.\n', Ebn0);
    
    % Add parallel support,
    N_parallel = 6;            
    assert(mod(min_errors, N_parallel) == 0, 'invalid N_parallel!');
    min_errors_each = min_errors / N_parallel;
    N_runs_each = zeros(1, N_parallel);
    N_trials_each = zeros(1, N_parallel);
    
    parfor p_iter = 1:N_parallel
        N_PCC_SCLF_errs = 0;
        N_runs = 0;
        N_trials = 0;
        decoder_info = PCC_CRC_polar_decoder_info;

        while N_PCC_SCLF_errs < min_errors_each
            random_bits = (rand([1,M])>0.5);
            CRC_aided_bits = [random_bits, logical(mod(random_bits*PCC_CRC_polar_config.Q_CRC, 2))];  % row vector.
            PCC_CRC_encoded_bits = PCC_polar_encoder(CRC_aided_bits, PCC_CRC_polar_config.PCC_conf);

            x_PCC_BPSK = 1-2*PCC_CRC_encoded_bits;
            noise_vec = sigma_sim * randn([1,N]);
            y_PCC_CRC = x_PCC_BPSK + noise_vec;   % add AWGN.

            llr_PCC_CRC = 2*y_PCC_CRC/(sigma_sim^2);

            % CRC-PCC-Polar SCLF Decoding.
            [polar_info_esti_PCC_CRC, decoder_info] = PCC_CRC_SCLF_decoder(...
                llr_PCC_CRC, PCC_CRC_polar_config, decoder_info);
            
            % Check correctness of PCC-SCL decoding result.
            if any(random_bits ~= polar_info_esti_PCC_CRC)
                N_PCC_SCLF_errs = N_PCC_SCLF_errs+1;  % BLER.
            end

            N_runs = N_runs + 1;
            N_trials = N_trials + 1 + decoder_info.num_flip_trials;
            
            if mod(N_runs, min_errors) == 0
                fprintf('Worker %d: Estimating BLER @ Eb/n0=%.2f dB, Complete: %.2f%%\n', ...
                p_iter, Ebn0, 100*(N_PCC_SCLF_errs / min_errors_each));
            end
        end
        
        N_runs_each(p_iter) = N_runs;
        N_trials_each(p_iter) = N_trials;
    end
    
    BLER = min_errors / sum(N_runs_each);
    if nargout == 1
        varargout = {BLER};
    elseif nargout == 2
        varargout = {BLER, sum(N_trials_each)/sum(N_runs_each)};
    else
        error('%s: invalid num of output arguments.', mfilename);
    end
end