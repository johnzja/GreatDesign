function BLER = sim_CASCL_given_construction(K_CRC, N, M, info_bits, Ebn0, min_errors, L)
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
    info_bits_logical = info_bits;
    
    %% Start Simulation.
    % model: Binomial(n,p), To estimate p as BLER.
    sigma_sim = 1/(sqrt(2*R))*10^(-Ebn0/20);
    fprintf('Estimating BLER @ Eb/n0=%.2f dB for CA-SCL decoder.\n', Ebn0);
    
    % Add parallel support,
    N_parallel = 2;
    min_errors_each = min_errors / N_parallel;
    N_runs_each = zeros(1, N_parallel);
    
    parfor p_iter = 1:N_parallel
        N_CASCL_errs = 0;
        N_runs = 0;
        while (N_CASCL_errs < min_errors_each)
            random_bits = (rand([1,M])>0.5);
            CRC_aided_bits = [random_bits, logical(mod(random_bits*CRC_Qmatrix, 2))];  % row vector.
            u = zeros(1,N);
            u(info_bits_logical) = CRC_aided_bits;
            x = my_polar_encoder(u,n);

            x_bpsk = 1-2*x;
            noise_vec = sigma_sim * randn([1,N]);
            y = x_bpsk + noise_vec;   % add AWGN noise.

            llr = 2*y/(sigma_sim^2);

            % CRC-Polar SCLF decoding.
            polar_info_esti_CASCL = CASCL_decoder(llr, L, M+K_CRC, ~info_bits_logical, H_crc, lambda_offset, llr_layer_vec, bit_layer_vec);
            polar_info_esti_CASCL = polar_info_esti_CASCL(1:M).';
            % Check correctness of PCC-SCL decoding result.
            if any(random_bits ~= polar_info_esti_CASCL)
                N_CASCL_errs = N_CASCL_errs+1;    % BLER.
            end

            N_runs = N_runs + 1;
            if mod(N_runs, min_errors) == 0
                fprintf('worker %d: Estimating BLER @ Eb/n0=%.2f dB, Complete: %.2f%%\n', ...
                p_iter, Ebn0, 100*(N_CASCL_errs / min_errors_each));
            end
        end 
        N_runs_each(p_iter) = N_runs;
    end

    BLER = min_errors / sum(N_runs_each);
end

