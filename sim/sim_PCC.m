function BLER = sim_PCC(PCC_struct, Ebn0, min_errors, L, force_rate)
    %% Gather information.
    K = PCC_struct.parity_bits_cnt;
    M = PCC_struct.info_bits_cnt;
    N = PCC_struct.N;
    
    n = log2(N);
    lambda_offset = 2.^(0:n);
    llr_layer_vec = get_llr_layer(N);
    bit_layer_vec = get_bit_layer(N);
    
    if ~exist('force_rate', 'var') || isempty(force_rate)
        R = M/N;    % code rate.
    else
        R = force_rate;
    end
    
    %% Start Simulation.
    % model: Binomial(n,p), To estimate p.
    
    sigma_sim = 1/(sqrt(2*R))*10^(-Ebn0/20);
    fprintf('Estimating BLER @ Eb/n0=%.2f dB for PCC-Polar SCL decoder.\n', Ebn0);
    
    N_parallel = 10;
    assert(mod(min_errors, N_parallel) == 0, 'invalid N_parallel!');
    min_errors_each = min_errors / N_parallel;
    N_runs_each = zeros(1, N_parallel);
    
    parfor p_iter = 1:N_parallel
        N_runs = 0;
        N_PCC_errs = 0;
        while N_PCC_errs < min_errors_each
            random_bits = (rand([1,M])>0.5);
            PCC_encoded_bits = PCC_polar_encoder(random_bits, PCC_struct);
            x_PCC_BPSK = 1-2*PCC_encoded_bits;

            noise_vec = sigma_sim * randn([1,N]);
            y_PCC = x_PCC_BPSK + noise_vec;   % add AWGN.
            llr_PCC = 2*y_PCC/(sigma_sim^2);

            % PCC-SCL decoder.
            polar_info_esti_PCC = PCC_SCL_decoder(llr_PCC,L,PCC_struct,...
                lambda_offset,llr_layer_vec,bit_layer_vec);

            % Check correctness of PCC-SCL decoding result.
            % err_cnt_PCC = sum(xor(random_bits, polar_info_esti_PCC));
            if any(random_bits ~= polar_info_esti_PCC)
                N_PCC_errs = N_PCC_errs+1;  % BLER.
            end
            
            N_runs = N_runs + 1;
            if mod(N_runs, min_errors) == 0
                fprintf('Job %d: Estimating BLER @ Eb/n0=%.2f dB, Complete: %.2f%%\n', ...
                p_iter, Ebn0, 100*(N_PCC_errs / min_errors));
            end
        end
        N_runs_each(p_iter) = N_runs;
    end
    
    BLER = min_errors / sum(N_runs_each);
end