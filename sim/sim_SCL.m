function BLER = sim_SCL(N, M, info_bits_logical, Ebn0, min_errors, L)
% addpath(strrep('D:\iTsinghua\Major\github_syncs\Encoding\PolarCpp\PolarCpp\x64\Release', '\', '/'));
addpath('../PolarCpp/PolarCpp/Release/');
addpath('codes/polar/');

%% Setup simulation parameters.
R = M/N;
n = log2(N);
sigma_sim = 1/sqrt(2*R)*(10^(-Ebn0/20));
frozen_bits = ~info_bits_logical;

fprintf('Estimating BLER @ Eb/n0=%.2f dB for Binary SCL decoder.\n', Ebn0);

% Setup mex-decoder.
decoder_config.partially_frozen = false;
decoder_config.is_qary          = false;
decoder_config.is_LLR           = true;
decoder_config.is_Genie         = false;
decoder_config.update_decoder   = true;       % it can be automatically modified into false.
decoder_config.is_list          = true;
decoder_config.L                = L;

% Setup parallel parameters.
N_parallel = 6;
assert(mod(min_errors, N_parallel) == 0, 'invalid N_parallel!');
min_errors_each = min_errors / N_parallel;
N_runs_each = zeros(1, N_parallel);

%% Simulation Loop
parfor p_iter = 1:N_parallel
    N_runs = 0;
    N_SCL_errs = 0;
    while N_SCL_errs < min_errors_each
        random_bits = rand([1,M])>0.5;
        u = zeros(1, N);
        u(info_bits_logical) = random_bits;
        x = my_polar_encoder(u,n);

        % add AWGN noise.
        x_bpsk = 1-2*x;
        noise_vec = sigma_sim * randn([1,N]);
        y = x_bpsk + noise_vec;
        llr = 2*y/sigma_sim^2;

        % call mex decoder.
        polar_info_esti = Qary_SC_Decoder(llr, N, 1, frozen_bits, 0, decoder_config);
        e = any(polar_info_esti ~= random_bits);
            
        if e
            N_SCL_errs = N_SCL_errs+1;
        end
            
        N_runs = N_runs + 1;
        if mod(N_runs, min_errors) == 0
            fprintf('Job %d: Estimating BLER @ Eb/n0=%.2f dB, Complete: %.2f%%\n', ...
            p_iter, Ebn0, 100*(N_SCL_errs / min_errors_each));
        end
    end
    N_runs_each(p_iter) = N_runs;
end

BLER = min_errors / sum(N_runs_each);

end

