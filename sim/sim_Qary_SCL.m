function BLER = sim_Qary_SCL(N, M, partially_frozen, info_syms, Ebn0, min_errors, L, GF_info)
% addpath(strrep('D:\iTsinghua\Major\github_syncs\Encoding\PolarCpp\PolarCpp\x64\Release', '\', '/'));
addpath('../PolarCpp/PolarCpp/Release/');
addpath('codes/polar/');

%% Setup simulation parameters.
R = M/N;
n = log2(N);

Esn0 = Ebn0 + 10*log10(GF_info.m*R);
sigma = (10^(-Esn0/20));
fprintf('Estimating BLER @ Eb/n0=%.2f dB for Q-ary SCL decoder.\n', Ebn0);

bpsk_mat = [1,1; -1,1; 1,-1; -1,-1]; % (I,Q).
m = GF_info.m;
pwrs = 2.^(0:m-1);

if ~partially_frozen
    frozen_bits = ~info_syms;
else
    error('partially-frozen SCL is not supported yet!');
end

% Setup mex-decoder.
decoder_config.partially_frozen = false;
decoder_config.is_qary          = true;
decoder_config.is_LLR           = true;
decoder_config.is_Genie         = false;
decoder_config.update_decoder   = true;       % it can be automatically modified into false.
decoder_config.is_list          = true;
decoder_config.L                = L;

% Setup parallel parameters.
N_parallel = 25;
assert(mod(min_errors, N_parallel) == 0, 'invalid N_parallel!');
min_errors_each = min_errors / N_parallel;
N_runs_each = zeros(1, N_parallel);

%% Simulation Loop
parfor p_iter = 1:N_parallel
    N_runs = 0;
    N_SCL_errs = 0;
    while N_SCL_errs < min_errors_each
        random_syms = randi([0, 2^m-1], [1, M]);
        u = zeros(1, N);
        u(info_syms) = random_syms;
        x = GF_polar_encoder(u, n, GF_info);

        x_bpsk = bpsk_mat(x+1,:).';
        x_bpsk = x_bpsk(:).';

        % add AWGN noise.
        noise_vec = randn([1, N*m]);
        y = x_bpsk + noise_vec*sigma;
        llr = 2*y/sigma^2;

        % call mex decoder.
        polar_info_esti = Qary_SC_Decoder(llr, N, m, frozen_bits, GF_info.alpha, decoder_config);

        % convert into symbols.
        info_esti = zeros(1, M);
        for idx = 1:M
            sym_bits = polar_info_esti(1+(idx-1)*m:2+(idx-1)*m);
            info_esti(idx) = sym_bits * pwrs.';
        end

        e = any(info_esti ~= random_syms);
            
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

