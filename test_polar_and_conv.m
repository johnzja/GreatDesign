clear;
addpath('codes/');
addpath('codes/polar/');
addpath('codes/polar/GA/');
addpath('codes/conv');

setup_mapper;
setup_encoder;

%% Get CRC16 objects.
crc_length = 16;
[gen, det, g] = get_crc_objective(crc_length);

%% Setup simulation parameters.
n = 11;
N = 2^n;                % (2048, 1024+16) Polar code.
K = N/2 + crc_length;   % Final polar coding rate = 1/2.

N_sim = 5000;
N_info_bits = 1024;                                                                         % 1kb file to transmit.
Ebn0_arr = [8, 6, 5, 4, 3, 2.5, 2, 1.75, 1.5, 1.25, 1, 0.75, 0.5, 0, -0.5, -1, -1.5];       % Eb/n0 in dB.

Es = 1;                                     
Eb = Es * 2;                                % BPSK with R=1/2 rate.
n0_arr = Eb ./ (10.^(Ebn0_arr/10));         % n0 in linear scale.

N_n0s = length(n0_arr);


%% Construct Code using GA methods.
design_snr = 2.5;
R=1/2;
sigma_cc = 1/sqrt(2 * R) * 10^(-design_snr/20);
[channels, ~] = GA(sigma_cc, N);
[~, channel_ordered] = sort(channels, 'descend');
info_bits = sort(channel_ordered(1 : K), 'ascend');
frozen_bits = ones(N , 1);
frozen_bits(info_bits) = 0;
info_bits_logical = logical(mod(frozen_bits + 1, 2));

%% Setup some handy variables.
%codes parameters to avoid redundant calculations
lambda_offset = 2.^(0 : log2(N));
llr_layer_vec = get_llr_layer(N);
bit_layer_vec = get_bit_layer(N);


%Self-made CRC check
[G_crc, H_crc] = crc_generator_matrix(g, K - crc_length);
crc_parity_check = G_crc(:, K - crc_length + 1 : end)';

%% Co-sim with conv code (R=1/2), BPSK Modulation. Es = 1, Eb = 2.
err_bits_cnt_conv = zeros(N_n0s, 1);
err_bits_cnt_polar = zeros(N_n0s, 1);

parfor n0_iter = 1:N_n0s
    % calculate sigma using equation: n0 = 2*sigma^2.
    sigma = sqrt(n0_arr(n0_iter)/2);
    for sim_iter = 1:N_sim
        random_bits = (rand([1, N_info_bits])>0.5);
        info_with_crc = [random_bits.'; mod(crc_parity_check * (random_bits.'), 2)];
        u = zeros(N,1);
        u(info_bits_logical) = info_with_crc;
        
        conv_encoded_bits = conv_encode(random_bits, conv_encoder_conf);
        polar_encoded_bits = my_polar_encoder(u, log2(N));

        % simulation with BPSK channel.
        bpsk_conv = 1 - 2 * conv_encoded_bits;
        bpsk_polar = 1 - 2 * polar_encoded_bits.';

        conv_noise = randn(1, length(bpsk_conv));
        polar_noise = randn(1, length(bpsk_polar));
        y_conv = bpsk_conv + sigma * conv_noise;
        y_polar = bpsk_polar + sigma * polar_noise;

        llr_polar = 2/sigma^2*y_polar.';
        % process y_conv into complex symbols with energy 1.
        complex_y_conv = (y_conv(1:2:end) + 1j*y_conv(2:2:end))/sqrt(2); % complex syms with energy 1. Eb = Es = 1 here.
        ch = ones(1, length(complex_y_conv));
        pred_L2 = bit_demapping(complex_y_conv, length(y_conv), mapping_conf, ch, ch_conf, 1);

        % decoders.
        conv_decoded_bits_soft = fast_conv_decode(pred_L2, conv_encoder_conf, true); % Perform soft-decode.
        polar_decoded_bits = CASCL_decoder(llr_polar, 16, K, frozen_bits, H_crc, lambda_offset, llr_layer_vec, bit_layer_vec);
        polar_decoded_bits = polar_decoded_bits(1:N_info_bits).';   % Restore to row vectors.
        % BER comparison.
        err_bits_cnt_conv(n0_iter) = err_bits_cnt_conv(n0_iter) + sum(xor(conv_decoded_bits_soft, random_bits));
        err_bits_cnt_polar(n0_iter) = err_bits_cnt_polar(n0_iter) + sum(xor(polar_decoded_bits, random_bits));
        
       % Display running info.
        if mod(sim_iter, floor(N_sim/10))==0
            disp(['SOFT Eb/n0=', num2str(Ebn0_arr(n0_iter)),'dB: ',num2str(sim_iter/N_sim*100),'% complete']);
        end
    end
end

%% Display.
figure;
hold on;
plot(Ebn0_arr, (err_bits_cnt_conv/(N_info_bits*N_sim)).', '-x');
plot(Ebn0_arr, (err_bits_cnt_polar/(N_info_bits*N_sim)).', '-*');

legend('conv', 'polar');
set(gca, 'yscale', 'log');
title('BER-Eb/n_0 Curve');
xlabel('Eb/n0(dB)');
ylabel('BER');
grid on;
