%% Study 3-D Blackwell Measure.
clc;clear;
addpath('codes/polar/GA/');

P = 1;                  % The constellation power is 2P.
A = sqrt(P);            % Amplitude on I or Q.

A0 = A*[1,1];
A1 = A*[-1,1];
A2 = A*[-1,-1];
A3 = A*[1, -1];

Esn0 = 1.5;             % unit: dB

% sigma = 0.8;
sigma = A * (10^(-Esn0/20));


%% Do QPSK sampling.
N_sample = 10000;
bm_posterior_coord = zeros(N_sample, 3);

for idx = 1:N_sample
    b0 = rand() > 0.5;
    b1 = rand() > 0.5;
    
    noise = sigma*randn([1,2]);
    y = (1-2*[b0, b1]) + noise;
    r0 = norm(y-A0);
    r1 = norm(y-A1);
    r2 = norm(y-A2);
    r3 = norm(y-A3);
    
    s0 = exp(-r0^2/(2*sigma^2));
    s1 = exp(-r1^2/(2*sigma^2));
    s2 = exp(-r2^2/(2*sigma^2));
    s3 = exp(-r3^2/(2*sigma^2));
    
    s = sum([s0, s1, s2, s3]);
    s0 = s0 / s;
    s1 = s1 / s;
    s2 = s2 / s;
    s3 = s3 / s;
    
    bm_posterior_coord(idx, :) = [s0, s1, s2];
end

scatter3(bm_posterior_coord(:,1), bm_posterior_coord(:,2), bm_posterior_coord(:,3), '*');
xlabel('s0'); ylabel('s1'); zlabel('s2');
axis equal;

% Tetrahedron Representation.
% symbols: [A0, A1, A2, A3].
% bit mapping: [b0b1] = [00, 10, 01, 11].
% Posterioiri proba: [s0, s1, s2, s3].
% [s0, s1, s2] in use.
v0 = [1, 0, 0];
v1 = [0.5, sqrt(3)/2, 0];
v2 = [0.5, 0.5/sqrt(3), sqrt(6)/3];

bm_display_coord = zeros(size(bm_posterior_coord));
for idx = 1:N_sample
    s0 = bm_posterior_coord(idx, 1);
    s1 = bm_posterior_coord(idx, 2);
    s2 = bm_posterior_coord(idx, 3);
    bm_display_coord(idx, :) = s0 * v0 + s1 * v1 + s2 * v2;
end

figure(2);
scatter3(bm_display_coord(:,1), bm_display_coord(:,2), bm_display_coord(:,3), '*');
xlabel('x'); ylabel('y'); zlabel('z');
axis equal;

%% Calculate the exact Blackwell Measure for QPSK-AWGN channel.
N_bins_each_dim = 32;
bin_centers = linspace(0, 1, N_bins_each_dim);
% index of the bins: 1 ~ N_bins_each_dim.
% index0 + index1 <= N + 1.
bm_dist = zeros(N_bins_each_dim, N_bins_each_dim, N_bins_each_dim);

% Enumerate all the bins: use the following code:
% for idx0 = 1:N_bins_each_dim
%     for idx1 = 1:(N_bins_each_dim+1-idx0)
%         for idx2 = 1:(N_bins_each_dim+2-idx0-idx1)
%             s0 = bin_centers(idx0);
%             s1 = bin_centers(idx1);
%             s2 = bin_centers(idx2);
%             s3 = 1-s0-s1-s2;
%             
%         end
%     end
% end

Nb = 2000;
p = 1/Nb;

proba_bins = linspace(0,1,Nb+1);
confidence_bounds = zeros(2, Nb);
expected_value = zeros(1, Nb);
confidence_bounds(1, 1) = -Inf;
confidence_bounds(2, Nb) = +Inf;

for idx = 1:Nb-1
    t = norminv(proba_bins(idx+1));
    confidence_bounds(1, idx+1) = t;
    confidence_bounds(2, idx) = t;
end

for idx = 1:Nb
    L = confidence_bounds(1, idx);
    U = confidence_bounds(2, idx);
    if L == -Inf
        E = -1/sqrt(2*pi)*exp(-0.5*U^2);
    elseif U == +Inf
        E = 1/sqrt(2*pi)*exp(-0.5*L^2);
    else
        E = -1/sqrt(2*pi)*(exp(-0.5*U^2)-exp(-0.5*L^2));
    end
    expected_value(idx) = E/p;
end

% Calculate the QPSK BM.
for idx_i = 1:Nb
    for idx_q = 1:Nb
        Ei = A + sigma * expected_value(idx_i);
        Eq = A + sigma * expected_value(idx_q);

        y = [Ei, Eq];                   % Received signal.
        r0 = norm(y-A0);
        r1 = norm(y-A1);
        r2 = norm(y-A2);
        r3 = norm(y-A3);

        s0 = exp(-r0^2/(2*sigma^2));
        s1 = exp(-r1^2/(2*sigma^2));
        s2 = exp(-r2^2/(2*sigma^2));
        s3 = exp(-r3^2/(2*sigma^2));

        s = sum([s0, s1, s2, s3]);
        s0 = s0 / s;
        s1 = s1 / s;
        s2 = s2 / s;
        s3 = s3 / s;
       
        [idx0, idx1, idx2, idx3] = convert_dist_into_index([s0, s1, s2, s3], N_bins_each_dim);
        
        bm_dist(idx0, idx1, idx2) = bm_dist(idx0, idx1, idx2) + (p^2)/4;
        bm_dist(idx1, idx2, idx3) = bm_dist(idx1, idx2, idx3) + (p^2)/4;
        bm_dist(idx2, idx3, idx0) = bm_dist(idx2, idx3, idx0) + (p^2)/4;
        bm_dist(idx3, idx0, idx1) = bm_dist(idx3, idx0, idx1) + (p^2)/4;
    end
end

%% Display bm_dist using 3D density plot...

I_theoretical = log2(1+(A/sigma)^2);        % Maximum value achievable.
fprintf('Theoretical AWGNC Gaussian-input I = %f bits/ch.use\n\n', I_theoretical);

I_QPSK = 2*Capacity_Binary_AWGN(A, sigma);
fprintf('Theoretical QPSK channel I = %f bits/ch.use\n\n', I_QPSK);

I = get_I_4D(bm_dist, bin_centers);
fprintf('BM I = %f bits\n\n', I);


%% Plot!
plot_BM(bm_dist, bin_centers);




