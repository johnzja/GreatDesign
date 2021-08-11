%% Study 3-D Blackwell Measure.
clc;clear;
addpath('codes/polar/GA/');
addpath('blackwell_measure/');

P = 1;                  % The constellation power is 2P.
A = sqrt(P);            % Amplitude on I or Q.

A0 = A*[ 1,  1];
A1 = A*[-1,  1];
A2 = A*[ 1, -1];
A3 = A*[-1, -1];

Esn0 = 0.187;             % unit: dB

sigma = A * (10^(-Esn0/20));
I_theoretical = log2(1+(A/sigma)^2);            % Maximum value achievable.
fprintf('Theoretical AWGNC Gaussian-input I = %f bits/ch.use\n', I_theoretical);

I_QPSK = 2*Capacity_Binary_AWGN(A, sigma);
fprintf('Theoretical QPSK channel I = %f bits/ch.use\n', I_QPSK);

%% Do QPSK sampling.
% N_sample = 10000;
% bm_posterior_coord = zeros(N_sample, 3);
% 
% for idx = 1:N_sample
%     b0 = rand() > 0.5;
%     b1 = rand() > 0.5;
%     
%     noise = sigma*randn([1,2]);
%     y = (1-2*[b0, b1]) + noise;
%     r0 = norm(y-A0);
%     r1 = norm(y-A1);
%     r2 = norm(y-A2);
%     r3 = norm(y-A3);
%     
%     s0 = exp(-r0^2/(2*sigma^2));
%     s1 = exp(-r1^2/(2*sigma^2));
%     s2 = exp(-r2^2/(2*sigma^2));
%     s3 = exp(-r3^2/(2*sigma^2));
%     
%     s = sum([s0, s1, s2, s3]);
%     s0 = s0 / s;
%     s1 = s1 / s;
%     s2 = s2 / s;
%     s3 = s3 / s;
%     
%     bm_posterior_coord(idx, :) = [s0, s1, s2];
% end
% 
% scatter3(bm_posterior_coord(:,1), bm_posterior_coord(:,2), bm_posterior_coord(:,3), '*');
% xlabel('s0'); ylabel('s1'); zlabel('s2');
% axis equal;
% 
% % Tetrahedron Representation.
% % symbols: [A0, A1, A2, A3].
% % bit mapping: [b0b1] = [00, 10, 01, 11].
% % Posterioiri proba: [s0, s1, s2, s3].
% % [s0, s1, s2] in use.
% v0 = [1, 0, 0];
% v1 = [0.5, sqrt(3)/2, 0];
% v2 = [0.5, 0.5/sqrt(3), sqrt(6)/3];
% 
% bm_display_coord = zeros(size(bm_posterior_coord));
% for idx = 1:N_sample
%     s0 = bm_posterior_coord(idx, 1);
%     s1 = bm_posterior_coord(idx, 2);
%     s2 = bm_posterior_coord(idx, 3);
%     bm_display_coord(idx, :) = s0 * v0 + s1 * v1 + s2 * v2;
% end
% 
% figure(2);
% scatter3(bm_display_coord(:,1), bm_display_coord(:,2), bm_display_coord(:,3), '*');
% xlabel('x'); ylabel('y'); zlabel('z');
% axis equal;

%% Calculate the exact Blackwell Measure for QPSK-AWGN channel.
N_bins_each_dim = 150;
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

Nb = 2500;
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
        bm_dist(idx1, idx0, idx3) = bm_dist(idx1, idx0, idx3) + (p^2)/4;
        bm_dist(idx2, idx3, idx0) = bm_dist(idx2, idx3, idx0) + (p^2)/4;
        bm_dist(idx3, idx2, idx1) = bm_dist(idx3, idx2, idx1) + (p^2)/4;
    end
end
fprintf('Construct QPSK BM complete!\n');

%% Display bm_dist using 3D density plot...


I = get_I_4D(bm_dist, bin_centers);
fprintf('BM I = %f bits\n\n', I);

%% Decomposition: Decompose into 4-ary symmetric channels and evaluate the capacity.
bm_dist_visited = false(N_bins_each_dim,N_bins_each_dim,N_bins_each_dim);
I = 0;
total_mass = 0;

for idx0 = 1:N_bins_each_dim
    for idx1 = 1:(N_bins_each_dim+1-idx0)
        for idx2 = 1:(N_bins_each_dim+2-idx0-idx1)
            idx3 = N_bins_each_dim+3-idx0-idx1-idx2;
            
            b0 = bm_dist_visited(idx0, idx1, idx2);
            b1 = bm_dist_visited(idx1, idx0, idx3);
            b2 = bm_dist_visited(idx2, idx3, idx0);
            b3 = bm_dist_visited(idx3, idx2, idx1);
            if any([b0,b1,b2,b3])
                assert(all([b0,b1,b2,b3]));
                continue;
            end
            
            t0 = bm_dist(idx0, idx1, idx2);
            t1 = bm_dist(idx1, idx0, idx3);
            t2 = bm_dist(idx2, idx3, idx0);
            t3 = bm_dist(idx3, idx2, idx1);
            
            if ~(t0 == t1 && t1 == t2 && t2 == t3)
                warning('Symmetry violated');
            end
            
            bm_dist_visited(idx0, idx1, idx2)=true;
            bm_dist_visited(idx1, idx0, idx3)=true;
            bm_dist_visited(idx2, idx3, idx0)=true;
            bm_dist_visited(idx3, idx2, idx1)=true;
           
            proba = t0+t1+t2+t3;
            if proba == 0
                continue;
            end
            
            Iq = 2;
            p = [bin_centers(idx0), bin_centers(idx1),bin_centers(idx2),bin_centers(idx3)].';
            for k = 1:4
                if p(k) > 0
                    Iq = Iq + p(k)*log2(p(k));
                end
            end
            I = I + Iq * proba;
            total_mass = total_mass + proba;
        end
    end
end

fprintf('Sym BM I = %f\n', I);


%% Plot!
plot_BM(bm_dist, bin_centers);

%% Study the symmetric properties and the symmetric kernel.
bm_dist_visited = false(N_bins_each_dim,N_bins_each_dim,N_bins_each_dim);
bm_dist_kernel = zeros(N_bins_each_dim,N_bins_each_dim,N_bins_each_dim);

N_kernel_cnt = 0;
for idx0 = 1:N_bins_each_dim
    for idx1 = 1:(N_bins_each_dim+1-idx0)
        for idx2 = 1:(N_bins_each_dim+2-idx0-idx1)
            idx3 = N_bins_each_dim+3-idx0-idx1-idx2;
            
            b0 = bm_dist_visited(idx0, idx1, idx2);
            b1 = bm_dist_visited(idx1, idx0, idx3);
            b2 = bm_dist_visited(idx2, idx3, idx0);
            b3 = bm_dist_visited(idx3, idx2, idx1);
            if any([b0,b1,b2,b3])
                assert(all([b0,b1,b2,b3]));
                continue;
            end
            
            bm_dist_visited(idx0, idx1, idx2)=true;
            bm_dist_visited(idx1, idx0, idx3)=true;
            bm_dist_visited(idx2, idx3, idx0)=true;
            bm_dist_visited(idx3, idx2, idx1)=true;
            
            N_kernel_cnt = N_kernel_cnt + 1;
            [~, idx_max] = max([idx0, idx1, idx2, idx3]);
            if idx_max == 1
                bm_dist_kernel(idx3, idx2, idx1) = 1.0;
            elseif idx_max == 2
                bm_dist_kernel(idx2, idx3, idx0) = 1.0;
            elseif idx_max == 3
                bm_dist_kernel(idx1, idx0, idx3) = 1.0;
            else
                bm_dist_kernel(idx0, idx1, idx2) = 1.0;
            end
                
        end
    end
end

plot_BM(bm_dist_kernel, bin_centers);
N_bins_3D = N_bins_each_dim*(N_bins_each_dim+1)*(N_bins_each_dim+2)/6;
assert(N_bins_3D == N_kernel_cnt*4);

% Represent the kernel with index.
kernel_index_vec0 = false(N_bins_each_dim,1);   % axis: idx0
kernel_index_vec1 = zeros(N_bins_each_dim,2);   % axis: idx1, lower & upper bound.
kernel_index_mat0 = zeros(N_bins_each_dim, N_bins_each_dim);    % axis: idx2, lower bound.
kernel_index_mat1 = zeros(N_bins_each_dim, N_bins_each_dim);    % axis: idx2, upper bound.

for idx0 = 1:N_bins_each_dim
    flag_v = 0;
    flag_exist_1 = 0;
    for idx1 = 1:(N_bins_each_dim+1-idx0)
        flag_u = 0;
        flag_exist_2 = 0;
        for idx2 = 1:(N_bins_each_dim+2-idx0-idx1)
            k = bm_dist_kernel(idx0, idx1, idx2);
            if flag_u == 0 && k > 0
                kernel_index_mat0(idx0, idx1) = idx2;   % setup lower bound.
                flag_u = 1;
                flag_exist_2 = 1;
            elseif flag_u == 1 && k == 0
                kernel_index_mat1(idx0, idx1) = idx2-1; % setup upper bound.
                flag_u = 0;
                break;
            end
        end
        
        if flag_u == 1
            kernel_index_mat1(idx0, idx1) = N_bins_each_dim+2-idx0-idx1;    % max. possible upper bound.
        end
        
        if flag_v == 0 && flag_exist_2 == 1
            kernel_index_vec1(idx0, 1) = idx1;      % setup lower bound;
            flag_v = 1;
            flag_exist_1 = 1;
        elseif flag_v == 1 && flag_exist_2 == 0
            kernel_index_vec1(idx0, 2) = idx1 - 1;
            flag_v = 0;
            break;
        end
    end
    if flag_v == 1
        kernel_index_vec1(idx0, 2) = N_bins_each_dim+1-idx0;
    end
    if flag_exist_1 == 1
        kernel_index_vec0(idx0) = true;
    end
end

% Count the total bins included in kernel.
N_kernel_cnt_new = 0;
for idx0 = 1:N_bins_each_dim
    if kernel_index_vec0(idx0)
        lb = kernel_index_vec1(idx0, 1);
        ub = kernel_index_vec1(idx0, 2);
        for idx1 = lb:ub
            lb2 = kernel_index_mat0(idx0, idx1);
            ub2 = kernel_index_mat1(idx0, idx1);
            N_kernel_cnt_new = N_kernel_cnt_new + (ub2-lb2+1);
        end
    end
end

assert(N_kernel_cnt == N_kernel_cnt_new);
fprintf('Kernel construction complete!\n');


