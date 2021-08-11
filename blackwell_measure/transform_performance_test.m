% study_BM_QPSK;
% addpath(strrep('D:\iTsinghua\Major\github_syncs\Encoding\PolarCpp\PolarCpp\x64\Debug', '\', '/'));
addpath(strrep('..\PolarCpp\PolarCpp\x64\Release', '\', '/'));

GF_info.alpha = 2;
GF_info.mult_table  = [0, 0, 0, 0; 0, 1, 2, 3; 0, 2, 3, 1; 0, 3, 1, 2];
GF_info.add_table   = [0, 1, 2, 3; 1, 0, 3, 2; 2, 3, 0, 1; 3, 2, 1, 0];
GF_info.kernel_index_vec0 = kernel_index_vec0;          % boolean vector.
GF_info.kernel_index_vec1 = int32(kernel_index_vec1);   % int32 vector.
GF_info.kernel_index_mat0 = int32(kernel_index_mat0);   % int32 matrix.
GF_info.kernel_index_mat1 = int32(kernel_index_mat1);   % int32 matrix.

%% Pure MATLAB.


%% CPP-MEX.

% N_bins = 4;
% bin_centers = linspace(0, 1, N_bins);
% bm_dist = zeros(N_bins, N_bins, N_bins);
% for idx_1 = 1:N_bins
%     for idx_2 = 1:N_bins
%         for idx_3 = 1:N_bins
%             bm_dist(idx_1, idx_2, idx_3) = (idx_1-1) + N_bins*(idx_2-1) + N_bins^2*(idx_3-1);
%         end
%     end
% end
tic;
Wn = bm_polar_transform(bm_dist, bm_dist, bin_centers, GF_info, true);
toc;

% load('data/bm_qpsk_64.mat', 'Synth_channels');
% Wn_r = Synth_channels{2}{2};
tic;
Wn_r = down_transform_sym(bm_dist, bm_dist, bin_centers, GF_info);
toc;

%% TEST.
d = sum(sum(sum(abs(Wn-Wn_r))))

%% Test whether the bm is right!!

N_bins = N_bins_each_dim;
i0 = 3; i1 = 2; i2 = 2; 
i3 = N_bins + 3 - i0 - i1 - i2;

j0 = 2; j1 = 2; j2 = 1; 
j3 = N_bins + 3 - j0 - j1 - j2;

bm_1 = zeros(N_bins,N_bins,N_bins);
bm_2 = zeros(N_bins,N_bins,N_bins);

bm_1(i0, i1, i2) = 1/4;
bm_1(i1, i0, i3) = 1/4;
bm_1(i2, i3, i0) = 1/4;
bm_1(i3, i2, i1) = 1/4;

bm_2(j0, j1, j2) = 1/4;
bm_2(j1, j0, j3) = 1/4;
bm_2(j2, j3, j0) = 1/4;
bm_2(j3, j2, j1) = 1/4;

% use mex-bm.
wp = bm_polar_transform(bm_1, bm_2, bin_centers, GF_info, true);
wp2 = down_transform_sym(bm_1, bm_2, bin_centers, GF_info);

d = sum(sum(sum(abs(wp-wp2))))

% Find the ~=0 elements.
z1 = find(wp ~= 0);
for idx = 1:length(z1)
    z = z1(idx)-1;
    k2 = floor(z/(N_bins*N_bins)) + 1;
    z = z - (k2-1)*N_bins*N_bins;
    k1 = floor(z/N_bins)+1;
    z = z - (k1-1)*N_bins;
    k0 = z+1;
    fprintf('(%d, %d, %d), %f\n', k0, k1, k2, wp(k0,k1,k2));
end

% conclusion: Quantization problem HERE.