study_BM_QPSK;
% addpath(strrep('D:\iTsinghua\Major\github_syncs\Encoding\PolarCpp\PolarCpp\x64\Debug', '\', '/'));
addpath(strrep('D:\iTsinghua\Major\github_syncs\Encoding\PolarCpp\PolarCpp\x64\Release', '\', '/'));

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

load('data/bm_qpsk_64.mat', 'Synth_channels');
Wn_r = Synth_channels{2}{2};
% tic;
% Wn_r = down_transform_sym(bm_dist, bm_dist, bin_centers, GF_info);
% toc;

%% TEST.
d = sum(sum(sum(abs(Wn-Wn_r))))




