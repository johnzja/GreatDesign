
% addpath(strrep('D:\iTsinghua\Major\github_syncs\Encoding\PolarCpp\PolarCpp\x64\Debug', '\', '/'));
addpath(strrep('D:\iTsinghua\Major\github_syncs\Encoding\PolarCpp\PolarCpp\x64\Debug', '\', '/'));

GF_info.alpha = 2;
GF_info.mult_table  = [0, 0, 0, 0; 0, 1, 2, 3; 0, 2, 3, 1; 0, 3, 1, 2];
GF_info.add_table   = [0, 1, 2, 3; 1, 0, 3, 2; 2, 3, 0, 1; 3, 2, 1, 0];

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

% tic;
% Wn_r = up_transform_4D(bm_dist, bm_dist, bin_centers, GF_info);
% toc;






