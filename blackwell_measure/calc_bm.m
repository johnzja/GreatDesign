GF_info.m = 2;          % GF(2^2).
GF_info.alpha = 2;
GF_info.mult_table  = [0, 0, 0, 0; 0, 1, 2, 3; 0, 2, 3, 1; 0, 3, 1, 2];
GF_info.add_table   = [0, 1, 2, 3; 1, 0, 3, 2; 2, 3, 0, 1; 3, 2, 1, 0];
GF_info.kernel_index_vec0 = kernel_index_vec0;          % boolean vector.
GF_info.kernel_index_vec1 = int32(kernel_index_vec1);   % int32 vector.
GF_info.kernel_index_mat0 = int32(kernel_index_mat0);   % int32 matrix.
GF_info.kernel_index_mat1 = int32(kernel_index_mat1);   % int32 matrix.

GF_info.num_threads = 6;
GF_info.jobs_each_thread = 4;

n = 7;
N = 2^n;
Nb = 2*N;       % Binary code length.

%% Polarization process.
addpath(strrep('D:\iTsinghua\Major\github_syncs\Encoding\PolarCpp\PolarCpp\x64\Release', '\', '/'));

Synth_channels = cell(n+1,1);
Synth_channels{1} = cell(1,1);
Synth_channels{1}{1} = bm_dist;

for layer = 1:n
    last_layer = Synth_channels{layer};
    
    this_layer_1 = cell(2^(layer-1), 1);
    this_layer_2 = cell(2^(layer-1), 1);
    
    this_layer_time_n = zeros(2^(layer-1),1);
    this_layer_time_p = zeros(2^(layer-1),1);
    
    for idx = 1:(2^(layer-1))
        bm_temp = last_layer{idx};
        tic;
        Wn = bm_polar_transform(bm_temp, bm_temp, bin_centers, GF_info, false);
        t = toc;
        this_layer_time_n(idx) = this_layer_time_n(idx) + t;
        tic;
        Wp = bm_polar_transform(bm_temp, bm_temp, bin_centers, GF_info, true);
        t = toc;
        this_layer_time_p(idx) = this_layer_time_p(idx) + t;
        
        this_layer_1{idx} = Wn;
        this_layer_2{idx} = Wp;
    end
    % Gather all the running times.
    mean_time_n = mean(this_layer_time_n);
    mean_time_p = mean(this_layer_time_p);
    total_time_cost = sum(this_layer_time_n) + sum(this_layer_time_p);
    Synth_channels{layer+1} = {this_layer_1{:}, this_layer_2{:}};
    fprintf('Layer %d polarization complete with mean Wn time cost = %.2fs and Wp time cost = %.2fs, TOTAL time = %.2fs\n', ...
        layer, mean_time_n, mean_time_p, total_time_cost);
end


%% Get the capacities.
% IWn = get_I_4D(Wn, bin_centers);
% IWp = get_I_4D(Wp, bin_centers);


%% Utils.

% All the utils have been written as .m files.