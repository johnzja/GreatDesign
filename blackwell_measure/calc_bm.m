GF_info.alpha = 2;
GF_info.mult_table  = [0, 0, 0, 0; 0, 1, 2, 3; 0, 2, 3, 1; 0, 3, 1, 2];
GF_info.add_table   = [0, 1, 2, 3; 1, 0, 3, 2; 2, 3, 0, 1; 3, 2, 1, 0];

Wn = up_transform_sym(bm_dist, bm_dist, bin_centers, GF_info);
fprintf('sym up transform complete\n');
% Wn = up_transform_4D(bm_dist, bm_dist, bin_centers, GF_info);
% fprintf('original up transform complete\n');


Wp = down_transform_sym(bm_dist, bm_dist, bin_centers, GF_info);
fprintf('sym down transform complete\n');

%% Get the capacities.
IWn = get_I_4D(Wn, bin_centers);
IWp = get_I_4D(Wp, bin_centers);


%% Utils.

% All the utils have been written as .m files.