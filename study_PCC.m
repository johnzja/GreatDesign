clear;clc;
addpath('sim/');

%% Simulate PCC code.
% Compare: K_CRC=4, K_PCC=8 code+SCLF v.s. K_PCC=12 code.
N = 512;
M = 256;
K_PCC = 12;

design_Ebn0 = 1.5;
[PCC_structs, ~] = get_standard_PCC(N,M,design_Ebn0);

PCC = PCC_structs(K_PCC);

% Ebn0_arr = 1.0:0.2:2.4;
Ebn0_arr = [2.25, 2.5, 2.75, 3.0];           % All under this Eb/n0.
min_errors = 1600;
N_ebn0 = length(Ebn0_arr);
L = 16;

BLERs = zeros(1, N_ebn0);


%% Simulation Loop.
for ebn0_iter = 1:N_ebn0
    Ebn0 = Ebn0_arr(ebn0_iter);
    BLERs(ebn0_iter) = sim_PCC(PCC, Ebn0, min_errors, L);
end

save('data/pccscl_ext.mat');
disp('Save file(s) complete!');
