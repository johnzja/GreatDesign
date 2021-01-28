load data/first_error_pattern.mat;

%% Step1: Map all the statistics and PCC equations to Polarized Channels.
first_error_pattern = [first_error_pattern, zeros(1,K_CRC)];
PCC_conf = PCC_CRC_polar_config.PCC_conf;
% find parity bits.
K_PCC = PCC_conf.parity_bits_cnt;
unfrozen_bits_cnt = PCC_conf.unfrozen_bits_cnt;
parity_bits_index = PCC_conf.parity_bits_index;
nonfrozen_bits_logical = PCC_conf.nonfrozen_bits_logical;
parity_bits_loc = zeros(K_PCC, 1);

parity_checked = false(1,unfrozen_bits_cnt);

for k = 1:K_PCC
   parity_bits_loc(k) =  parity_bits_index{k}(end);
   parity_checked(parity_bits_index{k}) = true;
end

info_bits_wrt_nonfrozen = true(unfrozen_bits_cnt,1);
info_bits_wrt_nonfrozen(parity_bits_loc) = false;

statistics_PCC = zeros(1, unfrozen_bits_cnt);   % 268.
statistics_PCC(info_bits_wrt_nonfrozen) = first_error_pattern;
statistics_Polar = zeros(1,N);
statistics_Polar(nonfrozen_bits_logical) = statistics_PCC;

% Map the check equations to Polarized Channels.
parity_checked_Polar = false(1,N);
parity_checked_Polar(nonfrozen_bits_logical) = parity_checked;

% the Critical Set is previously prepared.

%% Analyze.
percentage_error_in_PCCEqns = sum(statistics_Polar.*parity_checked_Polar)/sum(statistics_Polar)
percentage_error_in_CS = sum(statistics_Polar .* (CS.'))/sum(statistics_Polar)


