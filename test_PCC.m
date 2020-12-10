clear;
load data/PCC_config.mat;

N = 512;
M = 256;

N_sim = 1000;

for sim_iter = 1:N_sim
    random_bits = (rand([1,M])>0.5);
    encoded_bits = PCC_polar_encoder(random_bits, PCC);
    
end


% SC decoding.
nonfrozen_bits_logical = PCC.nonfrozen_bits_logical;




