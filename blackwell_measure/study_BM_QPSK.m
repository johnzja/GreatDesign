%% Study 3-D Blackwell Measure.
P = 1;
A0 = sqrt(P)*[1,1];
A1 = sqrt(P)*[-1,1];
A2 = sqrt(P)*[-1,-1];
A3 = sqrt(P)*[1, -1];

sigma = 0.8;

% Do QPSK sampling.
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
    
    s0 = s0 / sum([s0, s1, s2, s3]);
    s1 = s1 / sum([s0, s1, s2, s3]);
    s2 = s2 / sum([s0, s1, s2, s3]);
    s3 = s3 / sum([s0, s1, s2, s3]);
    
    bm_posterior_coord(idx, :) = [s1, s2, s3];
end

scatter3(bm_posterior_coord(:,1), bm_posterior_coord(:,2), bm_posterior_coord(:,3), '*');
xlabel('s1'); ylabel('s2'); zlabel('s3');
axis equal;
