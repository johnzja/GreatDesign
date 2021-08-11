function Pe = get_pe_4D(bm_dist, bin_centers)
    [~, ~, Nbins] = size(bm_dist);
    assert(Nbins == length(bin_centers));
    
    Pe = 0;
    p_total = 0;
    for idx0 = 1:Nbins
        for idx1 = 1:(Nbins+1-idx0)
            p_total_layer1 = 0;
            Pe_layer1 = 0;
            
            for idx2 = 1:(Nbins+2-idx0-idx1)
                s0 = bin_centers(idx0);
                s1 = bin_centers(idx1);
                s2 = bin_centers(idx2);
                s3 = 1-s0-s1-s2;
                
                [sm, ~] = max([s0, s1, s2, s3]);
                
                p_total_layer1 = p_total_layer1 + bm_dist(idx0, idx1, idx2);
                Pe_layer1 = Pe_layer1 +  (1-sm) * bm_dist(idx0, idx1, idx2);
            end
            
            p_total = p_total + p_total_layer1;
            Pe = Pe + Pe_layer1;
        end
    end
    
    eps = 1e-7;
    if abs(1-p_total) > eps
        fprintf('Check BM: Total mass = %f\n', p_total);
    end
end