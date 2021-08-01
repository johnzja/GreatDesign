function I = get_I_4D(bm_dist, bin_centers)
    [~, ~, Nbins] = size(bm_dist);
    assert(Nbins == length(bin_centers));
    
    I = 0;
    p_total = 0;
    for idx0 = 1:Nbins
        for idx1 = 1:(Nbins+1-idx0)
            p_total_layer1 = 0;
            I_layer1 = 0;
            
            for idx2 = 1:(Nbins+2-idx0-idx1)
                s0 = bin_centers(idx0);
                s1 = bin_centers(idx1);
                s2 = bin_centers(idx2);
                s3 = 1-s0-s1-s2;
                
                Ip = 2;
                if s0 > 0
                    Ip = Ip + s0 * log2(s0);
                end
                if s1 > 0
                    Ip = Ip + s1 * log2(s1);
                end
                if s2 > 0
                    Ip = Ip + s2 * log2(s2);
                end
                if s3 > 0
                    Ip = Ip + s3 * log2(s3);
                end
                
                p_total_layer1 = p_total_layer1 + bm_dist(idx0, idx1, idx2);
                I_layer1 = I_layer1 + Ip * bm_dist(idx0, idx1, idx2);
            end
            
            p_total = p_total + p_total_layer1;
            I = I + I_layer1;
        end
    end

    fprintf('Check BM: Total mass = %f\n', p_total);
end