% down-transform is more complicated.
function ret_dist = down_transform_4D(dist_1, dist_2, bin_centers, GF_info)
    [~, ~, Nbins] = size(dist_1);
    assert(all(size(dist_1) == size(dist_2)));
    assert(Nbins == length(bin_centers));
    
    add_table = GF_info.add_table + 1;
    mult_table = GF_info.mult_table + 1;    % MATLAB representation.
    alpha = GF_info.alpha + 1;
    
    ret_dist = zeros(Nbins, Nbins, Nbins);
    for idx0_d1 = 1:Nbins
        for idx1_d1 = 1:(Nbins+1-idx0_d1)
            for idx2_d1 = 1:(Nbins+2-idx0_d1-idx1_d1)
                for idx0_d2 = 1:Nbins
                    for idx1_d2 = 1:(Nbins+1-idx0_d2)
                        for idx2_d2 = 1:(Nbins+2-idx0_d2-idx1_d2)
                            s0 = bin_centers(idx0_d1);
                            s1 = bin_centers(idx1_d1);
                            s2 = bin_centers(idx2_d1);
                            s3 = 1-s0-s1-s2;
                            s = [s0, s1, s2, s3];
                            
                            t0 = bin_centers(idx0_d2);
                            t1 = bin_centers(idx1_d2);
                            t2 = bin_centers(idx2_d2);
                            t3 = 1-t0-t1-t2;     
                            t = [t0, t1, t2, t3];
                            
                            for u1 = 1:4
                                r = [0, 0, 0, 0];
                                for u2 = 1:4
                                    % x1 = u1 + alpha*u2.
                                    x1 = add_table(idx_r, mult_table(u2, alpha));
                                    x2 = u2;
                                    r(u2) = s(x1) * t(x2);
                                end
                                
                                r = r / sum(r);
                                p1 = dist_1(idx0_d1, idx1_d1, idx2_d1);
                                p2 = dist_2(idx0_d2, idx1_d2, idx2_d2);

                                [idx0, idx1, idx2, ~] = convert_dist_into_index(r, Nbins);
                                ret_dist(idx0, idx1, idx2) = ret_dist(idx0, idx1, idx2) + (1/4) * p1 * p2;
                            end
                        end
                    end
                end
            end
        end
    end
    
end