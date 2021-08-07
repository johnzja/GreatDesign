function ret_dist = down_transform_sym(dist_1, dist_2, bin_centers, GF_info)
    [~, ~, Nbins] = size(dist_1);
    assert(all(size(dist_1) == size(dist_2)));
    assert(Nbins == length(bin_centers));
    bm_dist_1_visited = false(Nbins,Nbins,Nbins);
    
    ret_dist = zeros(Nbins, Nbins, Nbins);
    
    for idx1_0 = 1:Nbins
        for idx1_1 = 1:(Nbins+1-idx1_0)
            for idx1_2 = 1:(Nbins+2-idx1_0-idx1_1)
                idx1_3 = Nbins+3-idx1_0-idx1_1-idx1_2;
                
                
                % Access blackwell-measure distribution 1.
                b0 = bm_dist_1_visited(idx1_0, idx1_1, idx1_2);
                b1 = bm_dist_1_visited(idx1_1, idx1_0, idx1_3);
                b2 = bm_dist_1_visited(idx1_2, idx1_3, idx1_0);
                b3 = bm_dist_1_visited(idx1_3, idx1_2, idx1_1);
                if any([b0,b1,b2,b3])
                    % assert(all([b0,b1,b2,b3]));
                    continue;
                end

                % Fetch the distribution pdf.
                t0 = dist_1(idx1_0, idx1_1, idx1_2);
                t1 = dist_1(idx1_1, idx1_0, idx1_3);
                t2 = dist_1(idx1_2, idx1_3, idx1_0);
                t3 = dist_1(idx1_3, idx1_2, idx1_1);

                bm_dist_1_visited(idx1_0, idx1_1, idx1_2)=true;
                bm_dist_1_visited(idx1_1, idx1_0, idx1_3)=true;
                bm_dist_1_visited(idx1_2, idx1_3, idx1_0)=true;
                bm_dist_1_visited(idx1_3, idx1_2, idx1_1)=true;

                p1 = t0+t1+t2+t3;
                if p1 == 0
                    continue;
                end

                p = [bin_centers(idx1_0), bin_centers(idx1_1),bin_centers(idx1_2),bin_centers(idx1_3)].';
                pc2 = [p(1), p(3), p(4), p(2)].';
                
                bm_dist_2_visited = false(Nbins, Nbins,Nbins);
                
                % Go on to fetch blackwell-measure distribution 2.
                for idx2_0 = 1:Nbins
                    for idx2_1 = 1:(Nbins+1-idx2_0)
                        for idx2_2 = 1:(Nbins+2-idx2_0-idx2_1)
                            idx2_3 = Nbins+3-idx2_0-idx2_1-idx2_2;
                            
                            b0 = bm_dist_2_visited(idx2_0, idx2_1, idx2_2);
                            b1 = bm_dist_2_visited(idx2_1, idx2_0, idx2_3);
                            b2 = bm_dist_2_visited(idx2_2, idx2_3, idx2_0);
                            b3 = bm_dist_2_visited(idx2_3, idx2_2, idx2_1);
                            if any([b0,b1,b2,b3])
                                % assert(all([b0,b1,b2,b3]));
                                continue;
                            end
                            
                            % Fetch the distribution pdf.
                            t0 = dist_2(idx2_0, idx2_1, idx2_2);
                            t1 = dist_2(idx2_1, idx2_0, idx2_3);
                            t2 = dist_2(idx2_2, idx2_3, idx2_0);
                            t3 = dist_2(idx2_3, idx2_2, idx2_1);

                            bm_dist_2_visited(idx2_0, idx2_1, idx2_2)=true;
                            bm_dist_2_visited(idx2_1, idx2_0, idx2_3)=true;
                            bm_dist_2_visited(idx2_2, idx2_3, idx2_0)=true;
                            bm_dist_2_visited(idx2_3, idx2_2, idx2_1)=true;
                            
                            p2 = t0+t1+t2+t3;
                            if p2 == 0
                                continue;
                            end
                            
                            q = [bin_centers(idx2_0), bin_centers(idx2_1),bin_centers(idx2_2),bin_centers(idx2_3)].';
                            
                            % Perform down-polarization.
                            qab     = [q(2), q(1), q(4), q(3)].';
                            qab3    = [q(3), q(4), q(1), q(2)].';
                            qb2     = [q(4), q(3), q(2), q(1)].';
                            
                            v = zeros(4,4);
                            v(:,1) = pc2 .* q;
                            v(:,2) = pc2 .* qab;
                            v(:,3) = pc2 .* qab3;
                            v(:,4) = pc2 .* qb2;
                            
                            lambdas = sum(v);
                            ch_prob = p1 * p2/4;
                            
                            for idx = 1:4
                                if lambdas(idx) > 0
                                    pvec = v(:,idx) / lambdas(idx);
                                    [idx_Wp_0, idx_Wp_1, idx_Wp_2, idx_Wp_3] = convert_dist_into_index(pvec, Nbins);
                                    ret_dist(idx_Wp_0, idx_Wp_1, idx_Wp_2) = ret_dist(idx_Wp_0, idx_Wp_1, idx_Wp_2) + ch_prob*lambdas(idx);
                                    ret_dist(idx_Wp_1, idx_Wp_0, idx_Wp_3) = ret_dist(idx_Wp_1, idx_Wp_0, idx_Wp_3) + ch_prob*lambdas(idx);
                                    ret_dist(idx_Wp_2, idx_Wp_3, idx_Wp_0) = ret_dist(idx_Wp_2, idx_Wp_3, idx_Wp_0) + ch_prob*lambdas(idx);
                                    ret_dist(idx_Wp_3, idx_Wp_2, idx_Wp_1) = ret_dist(idx_Wp_3, idx_Wp_2, idx_Wp_1) + ch_prob*lambdas(idx);
                                    
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
