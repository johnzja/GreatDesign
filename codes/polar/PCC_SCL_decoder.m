function polar_info_esti = PCC_SCL_decoder(llr, L, PCC_conf, lambda_offset, llr_layer_vec, bit_layer_vec)
    % Parameters
    % llr: vector of length N.
    % L: list size.
    % M: cnt. of  bits.
    % H_crc: CRC-check matrix.
    % lambda_offset
    % llr_layer_vec
    % bit_layer_vec

    N = length(llr);
    n = log2(N);
    %% Extract info from PCC_conf structure.
    M = PCC_conf.info_bits_cnt;
    K = PCC_conf.parity_bits_cnt;
    PCEqns = PCC_conf.parity_bits_index;
    nonfrozen_bits_logical = PCC_conf.nonfrozen_bits_logical;
    parity_locs_wrt_nonfrozen_bits = zeros(1,K);
    info_bits_wrt_nonfrozen_logical = true(M+K,1);
    
    
    for k = 1:K
        t = PCEqns{k}(end);
        parity_locs_wrt_nonfrozen_bits(k) = t;
        info_bits_wrt_nonfrozen_logical(t) = false;
    end
    
    %% Memory declarations.
    lazy_copy = zeros(n, L);%If Lazy Copy is used, there is no data-copy in the decoding process. We only need to record where the data come from. Here,data refer to LLRs and partial sums.
    %Lazy Copy is a relatively sophisticated operation for new learners of polar codes. If you do not understand such operation, you can directly copy data.
    %If you can understand lazy copy and you just start learning polar codes
    %for just fews days, you are very clever,
    P = zeros(N - 1, L);        %Channel llr is public-used, so N - 1 is enough.
    C = zeros(N - 1, 2 * L);    %I do not esitimate (x1, x2, ... , xN), so N - 1 is enough.
    u = zeros(M+K, L);          %unfrozen bits that polar codes carry, including crc bits.
    PM = zeros(L, 1);           %Path metrics
    activepath = zeros(L, 1);   %Indicate if the path is active. '1'->active; '0' otherwise.
    cnt_u = 1;                  %information bit counter, from 1 to M+K.
    
    %% Initialize
    activepath(1) = 1;
    lazy_copy(:, 1) = 1;
    %decoding starts
    %default: in the case of path clone, the origianl path always corresponds to bit 0, while the new path bit 1.

    for phi = 0 : N - 1
        layer = llr_layer_vec(phi + 1);
        phi_mod_2 = mod(phi, 2);
        for l_index = 1 : L
            if activepath(l_index) == 0
                continue;
            end
            switch phi%Decoding bits u_0 and u_N/2 needs channel LLR, so the decoding of them is separated from other bits. 
                case 0
                    index_1 = lambda_offset(n);
                    for beta = 0 : index_1 - 1
                        P(beta + index_1, l_index) = sign(llr(beta + 1)) * sign(llr(beta + index_1 + 1)) * min(abs(llr(beta + 1)), abs(llr(beta + index_1 + 1)));
                    end
                    for i_layer = n - 2 : -1 : 0
                        index_1 = lambda_offset(i_layer + 1);
                        index_2 = lambda_offset(i_layer + 2);
                        for beta = 0 : index_1 - 1
                            P(beta + index_1, l_index) = sign(P(beta + index_2, l_index)) *...
                                sign(P(beta + index_1 + index_2, l_index)) * min(abs(P(beta + index_2, l_index)), abs(P(beta + index_1 + index_2, l_index)));
                        end
                    end
                case N/2
                    index_1 = lambda_offset(n);
                    for beta = 0 : index_1 - 1
                        x_tmp = C(beta + index_1, 2 * l_index - 1);
                        P(beta + index_1, l_index) = (1 - 2 * x_tmp) * llr(beta + 1) + llr(beta + 1 + index_1);
                    end
                    for i_layer = n - 2 : -1 : 0
                        index_1 = lambda_offset(i_layer + 1);
                        index_2 = lambda_offset(i_layer + 2);
                        for beta = 0 : index_1 - 1
                            P(beta + index_1, l_index) = sign(P(beta + index_2, l_index)) *...
                                sign(P(beta + index_1 + index_2, l_index)) * min(abs(P(beta + index_2, l_index)), abs(P(beta + index_1 + index_2, l_index)));
                        end
                    end
                otherwise
                    index_1 = lambda_offset(layer + 1);
                    index_2 = lambda_offset(layer + 2);
                    for beta = 0 : index_1 - 1
                        P(beta + index_1, l_index) = (1 - 2 * C(beta + index_1, 2 * l_index - 1)) * P(beta + index_2, lazy_copy(layer + 2, l_index)) +...
                            P(beta + index_1 + index_2, lazy_copy(layer + 2, l_index));
                    end
                    for i_layer = layer - 1 : -1 : 0
                        index_1 = lambda_offset(i_layer + 1);
                        index_2 = lambda_offset(i_layer + 2);
                        for beta = 0 : index_1 - 1
                            P(beta + index_1, l_index) = sign(P(beta + index_2, l_index)) *...
                                sign(P(beta + index_1 + index_2, l_index)) * min(abs(P(beta + index_2, l_index)),...
                                abs(P(beta + index_1 + index_2, l_index)));
                        end
                    end
            end
        end

        % HERE: down to the leaf of decoding-tree.
        pcc_eqn_index = find(parity_locs_wrt_nonfrozen_bits==cnt_u,1);

        if nonfrozen_bits_logical(phi + 1) && isempty(pcc_eqn_index)    
            % if now we decode an info bit
            PM_pair = realmax * ones(2, L);
            for l_index = 1 : L
                if activepath(l_index) == 0
                    continue;
                end

                % It is an ordinary info bit.
                if P(1, l_index) >= 0           % if this bit looks like 0:
                    PM_pair(1, l_index) = PM(l_index);                  % decision: 0
                    PM_pair(2, l_index) = PM(l_index) + P(1, l_index);  % decision: 1
                else
                    PM_pair(1, l_index) = PM(l_index) - P(1, l_index);  % decision: 0
                    PM_pair(2, l_index) = PM(l_index);                  % decision: 1
                end
            end

            % Perform path splitting.
            middle = min(2 * sum(activepath), L);
            PM_sort = sort(PM_pair(:));
            PM_cv = PM_sort(middle);
            compare = PM_pair <= PM_cv; 
            kill_index = zeros(L, 1);   % to record the index of the paths that are killed
            kill_cnt = 0;               % the total number of killed path
            % the above two variables consist of a stack
            for i = 1 : L
                if (compare(1, i) == 0)&&(compare(2, i) == 0)%which indicates that this path should be killed
                    activepath(i) = 0;
                    kill_cnt = kill_cnt + 1;    %push stack
                    kill_index(kill_cnt) = i;
                end
            end
            for l_index = 1 : L
                if activepath(l_index) == 0
                    continue;
                end
                path_state = compare(1, l_index) * 2 + compare(2, l_index);
                switch path_state%path_state can equal to 0, but in this case we do no operation.
                    case 1
                        u(cnt_u, l_index) = 1;                  % the second row of "compare" wins. 
                        C(1, 2 * l_index - 1 + phi_mod_2) = 1;  % the new decision is bit 1.
                        PM(l_index) = PM_pair(2, l_index);
                    case 2
                        u(cnt_u, l_index) = 0;                  % the new decision is bit 0.
                        C(1, 2 * l_index - 1 + phi_mod_2) = 0;
                        PM(l_index) = PM_pair(1, l_index);
                    case 3                                      % All of 0,1 are both promising.
                        index = kill_index(kill_cnt);
                        kill_cnt = kill_cnt - 1;                % pop stack
                        activepath(index) = 1;
                        % lazy copy.
                        lazy_copy(:, index) = lazy_copy(:, l_index);    % index: son; l_index:father.
                        u(:, index) = u(:, l_index);
                        u(cnt_u, l_index) = 0;
                        u(cnt_u, index) = 1;
                        C(1, 2 * l_index - 1 + phi_mod_2) = 0;
                        C(1, 2 * index - 1 + phi_mod_2) = 1;
                        PM(l_index) = PM_pair(1, l_index);
                        PM(index) = PM_pair(2, l_index);
                end
            end
            cnt_u = cnt_u + 1;
        else    % non-info bit operation.
            
            if nonfrozen_bits_logical(phi+1) == 1
                for l_index = 1:L
                    if activepath(l_index) == 0
                        continue;
                    end
                    % Calculate parity check using previously-decoded bits.
                    pcc_eqn = PCEqns{pcc_eqn_index}(1:end-1);
                    p = mod(sum(u(pcc_eqn,l_index)),2);

                    % LLR<0: like 1; LLR>0: like 0.
                    if ((1-2*p)*P(1, l_index)) < 0    % this bit looks different from p, then add loss to PM.
                        PM(l_index) = PM(l_index) + abs(P(1, l_index));
                    end

                    if phi_mod_2 == 0               % odd-indexed frozen bit.
                        C(1, 2 * l_index - 1) = p;  % Left.
                    else
                        C(1, 2 * l_index) = p;      % Right.
                    end 
                    u(cnt_u, l_index) = p;          
                end
                cnt_u = cnt_u + 1;
            else
                for l_index = 1:L
                    if activepath(l_index) == 0
                        continue;
                    end
                    if P(1, l_index) < 0    % this bit looks like 1, then add loss to PM.
                        PM(l_index) = PM(l_index) - P(1, l_index);
                    end
                    if phi_mod_2 == 0               % odd-indexed frozen bit.
                        C(1, 2 * l_index - 1) = 0;  % Left.
                    else
                        C(1, 2 * l_index) = 0;      % Right.
                    end 
                end
            end

        end 

        for l_index = 1 : L%partial-sum return
            if activepath(l_index) == 0
                continue
            end
            if (phi_mod_2  == 1) && (phi ~= N - 1)
                layer = bit_layer_vec(phi + 1);
                for i_layer = 0 : layer - 1
                    index_1 = lambda_offset(i_layer + 1);
                    index_2 = lambda_offset(i_layer + 2);
                    for beta = index_1 : 2 * index_1 - 1
                        C(beta + index_1, 2 * l_index) = mod(C(beta, 2 *  lazy_copy(i_layer + 1, l_index) - 1) + C(beta, 2 * l_index), 2);%Left Column lazy copy
                        C(beta + index_2, 2 * l_index) = C(beta, 2 * l_index);   
                    end
                end
                index_1 = lambda_offset(layer + 1);
                index_2 = lambda_offset(layer + 2);
                for beta = index_1 : 2 * index_1 - 1
                    C(beta + index_1, 2 * l_index - 1) = mod(C(beta, 2 * lazy_copy(layer + 1, l_index) - 1) + C(beta, 2 * l_index), 2);%Left Column lazy copy
                    C(beta + index_2, 2 * l_index - 1) = C(beta, 2 * l_index);
                end 
            end
        end
        %lazy copy
        if phi < N - 1
            for i_layer = 1 : llr_layer_vec(phi + 2) + 1
                for l_index = 1 : L
                    lazy_copy(i_layer, l_index) = l_index;
                end
            end
        end
    end

    %% All the source bits are decoded now.
    %path selection.
    [~, ordered_path_index] = sort(PM);
    polar_info_esti = u(:,ordered_path_index(1));
    % Remove parity bits.
    polar_info_esti = polar_info_esti(info_bits_wrt_nonfrozen_logical);
    polar_info_esti = logical(polar_info_esti.');
end