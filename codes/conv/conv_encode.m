function encoded_bits = conv_encode(input_bits, conv_encoder_conf)
%
    %% Step1: Determine the number of output bits.
    n = conv_encoder_conf.n;
    k = conv_encoder_conf.k;
    N = conv_encoder_conf.N;
    A = conv_encoder_conf.A;
    
    N_inner_state_bits = (N-1)*k;
    % Align input bits.
    n_input_bits = length(input_bits);
    temp = mod(n_input_bits, k);
    if temp
        n_input_bits = n_input_bits + k - temp;
        input_bits = [input_bits, false([1, k-temp])];  % align
    end
    
    % calculate output length.
    n_input_blocks = n_input_bits/k;
    if conv_encoder_conf.trailing
        n_output_blocks = n_input_blocks+N-1;
    else
        n_output_blocks = n_input_blocks;
    end
    n_output_bits = n_output_blocks*n;
    encoded_bits = false([1, n_output_bits]);
    
    % initialize registers. From left to right.
    registers = cell(1,N-1);
    for iter=1:N-1
        registers{iter} = false([1, k]);
    end
    
    %% Step2: Start convolution.
    % zero-padding.
    if conv_encoder_conf.trailing
        n_input_blocks = n_output_blocks;
        input_bits = [input_bits, false([1, N_inner_state_bits])];
    end
    
    for iter=1:n_input_blocks
        i_index = k*(iter-1)+1;
        i_block = flip(input_bits(i_index:i_index+k-1));
        % calculate convolution.
        i_vct = combine(i_block, registers);
        o_block = false([1,n]);
        for o_iter=1:n
            o_block(o_iter) = logical(mod(i_vct*(A{o_iter}.'),2));
        end
        % append encoder output to output sequence.
        o_index = n*(iter-1)+1;
        encoded_bits(o_index:o_index+n-1)=o_block;
        % shift register.
        for r_iter=N-1:-1:2
            registers{r_iter} = registers{r_iter-1};
        end
        registers{1} = i_block;
    end
    
end

function ret=combine(i_block, registers)
    L = length(registers);
    Li = length(i_block);
    Lret = (L+1)*Li;
    ret = false([1, Lret]);
    ret(1:Li) = i_block;
    for k=1:L
        index=Li*(k)+1;
        ret(index:index+Li-1) = registers{k};
    end
end