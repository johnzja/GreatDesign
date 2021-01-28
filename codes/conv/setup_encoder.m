conv_encoder_conf.n = 2;    % Output bits.
conv_encoder_conf.k = 1;    % Input bits.
conv_encoder_conf.N = 4;    
conv_encoder_conf.window_factor=6;
conv_encoder_conf.trailing = true;

A = cell(conv_encoder_conf.n,1);
%       a4,a3,a2,a1
 
A{1} = [1, 1, 0, 1];  % 15
A{2} = [1, 1, 1, 1];  % 17. Rate=1/2.

conv_encoder_conf.A=A;
conv_encoder_conf.loss_func = @hamming_distance;

function d = hamming_distance(x, y)
    assert((length(x)==length(y)));
    d=sum(xor(x,y));
end