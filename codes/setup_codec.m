channel_codec.type = 'polar';    % 'conv', 'polar', 'RS', ...

if strcmp(channel_codec.type, 'polar')
    N = 512;
    K = 256;
    
    polar_codec.N = 512;
    polar_codec.K = 256;
    polar_codec.info_bits_logical = logical(1,N);
    polar_codec.concatenate = 'CRC';
    polar_codec.concatenate_conf = '18';
elseif strcmp(channel_codec.type, 'conv')
    
    
end