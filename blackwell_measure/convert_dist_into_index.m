function [idx0, idx1, idx2, idx3] = convert_dist_into_index(s, Nbins)

idx0 = 1+round((Nbins-1)*s(1));
idx1 = 1+round((Nbins-1)*s(2));
idx2 = 1+round((Nbins-1)*s(3));
idx3 = 1+round((Nbins-1)*s(4));

% HERE: The 3-D quantization method has some bugs.
if idx0 + idx1 + idx2 + idx3 > 3 + Nbins
    d = idx0 + idx1 + idx2 + idx3 - 3 - Nbins;
    [~, M_index] = max([idx0, idx1, idx2, idx3]);
    switch M_index
        case 1
            idx0 = idx0 - d;
        case 2
            idx1 = idx1 - d;
        case 3
            idx2 = idx2 - d;
        case 4
            idx3 = idx3 - d;
    end
end

end

