function [idx0, idx1, idx2, idx3] = convert_dist_into_index(s, Nbins)

idx0 = 1+round((Nbins-1)*s(1));
idx1 = 1+round((Nbins-1)*s(2));
idx2 = 1+round((Nbins-1)*s(3));
idx3 = 1+round((Nbins-1)*s(4));
[~, max_index] = max(s);

switch max_index
    case 1
        idx0 = 3+Nbins-idx1-idx2-idx3;
    case 2
        idx1 = 3+Nbins-idx0-idx2-idx3;
    case 3
        idx2 = 3+Nbins-idx0-idx1-idx3;
    case 4
        idx3 = 3+Nbins-idx0-idx1-idx2; 
end


% HERE: The 3-D quantization method has some bugs.
% Four equivalent forms.
% (0, 1, 2, 3)
% (1, 0, 3, 2)
% (2, 3, 0, 1)
% (3, 2, 1, 0).


% if idx0 + idx1 + idx2 + idx3 > 3 + Nbins
%     [~, max_index] = max([idx0, idx1, idx2, idx3]);
%     
% end

end

