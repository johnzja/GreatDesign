function [idx0, idx1, idx2, idx3] = convert_dist_into_index(s, Nbins)

ret = zeros(1,4);
idx0 = 1+round((Nbins-1)*s(1));
idx1 = 1+round((Nbins-1)*s(2));
idx2 = 1+round((Nbins-1)*s(3));
idx3 = 1+round((Nbins-1)*s(4));

ret = [idx0, idx1, idx2, idx3];
% 1 <= idx0 <= Nbins is guaranteed.

% [~, max_index] = max(s);

d = sum(ret) - 3 - Nbins;  

% Quantization error, which may cause deviation from the equation:
% (idx0-1) + ... + (idx3-1) = Nbins-1.
if d == 0
    return
elseif d == 1
    % Find the number closest to Nbins/2.
    [~, idx_max] = max(ret);
    ret(idx_max) = ret(idx_max)-1;
elseif d == -1
    [~, idx_min] = min(ret);
    ret(idx_min) = ret(idx_min)+1;
else
    error('Unexpected.');
end

% switch max_index
%     case 1
%         idx0 = 3+Nbins-idx1-idx2-idx3;
%     case 2
%         idx1 = 3+Nbins-idx0-idx2-idx3;
%     case 3
%         idx2 = 3+Nbins-idx0-idx1-idx3;
%     case 4
%         idx3 = 3+Nbins-idx0-idx1-idx2; 
% end


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
idx0 = ret(1);
idx1 = ret(2);
idx2 = ret(3);
idx3 = ret(4);
end

