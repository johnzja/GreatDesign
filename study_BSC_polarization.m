syms a1 a2 a3 a4;
syms p1 p2 p3;
syms p;


W = [a1*[1-p,p;p,1-p],a2*[(1-p)^2,p^2;p^2,(1-p)^2],a3*[(1-p)^3,p^3;p^3,(1-p)^3], a4*[(1-p)^4,p^4;p^4,(1-p)^4]];

[~,N] = size(W);
% Wp = zeros(2,2*N^2);
for u2 = 0:1
    for u1 = 0:1
        for s1 = 1:N
            for s2 = 1:N
                Wp(u2+1,u1*(N^2) + N*(s2-1) + s1) = 1/2*W(mod(u1+u2,2)+1,s1)*W(u2+1,s2);
            end
        end
    end
end

% calculate pseudo-likelihood ratios.
for k = 1:72
    LR(k) = Wp(1,k)/Wp(2,k);
end
LRu = unique(LR);