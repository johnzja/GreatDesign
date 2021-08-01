function plot_BM(bm_dist, bin_centers)

[~, ~, Nbins] = size(bm_dist);
v0 = [1, 0, 0];
v1 = [0.5, sqrt(3)/2, 0];
v2 = [0.5, 0.5/sqrt(3), sqrt(6)/3];


N_3d_bins = 0;
for idx0 = 1:Nbins
    for idx1 = 1:(Nbins+1-idx0)
        for idx2 = 1:(Nbins+2-idx0-idx1)
            N_3d_bins = N_3d_bins + 1;
        end
    end
end

X = zeros(N_3d_bins, 1);
Y = zeros(N_3d_bins, 1);
Z = zeros(N_3d_bins, 1);
S = zeros(N_3d_bins, 1);
C = zeros(N_3d_bins, 1);

index = 1;
for idx0 = 1:Nbins
    for idx1 = 1:(Nbins+1-idx0)
        for idx2 = 1:(Nbins+2-idx0-idx1)
            s0 = bin_centers(idx0);
            s1 = bin_centers(idx1);
            s2 = bin_centers(idx2);
            s3 = 1-s0-s1-s2;
            
            V = s0 * v0 + s1 * v1 + s2 * v2;
            X(index) = V(1); Y(index) = V(2); Z(index) = V(3);
            S(index) = 1;
            C(index) = bm_dist(idx0, idx1, idx2);
            
            index = index + 1;
        end
    end
end

non_zero_logical = (C>0);
X = X(non_zero_logical);
Y = Y(non_zero_logical);
Z = Z(non_zero_logical);
S = S(non_zero_logical);
C = C(non_zero_logical);

figure;

cmap = hot(256);
cmap = cmap(256:-1:1,:);
colormap(cmap);
scatter3(X, Y, Z, S, C);
xlabel('x');
ylabel('y');
zlabel('z');

end