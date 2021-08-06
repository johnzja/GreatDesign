%%

p = [0.7, 0.3, 0, 0].';
q = [0.8, 0.2, 0, 0].';

%% Calculate W-.
qc      = [q(1),q(4),q(2),q(3)].';
qcab    = [q(4),q(1),q(3),q(2)].';
qcab3   = [q(2),q(3),q(1),q(4)].';
qcb2    = [q(3),q(2),q(4),q(1)].';

pWn = [p.'*qc, p.'*qcab, p.'*qcab3, p.'*qcb2].';
IWn = get_QSC_Capacity(pWn);
IW1 = get_QSC_Capacity(p);
IW2 = get_QSC_Capacity(q);

pc2     = [p(1), p(3), p(4), p(2)].';
qab     = [q(2), q(1), q(4), q(3)].';
qab3    = [q(3), q(4), q(1), q(2)].';
qb2     = [q(4), q(3), q(2), q(1)].';

v0 = pc2 .* q;
v1 = pc2 .* qab;
v2 = pc2 .* qab3;
v3 = pc2 .* qb2;

lambdas = zeros(4,1);
lambdas(1) = sum(v0);
lambdas(2) = sum(v1);
lambdas(3) = sum(v2);
lambdas(4) = sum(v3);

IWps = zeros(4,1);
IWps(1) = get_QSC_Capacity(v0/lambdas(1));
IWps(2) = get_QSC_Capacity(v1/lambdas(2));
IWps(3) = get_QSC_Capacity(v2/lambdas(3));
IWps(4) = get_QSC_Capacity(v3/lambdas(4));

IWp = lambdas.'*IWps;

%% Mutual Information.
function I = get_QSC_Capacity(p)
    I = 2;
    for idx = 1:4
        s = p(idx);
        if s>0
            I = I + s*log2(s);
        end
    end
end