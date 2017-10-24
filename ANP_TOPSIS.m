%initialize
n = 3;
m = 3;

%ANP
A = csvread('factor_rating.csv', 0, 0, [0, 0, n - 1, n - 1]);
[V, D] = eig(A);
[~, column] = max(diag(D));
v = V(:,column);
w_2 = v/sum(v);
B = csvread('factor_relationship.csv', 0, 0, [0, 0, n - 1, n - 1]);
w_c = B * w_2;
%w_c is the final weight distribution

%TOPSIS
D = csvread('alternatives.csv', 0, 0, [0, 0, m - 1, n - 1]);
pn = csvread('factor_attribution.csv', 0, 0, [0, 0, 0, n - 1]);
tmp = sum(D.*D);
D_n = D;
for i = 1 : n
    if pn(1, i) == 0
        D_n(:, i) = -D_n(:, i) / sqrt(tmp(1, i));
    else
        D_n(:, i) = D_n(:, i) / sqrt(tmp(1, i));
    end
end
PIS = max(D_n);
NIS = min(D_n);

D_ = zeros(m, 2);
C_ = zeros(m, 1);
for i = 1 : m
    dis_p = 0;
    dis_n = 0;
    for j = 1 : n
        dis_p = dis_p + w_c(j,1) * (D_n(i, j) - PIS(j))^2;
        dis_n = dis_n + w_c(j,1) * (D_n(i, j) - NIS(j))^2;
    end
    D_(i, 1) = sqrt(dis_p);
    D_(i, 2) = sqrt(dis_n);
    C_(i, 1) = D_(i, 2) / (D_(i, 1) + D_(i, 2));
end
[~, optimize_choice] = max(C_);
clear tmp i j column V dis_n dis_p v;