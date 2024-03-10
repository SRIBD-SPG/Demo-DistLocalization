%% Maximum-Minimum Eigenvalue (MME) Detector
% Reference: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=5089517
% X: n*p complex matrix, as the receiving data matrix

function  statistics = detcMmeMul(X)

n = size(X, 1);
S = 1/n * X' * X;

% extract eigenvalues of SCM
[~, D] = eig(S);
d = real(diag(D));

% compute statistics
statistics = max(d) / min(d);
end
