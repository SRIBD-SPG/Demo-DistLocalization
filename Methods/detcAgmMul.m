%% Arithmetic to Geometric Mean (AGM) Detector
% Reference: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=5397901
% X: n*p complex matrix, as the receiving data matrix

function  statistics = detcAgmMul(X)

n = size(X, 1);
S = 1/n * X' * X;

% extract eigenvalues of SCM
[~, D] = eig(S);
d = real(diag(D));

% compute statistics
statistics = mean(d) / geo_mean(d);
end
