%% Scaled Largest Eigenvalue (SLE) Detector
% Reference: https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=5403561
% X: n*p complex matrix, as the receiving data matrix

function  statistics = detcSleMul(X)

n = size(X, 1);
S = 1/n * X' * X;

% extract eigenvalues of SCM
[~, D] = eig(S);
d = real(diag(D));

% compute statistics
statistics = max(d) / mean(d);
end
