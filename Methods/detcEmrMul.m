%% Eigenvalue-Moment-Ratio (EMR) Detector
% Reference: https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=6905846
% X: n*p complex matrix, as the receiving data matrix

function  statistics = detcEmrMul(X)

[n, p] = size(X);
S = 1/n * X' * X;

% extract eigenvalues of SCM
M1 = 1/p * real(trace(S));
M2 = 1/p * norm(S, 'fro')^2;

% compute statistics
statistics = M2 / (M1^2);
end
