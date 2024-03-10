%% Generalized Likelihood Ratio (GLT) Detector
% Reference: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=5755211
% X: n*p complex matrix, as the receiving data matrix

function  statistics = GLR(X, r, noise_power_range)

n = size(X, 1);
S = 1/n * X' * X;

%MLE decomposition
if (nargin < 2) 
    r = 1;
end
max_iter = 1e2;
ptol = 1e-6;
ftol = 1e-6;
[H, Sgm, ~] = FAD(S, r, noise_power_range, max_iter, ptol, ftol);
% disp(objs);

% compute statistics
H1_M = H*H'+Sgm;
H0_N = diag(diag(S));
% aa=log(det(H0_N))+trace(H0_N\S);
% bb=log(det(H1_M))+trace(H1_M\S);
statistics = real((trace(H0_N\S)-trace(H1_M\S)+log(det(H0_N))-log(det(H1_M))));
end
