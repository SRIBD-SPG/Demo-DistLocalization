%% Energy Detection (ED)
% Reference: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1447503
% X: n*p complex matrix, as the receiving data matrix

function  statistics = detcEnergyMul(X)
    statistics = norm(X, 'fro')^2;
end