%% factor-analysis decomposition (FAD) via MLE method
%
%  maximize -logdet(R) - tr(inv(R) * S)
%      s.t. R = H*H' + Sgm
%
%  ref: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8911457 (Sec. II-A) 
%  author: Rui Zhou
%  S: the sample covariance matrix
%  r: the rank of low-rank part (number of factors)
%  noise_power_range: range of noise power, if set negative then no range restrict
%  max_iter/ptol/ftol: algorithm control options

function [H, Sgm, objs] = FAD(S, r, noise_power_range, max_iter, ptol, ftol)

% initialize Sgm
[V, D] = eigDesOrder(S, r);
Sgm = diag(diag(S - V*D*V'));
R = V*D*V' + Sgm;

% create space
objs = [];

for iter = 1:max_iter
%     record current status
    R_old = R;
    
%     update H
    [V, D] = eigDesOrder((Sgm^-0.5) * S * (Sgm^-0.5), r);
    H = (Sgm^0.5) * V * diag(max(diag(D)-1, 0))^0.5;
    
%     update Sgm
    Sgm = diag(S - H*H');
    if (all(noise_power_range >= 0))
        Sgm = min(Sgm, noise_power_range(:, 2));
        Sgm = max(Sgm, noise_power_range(:, 1));
    end
    Sgm = diag(Sgm);
    
%     assemble R and recompute objective
    R = H*H' + Sgm;
    objs = [objs; -log(det(R)) - trace(inv(R) * S)];
    
%     check convergence
    if (all(abs(R - R_old) <= ptol * abs(R), 'all') && abs(objs(end)-objs(end-1)) < ftol * abs(objs(end)))
%       if (all(abs(R) -abs( R_old) <= ptol * abs(R), 'all') && abs(objs(end))-abs(objs(end-1)) < ftol * abs(objs(end)))
        break;
    end
end
end

